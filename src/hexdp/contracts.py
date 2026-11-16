"""Plugin contracts: what a detector, strategy, op, discretiser, filler, exporter or interface must provide.

Every plugin is a plain object with a ``name``, a ``params`` schema and the methods of its kind. The
contracts are ``typing.Protocol`` classes, so a plugin needs no base class; ``validate`` checks an
object against its contract and is applied when the registry loads a plugin.

Two things make the plugins usable by an agent and a UI, not only by code:

* ``params`` is a schema (``ParamSpec``: default, unit, description, bounds), so a caller can list and
  explain the knobs without reading the source.
* ``applicable`` returns a ``Score`` with reasons in plain words. The pipeline and the agent write
  those reasons into the decision log (``Context.decide``), which is how a run explains itself.
"""

from __future__ import annotations

import json
import os
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Protocol, runtime_checkable

from hexdp.core.blockstructure import BlockStructure
from hexdp.core.geometry import Geometry


@dataclass(frozen=True)
class ParamSpec:
    default: Any
    description: str
    unit: str | None = None
    lower: float | None = None
    upper: float | None = None
    choices: tuple | None = None
    kind: str | None = None        # "path": a file or directory; must exist, made absolute by resolve_params


@dataclass
class Score:
    """How well a plugin fits. 0 = not applicable; reasons say why, ``missing`` what it would need."""

    value: float
    reasons: list[str] = field(default_factory=list)
    missing: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return self.value > 0

    @classmethod
    def no(cls, *reasons: str, missing: tuple[str, ...] = ()) -> "Score":
        return cls(0.0, list(reasons), list(missing))


@dataclass
class Feature:
    """A detector finding: what (kind), where (entities), how sure (confidence), and why (reasons).
    Kinds used so far: rotation_axis, sweep, thin_wall, blade, periodic_pair."""

    kind: str
    data: dict = field(default_factory=dict)
    entities: tuple[int, ...] = ()
    confidence: float = 1.0
    reasons: list[str] = field(default_factory=list)


@dataclass
class Decision:
    stage: str
    choice: str
    reasons: list[str]
    alternatives: dict[str, float] = field(default_factory=dict)
    time: float = field(default_factory=time.time)


@dataclass
class Context:
    """What a plugin may use while it runs: its own work directory, a core budget, a time limit and
    the decision log. Plugins write files only below ``workspace``."""

    workspace: Path
    cores: int = 1
    timeout_s: float | None = None
    cpu_list: str | None = None        # taskset CPU list for subprocesses, e.g. "0-7"
    decisions: list[Decision] = field(default_factory=list)

    def __post_init__(self):
        self.workspace = Path(self.workspace)
        self.workspace.mkdir(parents=True, exist_ok=True)

    def path(self, *parts: str) -> Path:
        p = self.workspace.joinpath(*parts)
        p.parent.mkdir(parents=True, exist_ok=True)
        return p

    def decide(self, stage: str, choice: str, reasons, alternatives: dict[str, float] | None = None) -> None:
        self.decisions.append(Decision(stage, choice, list(reasons), dict(alternatives or {})))

    def save_decisions(self, name: str = "decisions.json") -> Path:
        p = self.path(name)
        p.write_text(json.dumps([asdict(d) for d in self.decisions], indent=1))
        return p


# ------------------------------------------------------------------------------------------- contracts
@runtime_checkable
class Detector(Protocol):
    name: str
    params: dict[str, ParamSpec]

    def detect(self, geo: Geometry, params: dict, ctx: Context) -> list[Feature]: ...


@runtime_checkable
class Strategy(Protocol):
    """Geometry -> BlockStructure (topology and corner positions; discretisation and fill come later)."""

    name: str
    params: dict[str, ParamSpec]

    def applicable(self, geo: Geometry, features: list[Feature]) -> Score: ...

    def run(self, geo: Geometry, features: list[Feature], params: dict, ctx: Context) -> BlockStructure: ...


@runtime_checkable
class Op(Protocol):
    """BlockStructure -> BlockStructure (collapse, ogrid/pillow, split)."""

    name: str
    params: dict[str, ParamSpec]

    def applicable(self, bs: BlockStructure, geo: Geometry) -> Score: ...

    def apply(self, bs: BlockStructure, geo: Geometry, params: dict, ctx: Context) -> BlockStructure: ...


@runtime_checkable
class Discretiser(Protocol):
    """Sets edge_cells (sheet-consistent) and edge_grading."""

    name: str
    params: dict[str, ParamSpec]

    def apply(self, bs: BlockStructure, geo: Geometry, params: dict, ctx: Context) -> BlockStructure: ...


@runtime_checkable
class Filler(Protocol):
    """Builds the fine grids (``grids``, ``fine_points``) of a discretised structure."""

    name: str
    params: dict[str, ParamSpec]

    def fill(self, bs: BlockStructure, geo: Geometry, params: dict, ctx: Context) -> BlockStructure: ...


@runtime_checkable
class Exporter(Protocol):
    name: str
    params: dict[str, ParamSpec]
    formats: tuple[str, ...]

    def write(self, bs: BlockStructure, path: Path, params: dict, ctx: Context) -> list[Path]: ...


@runtime_checkable
class InterfacePlugin(Protocol):
    """Coupling between two components. Step 1 only asks whether a pairing is possible; patch and
    solver-dict export follow in roadmap step 2."""

    name: str
    params: dict[str, ParamSpec]

    def requirements(self, a: Geometry, b: Geometry, params: dict) -> Score: ...


CONTRACTS: dict[str, type] = {
    "detectors": Detector,
    "strategies": Strategy,
    "ops": Op,
    "discretisers": Discretiser,
    "fillers": Filler,
    "exporters": Exporter,
    "interfaces": InterfacePlugin,
}


def validate(kind: str, obj: Any) -> list[str]:
    """Contract violations of ``obj`` as plugin kind ``kind`` (empty list: fine)."""
    proto = CONTRACTS[kind]
    problems = []
    if not isinstance(getattr(obj, "name", None), str) or not obj.name:
        problems.append("missing a non-empty 'name'")
    params = getattr(obj, "params", None)
    if not isinstance(params, dict):
        problems.append("missing a 'params' dict of ParamSpec")
    else:
        bad = [k for k, v in params.items() if not isinstance(v, ParamSpec)]
        if bad:
            problems.append(f"params not given as ParamSpec: {bad}")
    for meth in _methods(proto):
        if not callable(getattr(obj, meth, None)):
            problems.append(f"missing method '{meth}'")
    if kind == "exporters" and not isinstance(getattr(obj, "formats", None), tuple):
        problems.append("missing 'formats' tuple")
    return problems


def resolve_params(spec: dict[str, ParamSpec], given: dict | None) -> dict:
    """Defaults filled in, unknown keys and out-of-range values rejected with a clear message."""
    given = dict(given or {})
    unknown = sorted(set(given) - set(spec))
    if unknown:
        raise ValueError(f"unknown parameter(s) {unknown}; known: {sorted(spec)}")
    out = {}
    for k, s in spec.items():
        v = given.get(k, s.default)
        if s.choices is not None and v not in s.choices:
            raise ValueError(f"{k}={v!r} not in {s.choices}")
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            if s.lower is not None and v < s.lower or s.upper is not None and v > s.upper:
                raise ValueError(f"{k}={v} outside [{s.lower}, {s.upper}]")
        if s.kind == "path" and v:
            # stages run in their own working directory: a relative path would point elsewhere there
            v = Path(os.path.expandvars(str(v))).expanduser()
            if not v.exists():
                raise FileNotFoundError(f"{k}: no such file or directory: {v}")
            v = str(v.absolute())
        out[k] = v
    return out


def describe(obj: Any) -> dict:
    """JSON-ready description of a plugin: name, doc, params (for the CLI, the MCP server, a UI)."""
    doc = (getattr(obj, "__doc__", None) or type(obj).__doc__ or "").strip()
    return {
        "name": obj.name,
        "doc": doc.split("\n\n")[0],
        "params": {k: {kk: vv for kk, vv in asdict(v).items() if vv is not None} for k, v in obj.params.items()},
    }


def _methods(proto) -> list[str]:
    return [m for m, v in vars(proto).items() if callable(v) and not m.startswith("_")]
