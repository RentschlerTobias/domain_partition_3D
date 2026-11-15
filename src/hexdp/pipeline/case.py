"""Case files: one YAML per blocking job.

    name: canada_dt_sweep
    geometry: ${HEXDP_DATA}/canada/domain.brep     # .msh / .step / .brep / .stl; ${VAR} and paths relative
    mesh_size: null                                 #   to the YAML file are resolved
    roles: {12: wall, RU_INLET: inlet}              # optional manual patch roles (surface tag or name)
    detectors: []                                   # [{name, params}]
    strategy: sweep                                 # a name, a list of candidates, or "auto" (all registered)
    params: {sweep: {h2d: 0.1}}                     # per-plugin parameters, by plugin name
    ops: [block_collapse]                           # names or {name, params}
    discretiser: null                               # optional
    filler: null                                    # optional
    exporters: [vtk]
    reuse: []                                       # files copied into the workspace before the run
    resources: {cores: 8, cpu_list: "0-7", timeout_s: null}
    workspace: runs/canada_dt_sweep                 # default runs/<name> next to the YAML
    reference: {blocks: 72, min_sj: 0.0588}         # optional, compared by `hexdp eval`

A plugin name may also be a direct reference "package.module:Class" (plugins not yet registered).
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from pathlib import Path

import yaml

KEYS = {"name", "geometry", "mesh_size", "roles", "detectors", "strategy", "params", "ops", "discretiser",
        "filler", "exporters", "reuse", "resources", "workspace", "reference", "description"}


@dataclass
class StageSpec:
    name: str
    params: dict = field(default_factory=dict)


@dataclass
class Case:
    name: str
    geometry: Path
    path: Path | None = None
    description: str = ""
    mesh_size: float | None = None
    roles: dict = field(default_factory=dict)
    detectors: list[StageSpec] = field(default_factory=list)
    strategies: list[StageSpec] = field(default_factory=list)   # empty = all registered ("auto")
    ops: list[StageSpec] = field(default_factory=list)
    discretiser: StageSpec | None = None
    filler: StageSpec | None = None
    exporters: list[StageSpec] = field(default_factory=list)
    reuse: list[Path] = field(default_factory=list)
    cores: int = 1
    cpu_list: str | None = None
    timeout_s: float | None = None
    workspace: Path | None = None
    reference: dict = field(default_factory=dict)


class CaseError(ValueError):
    pass


_VAR = re.compile(r"\$\{(\w+)\}|\$(\w+)")


def _expand(value: str, where: str) -> str:
    def sub(m):
        var = m.group(1) or m.group(2)
        if var not in os.environ:
            raise CaseError(f"{where}: environment variable ${{{var}}} is not set")
        return os.environ[var]
    return _VAR.sub(sub, value)


def _path(value, base: Path, where: str) -> Path:
    p = Path(_expand(str(value), where)).expanduser()
    return p if p.is_absolute() else (base / p).resolve()


def _expand_params(d, where):
    """${VAR} in string parameter values, recursively."""
    if isinstance(d, dict):
        return {k: _expand_params(v, f"{where}.{k}") for k, v in d.items()}
    if isinstance(d, list):
        return [_expand_params(v, where) for v in d]
    return _expand(d, where) if isinstance(d, str) else d


def _stages(value, params: dict, where: str) -> list[StageSpec]:
    if value in (None, [], "none"):
        return []
    items = value if isinstance(value, list) else [value]
    out = []
    for it in items:
        if isinstance(it, str):
            out.append(StageSpec(it, _expand_params(dict(params.get(it, {})), f"params.{it}")))
        elif isinstance(it, dict) and "name" in it:
            merged = dict(params.get(it["name"], {}))
            merged.update(it.get("params") or {})
            out.append(StageSpec(it["name"], _expand_params(merged, f"{where}.{it['name']}")))
        else:
            raise CaseError(f"{where}: expected a plugin name or {{name, params}}, got {it!r}")
    return out


def parse(data: dict, base: Path, path: Path | None = None) -> Case:
    unknown = sorted(set(data) - KEYS)
    if unknown:
        raise CaseError(f"unknown key(s) {unknown}; known: {sorted(KEYS)}")
    for k in ("name", "geometry"):
        if not data.get(k):
            raise CaseError(f"missing '{k}'")
    params = data.get("params") or {}
    strat = data.get("strategy", "auto")
    res = data.get("resources") or {}
    single = lambda key: (_stages(data.get(key), params, key) or [None])[0]  # noqa: E731
    case = Case(
        name=str(data["name"]),
        geometry=_path(data["geometry"], base, "geometry"),
        path=path,
        description=str(data.get("description", "")),
        mesh_size=data.get("mesh_size"),
        roles=dict(data.get("roles") or {}),
        detectors=_stages(data.get("detectors"), params, "detectors"),
        strategies=[] if strat == "auto" else _stages(strat, params, "strategy"),
        ops=_stages(data.get("ops"), params, "ops"),
        discretiser=single("discretiser"),
        filler=single("filler"),
        exporters=_stages(data.get("exporters"), params, "exporters"),
        reuse=[_path(p, base, "reuse") for p in data.get("reuse") or []],
        cores=int(res.get("cores", 1)),
        cpu_list=res.get("cpu_list"),
        timeout_s=res.get("timeout_s"),
        workspace=_path(data["workspace"], base, "workspace") if data.get("workspace")
        else (base / "runs" / str(data["name"])).resolve(),
        reference=dict(data.get("reference") or {}),
    )
    if not case.geometry.exists():
        raise CaseError(f"geometry not found: {case.geometry}")
    missing = [str(p) for p in case.reuse if not p.exists()]
    if missing:
        raise CaseError(f"reuse file(s) not found: {missing}")
    return case


def load_case(path) -> Case:
    path = Path(path).resolve()
    try:
        with open(path) as fh:
            data = yaml.safe_load(fh) or {}
    except yaml.YAMLError as exc:
        raise CaseError(f"{path}: not valid YAML: {exc}") from exc
    if not isinstance(data, dict):
        raise CaseError(f"{path}: a case file is a mapping")
    return parse(data, path.parent, path)
