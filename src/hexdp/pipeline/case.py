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
    reuse: []                                       # files copied into the workspace before the run;
                                                    #   {path, optional: true} is skipped when missing
    resources: {cores: 8, cpu_list: "0-7", timeout_s: null}
    workspace: runs/canada_dt_sweep                 # default runs/<name> next to the YAML
    reference: {blocks: 72, min_sj: 0.0588}         # optional, compared by `hexdp eval`

A plugin name may also be a direct reference "package.module:Class" (plugins not yet registered).

A workflow is a case with parameters: one file for a whole geometry family.

    parameters: {sample: tistos_s1_0131}            # name: default (null = required), or {default, doc}
    name: tistos_${sample}                          # ${name} of a parameter anywhere in a string value
    geometry: ${HEXDP_DATA}/tistos/${sample}/mesh.msh
    accept: {valid: true, blocks: [10, 17], min_sj: 0.02}   # acceptance of each run (blocks: n or [lo, hi])
    origin: {...}                                   # informational: where and how it was found
    notes: [...]                                    # the reasoning behind the choices

    hexdp run workflow.yaml --set sample=tistos_s1_0317
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from pathlib import Path

import yaml

KEYS = {"name", "geometry", "mesh_size", "roles", "detectors", "strategy", "params", "ops", "discretiser",
        "filler", "exporters", "reuse", "resources", "workspace", "reference", "description",
        "parameters", "accept", "origin", "notes"}


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
    parameters: dict = field(default_factory=dict)      # the values used (workflows)
    accept: dict = field(default_factory=dict)


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


_PARAM = re.compile(r"\$\{(\w+)\}")


def parameter_values(declared: dict, overrides: dict | None = None) -> dict:
    """Values of a workflow's parameters: --set overrides, else the declared default; unknown or missing ones fail."""
    declared = dict(declared or {})
    overrides = dict(overrides or {})
    unknown = sorted(set(overrides) - set(declared))
    if unknown:
        raise CaseError(f"unknown parameter(s) {unknown}; declared: {sorted(declared)}")
    out = {}
    for k, d in declared.items():
        default = d.get("default") if isinstance(d, dict) else d
        v = overrides.get(k, default)
        if v is None:
            raise CaseError(f"parameter '{k}' has no default; give it with --set {k}=...")
        out[k] = str(v)
    return out


def substitute(data, values: dict):
    """${name} of a parameter in every string value (environment variables are left for parse)."""
    if isinstance(data, dict):
        return {k: substitute(v, values) for k, v in data.items()}
    if isinstance(data, list):
        return [substitute(v, values) for v in data]
    if isinstance(data, str):
        return _PARAM.sub(lambda m: values.get(m.group(1), m.group(0)), data)
    return data


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


def _reuse(items, base: Path) -> list[Path]:
    out = []
    for it in items or []:
        spec = it if isinstance(it, dict) else {"path": it}
        p = _path(spec["path"], base, "reuse")
        if p.exists() or not spec.get("optional"):
            out.append(p)
    return out


def parse(data: dict, base: Path, path: Path | None = None, overrides: dict | None = None) -> Case:
    unknown = sorted(set(data) - KEYS)
    if unknown:
        raise CaseError(f"unknown key(s) {unknown}; known: {sorted(KEYS)}")
    values = parameter_values(data.get("parameters"), overrides)
    data = {k: v for k, v in substitute(data, values).items() if k not in ("parameters", "origin", "notes")}
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
        reuse=_reuse(data.get("reuse"), base),
        cores=int(res.get("cores", 1)),
        cpu_list=res.get("cpu_list"),
        timeout_s=res.get("timeout_s"),
        workspace=_path(data["workspace"], base, "workspace") if data.get("workspace")
        else (base / "runs" / str(data["name"])).resolve(),
        reference=dict(data.get("reference") or {}),
        parameters=values,
        accept=dict(data.get("accept") or {}),
    )
    if not case.geometry.exists():
        raise CaseError(f"geometry not found: {case.geometry}")
    missing = [str(p) for p in case.reuse if not p.exists()]
    if missing:
        raise CaseError(f"reuse file(s) not found: {missing}")
    return case


def load_case(path, overrides: dict | None = None, workspace_root=None) -> Case:
    """Load a case or workflow; overrides set workflow parameters, workspace_root puts the workspace at
    workspace_root/<name> (validation runs)."""
    path = Path(path).resolve()
    try:
        with open(path) as fh:
            data = yaml.safe_load(fh) or {}
    except yaml.YAMLError as exc:
        raise CaseError(f"{path}: not valid YAML: {exc}") from exc
    if not isinstance(data, dict):
        raise CaseError(f"{path}: a case file is a mapping")
    case = parse(data, path.parent, path, overrides)
    if workspace_root is not None:
        case.workspace = Path(workspace_root).resolve() / case.name
    return case
