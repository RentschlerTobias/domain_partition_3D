"""MCP server: the hexdp tools for an agent (stdio transport).

    hexdp mcp            (or: python -m hexdp.mcp.server)

Environment:
    HEXDP_DATA        input data root (datasets, geometries); readable
    HEXDP_MCP_ROOT    job workspaces (default ./hexdp_mcp_runs); readable and writable
    HEXDP_MCP_CPUS    CPU list for jobs (default "0-7"), HEXDP_MCP_CORES core count (default 8)
    HEXDP_MCP_TIMEOUT time limit per job stage in seconds (default 7200)

Workflow of an agent:
    list_inputs -> open_geometry -> inspect_geometry / render_geometry -> detect_features ->
    list_strategies -> run_strategy (background job) -> wait_job / job_status / job_log -> job_result ->
    evaluate / render_blocks / render_quality -> block_op -> export -> record_decision ->
    save_workflow -> validate_workflow (other samples of the family) -> wait_job -> job_result

Guards: inputs are read only below HEXDP_DATA, HEXDP_MCP_ROOT or the working directory; paths
containing 'charite' are refused (confidential data); jobs run with the configured CPU list and only one
computation (job or block_op) runs at a time; there is no network tool. Large data stay files: tools return ids, paths and small JSON.
"""

from __future__ import annotations

import collections
import functools
import json
import os
import re
import time
import uuid
from pathlib import Path

import numpy as np
import yaml
from mcp.server.mcpserver import Image, MCPServer
from mcp.server.mcpserver.exceptions import ToolError

from hexdp import __version__, registry
from hexdp.contracts import Context, describe, resolve_params
from hexdp.core.blockstructure import invalid
from hexdp.core.geometry import PatchRole, RoleSource

GEOMETRY_SUFFIXES = {".msh", ".step", ".stp", ".iges", ".igs", ".brep", ".stl", ".vtk"}

INSTRUCTIONS = """hexdp builds hexahedral block structures (domain partitions) for CFD meshes.
Goal for a geometry: a valid block structure (no inverted cells, conforming) with as few blocks as the
geometry allows. Typical order: open_geometry, inspect_geometry, detect_features (detectors
'rotation_axis', 'sweep'), list_strategies (scores with reasons), run_strategy (a background job;
wait_job / job_status), job_result, evaluate, optionally block_op ('block_collapse'), export. Explain every
choice with record_decision. Strategies: 'sweep' (rotational sectors, CAD input), 'cgrid' (sweep with
C-layers at thin walls, needs param wall_cad), 'algohex' (general topology, dtOO .msh input; AlgoHex can
take hours - reuse an existing hex_hex_x.ovm via the 'reuse' argument when one is available).
File parameters (kind 'path', e.g. wall_cad) may be given relative to the data root, as list_inputs shows them.
Only one computation runs at a time: run_strategy and block_op are refused while a job is running.
Wait for a job with wait_job (blocks up to max_s seconds) instead of polling job_status.
For a geometry family (a parametrised geometry, e.g. dataset samples): once a sample works, save_workflow turns
the job and the block_ops applied after it into a workflow file with placeholders (e.g. the sample name);
validate_workflow runs it on other samples; adapt and save again until it is accepted on all of them. The
workflow then runs without an agent: hexdp run workflow.yaml --set sample=..."""


class State:
    def __init__(self):
        self.data = Path(os.environ.get("HEXDP_DATA", os.getcwd())).resolve()
        self.root = Path(os.environ.get("HEXDP_MCP_ROOT", "hexdp_mcp_runs")).resolve()
        self.cpus = os.environ.get("HEXDP_MCP_CPUS", "0-7")
        self.cores = int(os.environ.get("HEXDP_MCP_CORES", "8"))
        self.timeout = float(os.environ.get("HEXDP_MCP_TIMEOUT", "7200"))      # per stage, seconds
        self.root.mkdir(parents=True, exist_ok=True)
        os.environ.setdefault("HEXDP_DATA", str(self.data))     # ${HEXDP_DATA} in saved workflows
        self.geos: dict[str, dict] = {}       # geo_id -> {"geo", "path", "features", "detectors"}
        self.blocks: dict[str, dict] = {}     # bs_id -> {"bs", "geo_id", "origin", "job", "ops", "exports"}

    def readable(self, p) -> Path:
        """The path as seen below an allowed root (symlinked data directories are allowed, as the
        administrator put them there); '..' is normalised away before the check."""
        p = Path(os.path.expandvars(str(p))).expanduser()
        p = Path(os.path.normpath(p if p.is_absolute() else self.data / p))
        real = p.resolve()
        if "charite" in str(p).lower() or "charite" in str(real).lower():
            raise PermissionError("Charité data are confidential and not available to agents")
        roots = [self.data, self.root, Path(os.getcwd()).resolve()]
        if not any(q == r or r in q.parents for r in roots for q in (p, real)):
            raise PermissionError(f"{p} is outside the allowed roots {[str(r) for r in roots]}")
        if not p.exists():
            raise FileNotFoundError(str(p))
        return p

    def geo(self, geo_id):
        if geo_id not in self.geos:
            raise KeyError(f"unknown geo_id '{geo_id}'; open_geometry first (known: {sorted(self.geos)})")
        return self.geos[geo_id]

    def bs(self, bs_id):
        if bs_id not in self.blocks:
            raise KeyError(f"unknown bs_id '{bs_id}' (known: {sorted(self.blocks)})")
        return self.blocks[bs_id]

    def running_jobs(self) -> list[str]:
        """Jobs under the root whose process is still alive."""
        from hexdp.pipeline.runner import status

        return [d.name for d in self.root.iterdir()
                if d.is_dir() and (d / "status.json").exists() and status(d).get("state") == "running"]

    def require_idle(self, what: str) -> None:
        """One computation at a time on the job CPU list (the core budget must not depend on the agent)."""
        busy = self.running_jobs()
        if busy:
            raise RuntimeError(f"{what} refused: job {busy[0]} is still running on CPUs {self.cpus}; "
                               "wait for it (job_status) before starting another computation")

    def job_dir(self, job_id) -> Path:
        d = (self.root / job_id).resolve()
        if self.root not in d.parents or not d.is_dir():
            raise KeyError(f"unknown job '{job_id}'")
        return d


S: State | None = None
server = MCPServer(name="hexdp", version=__version__, instructions=INSTRUCTIONS)


def tool(fn):
    """Register fn as an MCP tool; expected errors reach the agent as their message (the SDK hides the
    message of an unexpected exception)."""
    @functools.wraps(fn)
    def wrapped(*a, **kw):
        try:
            return fn(*a, **kw)
        except ToolError:
            raise
        except (PermissionError, FileNotFoundError, KeyError, ValueError, FileExistsError, RuntimeError) as exc:
            msg = exc.args[0] if isinstance(exc, KeyError) and exc.args else str(exc)
            raise ToolError(f"{type(exc).__name__}: {msg}") from exc
    return server.tool()(wrapped)


def state() -> State:
    global S
    if S is None:
        S = State()
    return S


def _new_id(prefix):
    return f"{prefix}_{uuid.uuid4().hex[:8]}"


def _geo_of(entry: dict):
    """The open geometry of a block-structure entry (for the boundary topology check), or None."""
    g = state().geos.get(entry.get("geo_id")) if entry.get("geo_id") else None
    return g["geo"] if g else None


def _paths(plugin, params: dict | None) -> dict:
    """File parameters (ParamSpec kind 'path') through the read guard: relative to HEXDP_DATA, inside the
    allowed roots, never Charite; absolute afterwards, so a stage in another working directory finds them."""
    params = dict(params or {})
    for k, spec in plugin.params.items():
        if spec.kind == "path" and params.get(k):
            params[k] = str(state().readable(params[k]))
    return params


def _jsonable(x, depth=0):
    if isinstance(x, dict):
        return {str(k): _jsonable(v, depth + 1) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        if len(x) > 60:
            return [_jsonable(v, depth + 1) for v in x[:60]] + [f"... {len(x) - 60} more"]
        return [_jsonable(v, depth + 1) for v in x]
    if isinstance(x, np.generic):
        return x.item()
    if isinstance(x, np.ndarray):
        return _jsonable(x.tolist(), depth)
    return x


# --------------------------------------------------------------------------------------------- inputs
@tool
def list_inputs(subdir: str = "", pattern: str = "*") -> dict:
    """List geometry files, parameter files of parametric designs (params*.json, e.g. dtOO runners: strategy
    'runner') and case YAMLs below HEXDP_DATA/subdir, recursively (max 200)."""
    s = state()
    base = s.readable(subdir) if subdir else s.data
    import fnmatch

    found = []
    for dirpath, dirnames, filenames in os.walk(base, followlinks=True):   # data dirs are often symlinks
        dirnames[:] = sorted(d for d in dirnames if not d.startswith(".") and d != "runs" and "charite" not in d.lower())
        found += [Path(dirpath) / f for f in sorted(filenames) if fnmatch.fnmatch(f, pattern)]
    out = []
    for p in found:
        if p.is_file() and (p.suffix.lower() in GEOMETRY_SUFFIXES or p.suffix in (".yaml", ".yml")
                            or (p.suffix == ".json" and p.stem.startswith("params"))):   # parametric designs
            if "charite" in str(p).lower():
                continue
            out.append({"path": str(p.relative_to(s.data)) if str(p).startswith(str(s.data) + os.sep) else str(p),
                        "size_mb": round(p.stat().st_size / 1e6, 2)})
            if len(out) >= 200:
                break
    return {"root": str(base), "files": out}


@tool
def open_geometry(path: str, mesh_size: float | None = None) -> dict:
    """Load a geometry (.msh, .step/.stp/.brep/.iges, .stl) -> geo_id and a summary."""
    from hexdp.io.gmsh_model import load

    s = state()
    p = s.readable(path)
    geo = load(p, mesh_size=mesh_size)
    gid = _new_id("geo")
    s.geos[gid] = {"geo": geo, "path": str(p), "features": []}
    return {"geo_id": gid, "path": str(p), "summary": _jsonable(geo.summary())}


@tool
def inspect_geometry(geo_id: str) -> dict:
    """Surfaces (tag, name, role and its source, area, bbox) and curves by kind."""
    g = state().geo(geo_id)["geo"]
    surfaces = []
    for t, srf in sorted(g.surfaces.items()):
        lo, hi = srf.nodes.min(0), srf.nodes.max(0)
        surfaces.append({"tag": t, "name": srf.name, "role": srf.role.value, "role_source": srf.role_source.value,
                         "area": round(srf.area, 6), "bbox": [np.round(lo, 4).tolist(), np.round(hi, 4).tolist()]})
    kinds = {}
    for c in g.curves.values():
        kinds.setdefault(c.kind.value, []).append(c.tag)
    return {"summary": _jsonable(g.summary()), "surfaces": surfaces, "curves_by_kind": kinds}


@tool
def set_roles(geo_id: str, roles: dict[str, str]) -> dict:
    """Set patch roles by surface tag or name: {"12": "inlet", "RU_BLADE": "wall"}. Roles: wall, inlet,
    outlet, periodic, interface, symmetry, unknown."""
    g = state().geo(geo_id)["geo"]
    done = {}
    for key, role in roles.items():
        tags = [t for t, srf in g.surfaces.items() if str(t) == str(key) or srf.name == key]
        if not tags:
            raise KeyError(f"no surface with tag or name '{key}'")
        g.set_role(tags, PatchRole(role), RoleSource.MANUAL)
        done[key] = tags
    return {"set": done, "roles": _jsonable(g.summary()["roles"])}


@tool
def render_geometry(geo_id: str, elev: float = 25, azim: float = -60) -> Image:
    """PNG of the boundary coloured by patch role, feature curves in black."""
    from hexdp.viz.render import render_geometry as rg

    return Image(data=rg(state().geo(geo_id)["geo"], elev=elev, azim=azim), format="png")


# --------------------------------------------------------------------------------------------- plugins
@tool
def list_plugins(kind: str | None = None) -> dict:
    """Registered plugins with their doc and parameter schema (kinds: detectors, strategies, ops,
    discretisers, fillers, exporters, interfaces)."""
    kinds = [kind] if kind else list(registry.KINDS)
    return {k: [describe(registry.load(k, n)) for n in registry.names(k)] for k in kinds}


@tool
def detect_features(geo_id: str, detectors: list[str] | None = None) -> dict:
    """Run detectors (default: all registered) and keep their features for list_strategies / run_strategy."""
    s = state()
    entry = s.geo(geo_id)
    names = detectors or registry.names("detectors")
    ctx = Context(s.root / "_detect" / geo_id, cores=1)
    feats = []
    for n in names:
        if not entry["geo"].surfaces:        # a parameter file: the strategy builds the geometry
            break
        d = registry.load("detectors", n)
        feats += d.detect(entry["geo"], resolve_params(d.params, {}), ctx)
    entry["features"] = feats
    entry["detectors"] = list(names)
    return {"features": [{"kind": f.kind, "confidence": f.confidence, "reasons": f.reasons,
                          "data": _jsonable({k: v for k, v in f.data.items() if k not in ("rms_rel_facet",)})}
                         for f in feats]}


@tool
def list_strategies(geo_id: str) -> dict:
    """Every registered strategy with its applicability score, reasons and what it is missing, best first
    (uses the features of the last detect_features)."""
    entry = state().geo(geo_id)
    rows = []
    for n in registry.names("strategies"):
        st = registry.load("strategies", n)
        sc = st.applicable(entry["geo"], entry["features"])
        rows.append({"name": n, "score": sc.value, "applicable": sc.ok, "reasons": sc.reasons, "missing": sc.missing})
    rows.sort(key=lambda r: -r["score"])
    return {"features_used": [f.kind for f in entry["features"]], "strategies": rows}


# --------------------------------------------------------------------------------------------- jobs
@tool
def run_strategy(geo_id: str, strategy: str, params: dict | None = None, ops: list[str] | None = None,
                 discretiser: str | None = None, discretiser_params: dict | None = None,
                 filler: str | None = None, filler_params: dict | None = None, reuse: list[str] | None = None,
                 detectors: list[str] | None = None, name: str | None = None) -> dict:
    """Start a blocking job in the background: strategy (name), its params, optional ops (e.g.
    ["block_collapse"]), filler (e.g. "tfi"), files to reuse in the workspace (e.g. an existing
    hex_hex_x.ovm). Returns job_id; poll job_status."""
    from hexdp.pipeline.jobs import start

    s = state()
    entry = s.geo(geo_id)
    params = _paths(registry.load("strategies", strategy), params)   # fails early on an unknown name
    if discretiser:
        discretiser_params = _paths(registry.load("discretisers", discretiser), discretiser_params)
    if filler:
        filler_params = _paths(registry.load("fillers", filler), filler_params)
    s.require_idle("run_strategy")
    job = name or _new_id(re.sub(r"\W+", "_", strategy.split(":")[-1]).lower())
    ws = s.root / job
    ws.mkdir(parents=True, exist_ok=False)
    case = {"name": job, "geometry": entry["path"], "strategy": [{"name": strategy, "params": params or {}}],
            "detectors": detectors if detectors is not None else entry.get("detectors", []),
            "ops": ops or [], "exporters": ["vtk"], "workspace": str(ws),
            "reuse": [str(s.readable(r)) for r in reuse or []],
            "resources": {"cores": s.cores, "cpu_list": s.cpus, "timeout_s": s.timeout}}
    if discretiser:
        case["discretiser"] = {"name": discretiser, "params": discretiser_params or {}}
    if filler:
        case["filler"] = {"name": filler, "params": filler_params or {}}
    roles = {str(t): srf.role.value for t, srf in entry["geo"].surfaces.items() if srf.role_source is RoleSource.MANUAL}
    if roles:
        case["roles"] = roles
    cf = ws / "case.yaml"
    cf.write_text(yaml.safe_dump(case, sort_keys=False))
    start(cf)
    return {"job_id": job, "workspace": str(ws), "case": case}


@tool
def job_status(job_id: str) -> dict:
    """State of a job (running / done / failed / died), elapsed time, last log line of the active stage."""
    from hexdp.pipeline.runner import status

    d = state().job_dir(job_id)
    st = status(d)
    logs = sorted((p for p in d.glob("*.log")), key=lambda p: p.stat().st_mtime)
    last = ""
    if logs:
        lines = [ln for ln in logs[-1].read_text(errors="replace").splitlines() if ln.strip()]
        last = lines[-1][:300] if lines else ""
    if (d / "rows.json").exists():                                  # a validation job: runs finished so far
        st["runs_done"] = len(json.loads((d / "rows.json").read_text())["rows"])
        st["runs_total"] = len(json.loads((d / "runs.json").read_text()))
    st["elapsed_s"] = round((st.get("end") or time.time()) - st["start"], 1) if st.get("start") else None
    return {"job_id": job_id, **_jsonable(st), "active_log": logs[-1].name if logs else None, "last_line": last}


@tool
def wait_job(job_id: str, max_s: float = 600) -> dict:
    """Wait until the job is no longer running, or at most max_s seconds (capped at 1800); returns
    job_status. Use this instead of polling job_status in a loop."""
    from hexdp.pipeline.runner import status

    d = state().job_dir(job_id)
    t_end = time.time() + min(max(float(max_s), 0.0), 1800.0)
    while status(d).get("state") == "running" and time.time() < t_end:
        time.sleep(min(5.0, max(t_end - time.time(), 0.0)))
    return job_status(job_id)


@tool
def job_log(job_id: str, log: str | None = None, tail: int = 40) -> dict:
    """Tail of a job log (default: the most recent one); lists the available logs."""
    d = state().job_dir(job_id)
    logs = sorted((p for p in d.glob("*.log")), key=lambda p: p.stat().st_mtime)
    pick = d / log if log else (logs[-1] if logs else d / "run.out")
    text = pick.read_text(errors="replace").splitlines()[-tail:] if pick.exists() else []
    return {"logs": [p.name for p in logs], "log": pick.name, "lines": text}


@tool
def job_result(job_id: str) -> dict:
    """Report of a finished job and a bs_id of its block structure (for evaluate / render / block_op)."""
    from hexdp.io.blocks_vtk import read_blocks_vtk

    s = state()
    d = s.job_dir(job_id)
    rep_file = d / "report.json"
    if not rep_file.exists():
        return {"job_id": job_id, "state": "running or not started", "hint": "poll job_status"}
    rep = json.loads(rep_file.read_text())
    if rep.get("kind") == "batch":
        return {"job_id": job_id, **_jsonable({k: rep.get(k) for k in ("state", "workflow", "runs", "accepted", "time_s")}),
                "rows": [_jsonable({k: r.get(k) for k in ("parameters", "state", "blocks", "valid", "min_sj", "issues",
                                                          "accept", "time_s", "error", "workspace")}) for r in rep["rows"]]}
    out = {k: rep.get(k) for k in ("state", "strategy", "scores", "timings", "blocks", "outputs", "error", "time_s")}
    if rep.get("state") == "done":
        vtk = (rep.get("outputs") or {}).get("vtk", [None])[0]
        bs = read_blocks_vtk(vtk) if vtk else None
        if bs is not None:
            bid = _new_id("bs")
            geo_id = next((g for g, e in s.geos.items() if e["path"] == rep.get("geometry")), None)
            s.blocks[bid] = {"bs": bs, "geo_id": geo_id, "origin": job_id, "job": job_id, "ops": [], "exports": []}
            out["bs_id"] = bid
            bad = invalid(bs.check(geo=_geo_of(s.blocks[bid])))
            out["valid"] = not bad
            if bad:
                out["invalid_issues"] = dict(collections.Counter(i.kind for i in bad))
    return {"job_id": job_id, **_jsonable(out)}


# --------------------------------------------------------------------------------------------- blocks
@tool
def evaluate(bs_id: str, max_issues: int = 20) -> dict:
    """Validity, block-structure summary (blocks, cells, min scaled Jacobian, sheets) and the first issues.
    Invalid: inverted cells, a face shared by more than two blocks, hanging vertices, sheet count
    mismatches, degenerate blocks, a boundary whose topology differs from the geometry's (a hole or a
    missing block; checked when the geometry boundary is closed). A t_junction (1:n block layout over a conforming fine mesh) is a warning."""
    entry = state().bs(bs_id)
    bs, geo = entry["bs"], _geo_of(entry)
    issues = bs.check(geo=geo)
    return {"valid": not invalid(issues), "summary": _jsonable(bs.summary(geo)),
            "issues": [{"kind": i.kind, "message": i.message, "where": list(i.where)} for i in issues[:max_issues]],
            "issues_total": len(issues)}


@tool
def render_blocks(bs_id: str, elev: float = 25, azim: float = -60, highlight: list[int] | None = None) -> Image:
    """PNG of the block edges coloured by block (highlighted blocks in red)."""
    from hexdp.viz.render import render_blocks as rb

    return Image(data=rb(state().bs(bs_id)["bs"], elev=elev, azim=azim, highlight=tuple(highlight or ())),
                 format="png")


@tool
def render_quality(bs_id: str) -> Image:
    """PNG histogram of the fine-cell scaled Jacobian."""
    from hexdp.viz.render import render_quality as rq

    return Image(data=rq(state().bs(bs_id)["bs"]), format="png")


@tool
def block_op(bs_id: str, op: str, params: dict | None = None) -> dict:
    """Apply an op (e.g. 'block_collapse') to a block structure, synchronously in its own workspace
    (block_collapse on the canada draft tube takes ~1 min). Returns the new bs_id and summary."""
    s = state()
    entry = s.bs(bs_id)
    if entry["geo_id"] is None:
        raise ValueError("this block structure has no open geometry; open_geometry with the job's geometry first")
    geo = s.geo(entry["geo_id"])["geo"]
    plugin = registry.load("ops", op)
    s.require_idle("block_op")
    sc = plugin.applicable(entry["bs"], geo)
    if not sc.ok:
        return {"applied": False, "reasons": sc.reasons, "missing": sc.missing}
    ctx = Context(s.root / _new_id(f"op_{op}"), cores=s.cores, cpu_list=s.cpus)
    new = plugin.apply(entry["bs"], geo, resolve_params(plugin.params, _paths(plugin, params)), ctx)
    bid = _new_id("bs")
    s.blocks[bid] = {"bs": new, "geo_id": entry["geo_id"], "origin": f"{op}({bs_id})", "job": entry.get("job"),
                     "ops": entry.get("ops", []) + [{"name": op, "params": _paths(plugin, params)}], "exports": []}
    return {"applied": True, "bs_id": bid, "workspace": str(ctx.workspace), "reasons": sc.reasons,
            "before": entry["bs"].n_blocks, "summary": _jsonable(new.summary(geo))}


@tool
def export(bs_id: str, format: str = "vtk", name: str | None = None, allow_invalid: bool = False) -> dict:
    """Write a block structure with an exporter ('vtk', 'foam'); returns the written paths. An invalid
    structure (see evaluate) is refused unless allow_invalid is set."""
    s = state()
    entry = s.bs(bs_id)
    bs = entry["bs"]
    bad = invalid(bs.check(geo=_geo_of(entry)))
    if bad and not allow_invalid:
        raise ValueError(f"block structure {bs_id} is invalid ({dict(collections.Counter(i.kind for i in bad))}); "
                         "fix it or pass allow_invalid=true")
    ex = registry.load("exporters", format)
    ctx = Context(s.root / "_export", cores=1, cpu_list=s.cpus)
    paths = ex.write(bs, ctx.path(name or bs_id), resolve_params(ex.params, {}), ctx)
    if format not in entry.setdefault("exports", []):
        entry["exports"].append(format)
    return {"written": [str(p) for p in paths]}


@tool
def record_decision(job_or_geo_id: str, text: str) -> dict:
    """Write the reasoning behind a choice into the job's (or geometry's) decision notes."""
    s = state()
    d = s.root / job_or_geo_id if (s.root / job_or_geo_id).is_dir() else s.root / "_notes" / job_or_geo_id
    d.mkdir(parents=True, exist_ok=True)
    f = d / "agent_decisions.md"
    with open(f, "a") as fh:
        fh.write(f"- {time.strftime('%Y-%m-%d %H:%M:%S')}: {text.strip()}\n")
    return {"file": str(f)}


# --------------------------------------------------------------------------------------------- workflows
def _workflow_file(name: str) -> Path:
    if not re.fullmatch(r"[A-Za-z0-9][\w.-]*", name):
        raise ValueError(f"workflow name {name!r}: letters, digits, '_', '-', '.' only")
    return state().root / "_workflows" / f"{name}.yaml"


@tool
def save_workflow(bs_id: str, name: str, placeholders: dict[str, str] | None = None, accept: dict | None = None,
                  notes: list[str] | None = None, overwrite: bool = False) -> dict:
    """Turn the job behind bs_id and the block_ops applied after it into a workflow file for a geometry
    family. placeholders {param: value}: every occurrence of value (e.g. the sample name 'tistos_s1_0131')
    becomes ${param}; choose a value that appears only where samples differ. accept: acceptance of each
    run, e.g. {"valid": true, "blocks": [10, 17], "min_sj": 0.02} (default {"valid": true}). notes: the
    reasoning (the job's record_decision notes are added). Returns the file and its text; run it on other
    samples with validate_workflow."""
    from hexdp.pipeline.workflow import build_workflow, dump_workflow

    s = state()
    entry = s.bs(bs_id)
    if not entry.get("job"):
        raise ValueError(f"{bs_id} has no job behind it")
    bad = invalid(entry["bs"].check(geo=_geo_of(entry)))
    if bad:
        raise ValueError(f"{bs_id} is invalid ({dict(collections.Counter(i.kind for i in bad))}); "
                         "a workflow is saved only from a valid block structure")
    f = _workflow_file(name)
    if f.exists() and not overwrite:
        raise FileExistsError(f"workflow {name} exists ({f}); pass overwrite=true to replace it")
    jd = s.job_dir(entry["job"])
    case = yaml.safe_load((jd / "case.yaml").read_text())
    for k, v in (placeholders or {}).items():
        if not re.fullmatch(r"\w+", k):
            raise ValueError(f"placeholder name {k!r}: letters, digits, '_' only")
        if not v or v not in yaml.safe_dump(case) + json.dumps(entry["ops"]):
            raise ValueError(f"placeholder {k}: value {v!r} does not occur in the job's case or ops")
    job_notes = (jd / "agent_decisions.md").read_text().splitlines() if (jd / "agent_decisions.md").exists() else []
    summ = entry["bs"].summary(_geo_of(entry))
    wf = build_workflow(case, name=name, ops=entry["ops"], placeholders=placeholders, data_root=s.data, accept=accept,
                        notes=[*job_notes, *(notes or [])], exporters=entry.get("exports", []),
                        origin={"job": entry["job"], "bs_id": bs_id, "geometry": case.get("geometry"),
                                "blocks": summ.get("blocks"), "min_sj": summ.get("min_sj_fine")})
    f.parent.mkdir(parents=True, exist_ok=True)
    dump_workflow(wf, f)
    return {"workflow": str(f), "name": name, "text": f.read_text()}


@tool
def validate_workflow(workflow: str, runs: list[dict[str, str]]) -> dict:
    """Run a workflow (a name from save_workflow or a path) for several parameter sets, e.g.
    [{"sample": "tistos_s1_0317"}, {"sample": "tistos_s1_0527"}], in one background job (sequentially);
    wait_job, then job_result lists per run: valid, blocks, min_sj, accept. Reuse files that a sample does
    not have are skipped; each validation job runs in fresh workspaces."""
    from hexdp.pipeline.case import load_case
    from hexdp.pipeline.jobs import start_batch

    s = state()
    f = _workflow_file(workflow) if not workflow.endswith((".yaml", ".yml")) else s.readable(workflow)
    if not f.exists():
        raise FileNotFoundError(f"no workflow {workflow} ({f})")
    if not runs:
        raise ValueError("runs: give at least one parameter set")
    for values in runs:                             # every input through the guards before anything starts
        c = load_case(f, values)
        s.readable(c.geometry)
        for r in c.reuse:
            s.readable(r)
        for kind, specs in (("strategies", c.strategies), ("ops", c.ops),
                            ("discretisers", [c.discretiser] if c.discretiser else []),
                            ("fillers", [c.filler] if c.filler else [])):
            for sp in specs:
                _paths(registry.load(kind, sp.name), sp.params)
    s.require_idle("validate_workflow")
    job = _new_id(f"validate_{f.stem}"[:40])
    start_batch(f, runs, s.root / job)
    return {"job_id": job, "workflow": str(f), "runs": len(runs)}


def main():
    server.run("stdio")


if __name__ == "__main__":
    main()
