"""Run a case: load the geometry and plugins, run the chain in the case workspace, write the report.

Workspace layout (``case.workspace``):

    status.json      state (running / done / failed), pid, start, end, error
    report.json      strategy, scores, timings, block-structure summary, outputs, reference comparison
    decisions.json   the decision log of the run
    <stage logs and outputs of the plugins>

A run reuses what the workspace already holds (each plugin skips stages whose output exists), so
running a case again resumes it; ``fresh=True`` refuses to run into a non-empty workspace.
"""

from __future__ import annotations

import json
import os
import shutil
import time
import traceback
from pathlib import Path

from hexdp import __version__, registry
from hexdp.contracts import Context
from hexdp.core.geometry import PatchRole, RoleSource
from hexdp.pipeline.case import Case
from hexdp.pipeline.chain import Stage, run_chain


def _write(path: Path, data: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=1, default=str))
    os.replace(tmp, path)


def status(workspace) -> dict:
    p = Path(workspace) / "status.json"
    if not p.exists():
        return {"state": "new"}
    st = json.loads(p.read_text())
    if st.get("state") == "running" and st.get("pid") and not _alive(st["pid"]):
        st["state"] = "died"
    return st


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
        return True
    except OSError:
        return False


def _stages(kind: str, specs) -> list[Stage]:
    return [Stage(registry.load(kind, s.name), s.params) for s in specs]


def compare_reference(summary: dict, ref: dict) -> dict:
    """Reference checks: blocks exact; min_sj within abs 1e-3 unless ref gives min_sj_tol."""
    out = {}
    if "blocks" in ref:
        out["blocks"] = {"ref": ref["blocks"], "got": summary["blocks"], "ok": summary["blocks"] == ref["blocks"]}
    if "min_sj" in ref:
        got = summary.get("min_sj_fine")
        tol = ref.get("min_sj_tol", 1e-3)
        out["min_sj"] = {"ref": ref["min_sj"], "got": got, "ok": got is not None and abs(got - ref["min_sj"]) <= tol}
    return out


def run_case(case: Case, fresh: bool = False) -> dict:
    ws = case.workspace
    if fresh and ws.exists() and any(ws.iterdir()):
        raise FileExistsError(f"workspace not empty: {ws}")
    ws.mkdir(parents=True, exist_ok=True)
    st = status(ws)
    if st.get("state") == "running":
        raise RuntimeError(f"case '{case.name}' is already running (pid {st.get('pid')}) in {ws}")
    t0 = time.time()
    _write(ws / "status.json", {"state": "running", "pid": os.getpid(), "start": t0, "case": str(case.path)})
    report = {"case": case.name, "case_file": str(case.path), "geometry": str(case.geometry),
              "hexdp": __version__, "workspace": str(ws)}
    try:
        for f in case.reuse:
            dst = ws / f.name
            if not dst.exists():
                (shutil.copytree if f.is_dir() else shutil.copy2)(f, dst)
        from hexdp.io.gmsh_model import load

        geo = load(case.geometry, mesh_size=case.mesh_size)
        for key, role in case.roles.items():
            tags = [t for t, s in geo.surfaces.items() if str(t) == str(key) or s.name == key]
            if not tags:
                raise ValueError(f"roles: no surface with tag or name '{key}'")
            geo.set_role(tags, PatchRole(role), RoleSource.MANUAL)
        report["geometry_summary"] = geo.summary()
        ctx = Context(ws, cores=case.cores, timeout_s=case.timeout_s, cpu_list=case.cpu_list)
        strategies = (_stages("strategies", case.strategies) if case.strategies
                      else [Stage(registry.load("strategies", n)) for n in registry.names("strategies")])
        res = run_chain(
            geo, ctx,
            detectors=_stages("detectors", case.detectors),
            strategies=strategies,
            ops=_stages("ops", case.ops),
            discretiser=_stages("discretisers", [case.discretiser])[0] if case.discretiser else None,
            filler=_stages("fillers", [case.filler])[0] if case.filler else None,
            exporters=_stages("exporters", case.exporters),
        )
        summary = res.blocks.summary()
        report.update(state="done", strategy=res.strategy,
                      scores={n: {"value": s.value, "reasons": s.reasons, "missing": s.missing} for n, s in res.scores.items()},
                      features=[{"kind": f.kind, "confidence": f.confidence, "reasons": f.reasons} for f in res.features],
                      timings=res.timings, outputs=res.outputs, blocks=summary,
                      reference=compare_reference(summary, case.reference))
        res.blocks.save(ws / "blocks.json")
    except Exception as exc:
        report.update(state="failed", error=f"{type(exc).__name__}: {exc}", traceback=traceback.format_exc())
    report["time_s"] = round(time.time() - t0, 2)
    _write(ws / "report.json", report)
    _write(ws / "status.json", {"state": report["state"], "pid": os.getpid(), "start": t0, "end": time.time(),
                                "case": str(case.path), "error": report.get("error")})
    return report
