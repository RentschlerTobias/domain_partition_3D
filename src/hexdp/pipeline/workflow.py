"""Workflows: a case found once (often by an agent) turned into a template for a geometry family.

    build_workflow(case, ...)   the case file of a job + the ops applied afterwards -> workflow dict
    run_batch(wf, runs, ...)    run the workflow for several parameter sets, one row per run

A workflow is a case file with ``parameters`` (see hexdp.pipeline.case): the geometry and the other
sample-specific paths carry ``${name}`` placeholders, and ``accept`` says when a run counts as good.
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path

import yaml

from hexdp import __version__
from hexdp.pipeline.case import CaseError, load_case
from hexdp.pipeline.evaluate import row
from hexdp.pipeline.runner import _write, run_case


def parse_cpu_list(spec: str) -> list[int]:
    out = set()
    for part in str(spec).split(","):
        a, _, b = part.partition("-")
        out |= set(range(int(a), int(b or a) + 1))
    return sorted(out)

HEADER = """# hexdp workflow '{name}': a block structure recipe for a geometry family.
#   hexdp run {file} --set {example}
#   hexdp batch {file} --each {first}=a,b,c
# Placeholders: ${{name}} = a parameter below, ${{HEXDP_DATA}} = the data root (environment).
"""


def _replace(data, pairs):
    if isinstance(data, dict):
        return {_replace(k, pairs) if isinstance(k, str) else k: _replace(v, pairs) for k, v in data.items()}
    if isinstance(data, list):
        return [_replace(v, pairs) for v in data]
    if isinstance(data, str):
        for old, new in pairs:
            data = data.replace(old, new)
    return data


def build_workflow(case: dict, *, name: str, ops=(), placeholders: dict | None = None, data_root=None,
                   accept: dict | None = None, notes=(), exporters=(), origin: dict | None = None) -> dict:
    """Workflow from a job's case dict and the ops applied to its result afterwards ([{name, params}]).

    placeholders {param: value}: every occurrence of value in a string becomes ${param} (choose values that
    occur only where the sample differs, e.g. 'tistos_s1_0131'); the data root becomes ${HEXDP_DATA}.
    Reuse files become optional (a new sample has none); the workspace defaults to runs/<name>."""
    placeholders = dict(placeholders or {})
    wf = {k: v for k, v in case.items() if k not in ("name", "workspace", "parameters", "origin", "notes")}
    wf["ops"] = list(case.get("ops") or []) + [o if isinstance(o, str) else {"name": o["name"], "params": o.get("params") or {}}
                                               for o in ops]
    ex = list(case.get("exporters") or [])
    wf["exporters"] = ex + [e for e in exporters if e not in ex]
    wf["reuse"] = [r if isinstance(r, dict) else {"path": r, "optional": True} for r in case.get("reuse") or []]
    wf["accept"] = dict(accept) if accept is not None else {"valid": True}
    pairs = []
    if data_root is not None:
        pairs.append((str(Path(data_root).resolve()), "${HEXDP_DATA}"))
    pairs += sorted(((v, f"${{{k}}}") for k, v in placeholders.items()), key=lambda p: -len(p[0]))
    wf = _replace(wf, pairs)
    out = {"parameters": {k: {"default": v} for k, v in placeholders.items()},
           "name": name + "".join(f"_${{{k}}}" for k in placeholders)}
    out.update(wf)
    out["origin"] = {"hexdp": __version__, "date": time.strftime("%Y-%m-%d %H:%M:%S"), **(origin or {})}
    if notes:
        out["notes"] = [str(n) for n in notes]
    return out


def dump_workflow(wf: dict, path) -> Path:
    path = Path(path)
    params = list(wf.get("parameters") or {})
    first = params[0] if params else "name"
    example = " ".join(f"{k}=..." for k in params) or "(no parameters)"
    head = HEADER.format(name=path.stem, file=path.name, example=example, first=first)
    path.write_text(head + yaml.safe_dump(wf, sort_keys=False, width=110))
    return path


def parse_each(spec: str) -> list[dict]:
    """'sample=a,b,c' or 'sample=@file' (one value per line) -> [{sample: a}, ...]."""
    if "=" not in spec:
        raise CaseError(f"--each expects name=v1,v2,... or name=@file, got {spec!r}")
    k, v = spec.split("=", 1)
    vals = ([ln.strip() for ln in Path(v[1:]).read_text().splitlines() if ln.strip() and not ln.startswith("#")]
            if v.startswith("@") else [x for x in v.split(",") if x])
    return [{k: x} for x in vals]


def run_batch(wf, runs: list[dict], *, common: dict | None = None, workspace_root=None, rerun: bool = False,
              status_dir=None, cpu_list: str | None = None) -> list[dict]:
    """Run the workflow once per parameter set (sequentially); a finished run is reused unless rerun.
    status_dir: write status.json / report.json there (a background validation job); cpu_list overrides the
    workflow's resources (the CPUs the user allowed for this run)."""
    sd = Path(status_dir) if status_dir else None
    t0 = time.time()
    if sd:
        sd.mkdir(parents=True, exist_ok=True)
        _write(sd / "status.json", {"state": "running", "pid": os.getpid(), "start": t0, "case": str(wf)})
    rows = []
    for values in runs:
        values = {**(common or {}), **values}
        try:
            case = load_case(wf, values, workspace_root)
            if cpu_list:
                case.cpu_list = cpu_list
                case.cores = len(parse_cpu_list(cpu_list))
        except (CaseError, FileNotFoundError) as exc:
            rows.append({"case": ",".join(f"{k}={v}" for k, v in values.items()), "state": "invalid",
                         "parameters": values, "error": str(exc)})
            continue
        rep_file = case.workspace / "report.json"
        if not rerun and rep_file.exists() and json.loads(rep_file.read_text()).get("state") == "done":
            report = json.loads(rep_file.read_text())
        else:
            report = run_case(case)
        rows.append({**row(report), "parameters": case.parameters, "workspace": str(case.workspace)})
        if sd:
            _write(sd / "rows.json", {"rows": rows})
    if sd:
        ok = sum(r.get("accept") == "ok" for r in rows)
        _write(sd / "report.json", {"kind": "batch", "state": "done", "workflow": str(wf), "rows": rows,
                                    "accepted": ok, "runs": len(rows), "time_s": round(time.time() - t0, 2)})
        _write(sd / "status.json", {"state": "done", "pid": os.getpid(), "start": t0, "end": time.time(),
                                    "case": str(wf)})
    return rows
