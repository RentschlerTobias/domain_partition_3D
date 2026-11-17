"""Background jobs: start `hexdp run CASE` detached and look at its status (used by the MCP server)."""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

from hexdp.pipeline.case import load_case
from hexdp.pipeline.runner import status


def start(case_file, overrides: dict | None = None) -> Path:
    """Start the case in its own process (survives the caller); returns the workspace."""
    case = load_case(case_file, overrides)
    st = status(case.workspace)
    if st.get("state") == "running":
        raise RuntimeError(f"already running (pid {st.get('pid')}): {case.workspace}")
    case.workspace.mkdir(parents=True, exist_ok=True)
    with open(case.workspace / "run.out", "ab") as out:
        sets = [a for k, v in (overrides or {}).items() for a in ("--set", f"{k}={v}")]
        proc = subprocess.Popen([sys.executable, "-m", "hexdp.cli", "run", str(case.path), *sets], stdout=out,
                                stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, start_new_session=True)
    # running from this moment on, not only once the child has written its own status
    (case.workspace / "status.json").write_text(json.dumps(
        {"state": "running", "pid": proc.pid, "start": time.time(), "case": str(case.path)}))
    return case.workspace


def start_batch(workflow, runs: list[dict], job_dir) -> Path:
    """Start `hexdp batch` for the parameter sets in its own process; status and report in job_dir, the
    run workspaces below it."""
    job_dir = Path(job_dir)
    job_dir.mkdir(parents=True, exist_ok=True)
    (job_dir / "runs.json").write_text(json.dumps(runs))
    code = ("import json,sys; from hexdp.pipeline.workflow import run_batch; "
            "run_batch(sys.argv[1], json.load(open(sys.argv[2])), workspace_root=sys.argv[3], status_dir=sys.argv[3])")
    with open(job_dir / "run.out", "ab") as out:
        proc = subprocess.Popen([sys.executable, "-c", code, str(workflow), str(job_dir / "runs.json"), str(job_dir)],
                                stdout=out, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, start_new_session=True)
    (job_dir / "status.json").write_text(json.dumps(
        {"state": "running", "pid": proc.pid, "start": time.time(), "case": str(workflow)}))
    return job_dir
