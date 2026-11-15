"""Background jobs: start `hexdp run CASE` detached and look at its status (used by the MCP server)."""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

from hexdp.pipeline.case import load_case
from hexdp.pipeline.runner import status


def start(case_file) -> Path:
    """Start the case in its own process (survives the caller); returns the workspace."""
    case = load_case(case_file)
    st = status(case.workspace)
    if st.get("state") == "running":
        raise RuntimeError(f"already running (pid {st.get('pid')}): {case.workspace}")
    case.workspace.mkdir(parents=True, exist_ok=True)
    with open(case.workspace / "run.out", "ab") as out:
        proc = subprocess.Popen([sys.executable, "-m", "hexdp.cli", "run", str(case.path)], stdout=out,
                                stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, start_new_session=True)
    # running from this moment on, not only once the child has written its own status
    (case.workspace / "status.json").write_text(json.dumps(
        {"state": "running", "pid": proc.pid, "start": time.time(), "case": str(case.path)}))
    return case.workspace
