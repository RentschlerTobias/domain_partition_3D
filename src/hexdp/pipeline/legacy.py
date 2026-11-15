"""Run migrated command-line stages as subprocesses inside a plugin's workspace.

The plugins wrap the migrated scripts without changing their numerics: each stage runs as
``python -m <module> <args>`` with the job's working directory, CPU list (taskset), thread count and
timeout, and its output goes to a log file in the workspace. Stages that need pythonOCC run in the
dtOO environment (``HEXDP_OCC_PYTHON``, ``HEXDP_OCC_LD_LIBRARY_PATH``).
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import hexdp
from hexdp.contracts import Context

SRC = str(Path(hexdp.__file__).resolve().parents[1])
OCC_PYTHON = os.environ.get("HEXDP_OCC_PYTHON", "/opt/venv-dtoo/bin/python")
OCC_LD = os.environ.get("HEXDP_OCC_LD_LIBRARY_PATH", "/dtOO-install/lib")


class StageError(RuntimeError):
    pass


def run_module(ctx: Context, module: str, args, log: str, *, cwd=None, threads: int = 1, occ: bool = False,
               env: dict | None = None) -> Path:
    """Run ``python -m module args`` and return the log path; raise StageError with the log tail on failure."""
    python = OCC_PYTHON if occ else sys.executable
    cmd = [python, "-m", module, *map(str, args)]
    if ctx.cpu_list and shutil.which("taskset"):
        cmd = ["taskset", "-c", ctx.cpu_list, *cmd]
    e = dict(os.environ)
    e["PYTHONPATH"] = SRC + (os.pathsep + e["PYTHONPATH"] if e.get("PYTHONPATH") else "")
    e.update(PYTHONDONTWRITEBYTECODE="1", MPLBACKEND="Agg", OMP_NUM_THREADS=str(threads),
             OPENBLAS_NUM_THREADS=str(threads), MKL_NUM_THREADS=str(threads))
    if occ:
        e["LD_LIBRARY_PATH"] = OCC_LD + (os.pathsep + e["LD_LIBRARY_PATH"] if e.get("LD_LIBRARY_PATH") else "")
    e.update(env or {})
    logp = ctx.path(log)
    with open(logp, "w") as fh:
        try:
            r = subprocess.run(cmd, cwd=str(cwd or ctx.workspace), stdout=fh, stderr=subprocess.STDOUT, env=e,
                               timeout=ctx.timeout_s)
        except subprocess.TimeoutExpired as exc:
            raise StageError(f"{module}: timeout after {ctx.timeout_s} s (log {logp})") from exc
    if r.returncode != 0:
        tail = "".join(open(logp).readlines()[-15:])
        raise StageError(f"{module} failed with exit code {r.returncode} (log {logp}):\n{tail}")
    return logp
