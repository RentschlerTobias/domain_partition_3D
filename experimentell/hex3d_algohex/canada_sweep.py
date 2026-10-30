#!/usr/bin/env python3
"""P3a -- the n-sweep of AlgoHex replays on the canada domain.

The question the sweep answers: does the target cell count `-n` decide
whether the extracted hex complex can be turned into a valid block
structure? Four replays of one parametrisation state at n = 2000, 8000,
20000 and 100000, each ranked by bad valency and then invalid tets, and
the best two pushed through beam collapse -> TFI -> export.

Nothing here re-derives an AlgoHex command line. The replay, the beam
collapse, the TFI and the export are the commands of
`datagen/run_sample_keepbl.sh` with the canada paths filled in; this
module only collects their output into `report.json` and checks that the
report is complete.

usage:
  canada_sweep.py run [--sweep DIR] [--report PATH] [--ns 2000 8000 ...]
                       [--beam K] [--force]
  canada_sweep.py --check PATH

`run` is idempotent: a lane whose hex mesh and metrics already exist is
only re-read, never re-meshed, and a beam directory that already holds a
log is left alone unless `--force`. Nothing under /work/canada is ever
deleted -- the sweep writes into its own directory only.
"""

import argparse
import ast
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent

# --- the canada sweep (fixed by the task, not re-derived here) -----------
SWEEP = Path("/work/canada/sweep")
STATE = Path("/work/canada/dt/cad1/state/param_state.ovm")
TET = Path("/work/canada/dt/cad1/tet.vtk")
ALGOHEX = Path("/work/tools/algohex-size/build/Build/bin/HexMeshing")
NS = (2000, 8000, 20000, 100000)

# beam / tfi / export flags, from datagen/run_sample_keepbl.sh
BEAM_WIDTH, BEAM_DEPTH = 32, 20
BEAM_GUARD, BEAM_HTOL, BEAM_WORKERS = "struct", 0.15, 4

PY = sys.executable
ENV = dict(os.environ)
ENV.update(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
           NUMEXPR_NUM_THREADS="1", COLUMNS="200", PYTHONUNBUFFERED="1")

REPLAY_TIMEOUT = 10800  # s, the replay timeout of the sweep command


# --------------------------------------------------------------------------
# reading one replay
# --------------------------------------------------------------------------
def flat_numbers(obj, prefix=""):
    """Every numeric leaf of the metrics json, keyed by its dotted path.

    The AlgoHex json is a nested dict whose leaves are floats (or `null`
    for stages that did not run). Dumping all of them keeps the report
    honest: a later question about a stage nobody looked at can still be
    answered from the report alone.
    """
    out = {}
    if isinstance(obj, dict):
        for k, v in obj.items():
            out.update(flat_numbers(v, prefix + str(k) + "."))
    elif isinstance(obj, bool):
        out[prefix[:-1]] = int(obj)
    elif isinstance(obj, (int, float)):
        out[prefix[:-1]] = obj
    return out


def count_ovm_cells(path):
    """Number of hexahedra in an OVM ASCII file.

    AlgoHex prints no hex count in its log, so the count is read straight
    off the mesh: the `Polyhedra` section header is followed by its count
    line. Streamed, so a 90 MB state file costs nothing.
    """
    try:
        with open(path) as f:
            prev = None
            for line in f:
                s = line.strip()
                if prev == "Polyhedra":
                    return int(s)
                prev = s
    except (OSError, ValueError):
        pass
    return None


DEFECT_RE = re.compile(r"#invalid_tets\s*=\s*(\d+)\s+and\s+"
                       r"#invalid_edge_valencies\s*=\s*(\d+)")
TOTAL_RE = re.compile(r"^\s*Total time\s+(\d+) ms", re.M)
SAVE_RE = re.compile(r"save json data file into (\S+)")
CLOSED_RE = re.compile(r"^\s*Closed: (.+)$", re.M)
CELLS_RE = re.compile(r"Input Tet Mesh has (\d+) Vertices")


def read_log(path):
    """What the AlgoHex log adds on top of the metrics json."""
    out = {"log": str(path), "log_exists": Path(path).exists(),
           "log_invalid_tets": None, "log_invalid_edge_valencies": None,
           "total_time_ms": None, "finished": False, "closed_at": None,
           "input_tets": None, "log_bytes": 0}
    if not out["log_exists"]:
        return out
    text = Path(path).read_text(errors="replace")
    defects = DEFECT_RE.findall(text)
    if defects:
        inv, val = defects[-1]
        out["log_invalid_tets"] = int(inv)
        out["log_invalid_edge_valencies"] = int(val)
    m = TOTAL_RE.search(text)
    if m:
        out["total_time_ms"] = int(m.group(1))
    m = SAVE_RE.search(text)
    closed = CLOSED_RE.findall(text)
    # this build closes the run with a `Closed:` line instead of naming the
    # json it wrote, so either one means the replay is over; a lane that was
    # restarted appends, so the LAST closing line is the one that counts
    out["finished"] = bool(m) or bool(closed)
    if closed:
        out["closed_at"] = closed[-1].strip()
    m = CELLS_RE.search(text)
    if m:
        out["input_tets"] = int(m.group(1))
    out["log_bytes"] = len(text)
    return out


def _lane_file(d, n, stem, suffix):
    """`hexmeshing_n20000.log`, but also the older `..._n20k.log` spelling."""
    for name in (f"{stem}_n{n}{suffix}", f"{stem}_n{n // 1000}k{suffix}"):
        p = d / name
        if p.exists():
            return p
    hits = sorted(d.glob(f"{stem}_n{n}*{suffix}"))
    return hits[0] if hits else d / f"{stem}_n{n}{suffix}"


def read_run(d, n):
    """One replay lane: the metrics json, the log and the hex mesh."""
    row = {"n": n, "dir": str(d), "metrics_file": None, "metrics": {},
           "raw_metrics": None, "valid_volume": None,
           "valid_volume_seamless": None, "invalid_tets": None,
           "bad_valency": None, "bad_valency_seamless": None,
           "hexes": None, "ovm": None}
    mj = _lane_file(d, n, "metrics", ".json")
    metrics = {}
    raw = None
    if mj.exists():
        try:
            raw = json.loads(mj.read_text())
            metrics = flat_numbers(raw)
        except (OSError, ValueError) as e:
            row["metrics_error"] = str(e)
    row["metrics"] = metrics
    row["metrics_file"] = str(mj)
    # the headline numbers, and where they come from
    P = (raw or {}).get("Parametrization") or {}
    row["valid_volume"] = P.get("valid_volume")
    row["invalid_tets"] = P.get("n_invalid_param_tets")
    row["bad_valency"] = P.get("n_invalid_valencies")
    row["parametric_volume"] = P.get("parametric_volume")
    row["final_energy"] = P.get("final_energy")
    row["seamless"] = {
        "valid_volume": P.get("valid_volume_seamless"),
        "invalid_tets": P.get("n_invalid_param_tets_seamless"),
        "bad_valency": P.get("n_invalid_valencies_seamless"),
    }
    row.update(read_log(_lane_file(d, n, "hexmeshing", ".log")))
    ovm = _lane_file(d, n, "hex", ".ovm")
    row["ovm"] = str(ovm)
    row["hexes"] = count_ovm_cells(ovm) if ovm.exists() else None
    # a lane counts as complete when it produced a mesh and a filled-in
    # metrics json; the log's own "save json data file" line is only
    # informational, and some builds close the run without printing it
    row["complete"] = bool(
        row["valid_volume"] is not None and row["valid_volume"] > 0
        and row["hexes"] is not None)
    return row


def rank_key(row):
    """Rank by bad valency, then invalid tets (the task's order).

    A lane that is not complete (no mesh, or the metrics json still holds the
    -1 skeleton AlgoHex writes at startup) sorts last, so a running or failed
    lane can never look like the best one.
    """
    bv = row.get("bad_valency")
    it = row.get("invalid_tets")
    return (0 if row.get("complete") else 1,
            bv if isinstance(bv, (int, float)) else 1 << 30,
            it if isinstance(it, (int, float)) else 1 << 30,
            row.get("n", 0))


# --------------------------------------------------------------------------
# running one replay
# --------------------------------------------------------------------------
def replay(n, d, force=False):
    """Run the AlgoHex replay for one n, unless the lane is already done.

    Never writes over a lane that already holds a hex mesh unless `force`:
    the sweep only ever gets added to.
    """
    d = Path(d)
    d.mkdir(parents=True, exist_ok=True)
    ovm = d / f"hex_n{n}.ovm"
    if ovm.exists() and not force:
        print(f"[replay] n={n}: {ovm} exists, leaving it alone (--force to redo)")
        return read_run(d, n)
    row = read_run(d, n)
    if row["complete"] and not force:
        print(f"[replay] n={n}: already complete in {d}")
        return row
    log = d / f"hexmeshing_n{n}.log"
    cmd = [str(ALGOHEX), "--load-param-state", str(STATE),
           "-o", str(d / f"hex_n{n}.ovm"), "-j", str(d / f"metrics_n{n}.json"),
           "-n", str(n)]
    print(f"[replay] n={n}: {' '.join(cmd)}")
    with open(log, "a") as fh:
        fh.write(f"\n===== replay n={n} {time.strftime('%c')} =====\n")
        fh.flush()
        try:
            subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT,
                           env=ENV, timeout=REPLAY_TIMEOUT, check=False)
        except subprocess.TimeoutExpired:
            with open(log, "a") as fh:
                fh.write(f"\n[replay] TIMEOUT after {REPLAY_TIMEOUT} s\n")
            print(f"[replay] n={n}: TIMEOUT after {REPLAY_TIMEOUT} s")
    return read_run(d, n)


# --------------------------------------------------------------------------
# beam collapse -> tfi -> export
# --------------------------------------------------------------------------
BEAM_START_RE = re.compile(r"^\[beam\] start (\{.*\})", re.M)
BEAM_RESULT_RE = re.compile(r"^\[beam\] result: (\{.*\}) via sheets (\[.*\])", re.M)
BEAM_REJECT_RE = re.compile(r"^\[beam\] REJECT: (.*?) \((\d+) blocks, "
                            r"(\d+) cells\)", re.M)
# the reject message names the raw Hausdorff that broke the gate
REJECT_HD_RE = re.compile(r"Hausdorff ([\d.]+) > ([\d.]+)")
BEAM_DEPTH_RE = re.compile(r"^\[beam\] depth \d+: .*?best (\d+) blocks", re.M)
HAUS_TOL_RE = re.compile(r"-> Hausdorff tol ([\d.eE+-]+) \(flag ([\d.eE+-]+)\)")
HAUS_BAD_RE = re.compile(r"Hausdorff ([\d.eE+-]+) > ([\d.eE+-]+)")
HAUS_TO_SURF_RE = re.compile(r"boundary Hausdorff to input surface: ([\d.eE+-]+)")
BAF_RE = re.compile(r"^\s{2}(\w+)\s+(.*?)\s+->\s+(.*?)\s*$", re.M)
GUARD_RE = re.compile(r"^\[beam\] guard=(\w+): .*?width (\d+), depth (\d+)", re.M)
VALIDATOR_RE = re.compile(r"----- HexBlockValidator -----\n"
                          r"(.*?)(?=\n\n|\n(?![ !])|\Z)", re.S)


def _dict(text):
    try:
        return ast.literal_eval(text)
    except (ValueError, SyntaxError):
        return {}


def parse_beam_log(path):
    """The beam log, reduced to the numbers the report needs."""
    out = {"log": str(path), "log_bytes": Path(path).stat().st_size
           if Path(path).exists() else None}
    if out["log_bytes"] is None:
        return out
    out.update(status=None, reject_reason=None, reject_blocks=None,
               reject_cells=None, blocks=None, cells=None, inverted=None,
               min_sj=None, hausdorff=None, max_base_hausdorff=None,
               cuboids=None, valid=None, validator=[], before_after={},
               depths_done=0, blocks_best_so_far=None, hausdorff_tol=None,
               hausdorff_flag=None, hausdorff_to_surface=None, guard=None,
               width=None, depth=None)
    text = Path(path).read_text(errors="replace")
    for m in BEAM_START_RE.finditer(text):
        out["start"] = _dict(m.group(1))
    res = list(BEAM_RESULT_RE.finditer(text))
    if res:
        out["result"] = _dict(res[-1].group(1))
    rej = BEAM_REJECT_RE.search(text)
    if rej:
        out["status"] = "REJECT"
        out["reject_reason"] = rej.group(1).strip()
        out["reject_blocks"] = int(rej.group(2))
        out["reject_cells"] = int(rej.group(3))
        # the reject message names the raw Hausdorff that broke the gate
        m = REJECT_HD_RE.search(rej.group(1))
        if m:
            out["hausdorff"] = float(m.group(1))
            out["max_base_hausdorff"] = float(m.group(2))
    else:
        out["status"] = "OK"
    if out["status"] == "OK" and not BEAM_RESULT_RE.search(text):
        # a start line without a result or a reject: the lane is still going
        out["status"] = "RUNNING"
    bsf = BEAM_DEPTH_RE.findall(text)
    if bsf:
        out["blocks_best_so_far"] = max(int(x) for x in bsf)
    m = HAUS_TOL_RE.search(text)
    if m:
        out["hausdorff_tol"] = float(m.group(1))
        out["hausdorff_flag"] = float(m.group(2))
    m = HAUS_BAD_RE.search(text)
    if m:
        out["hausdorff_max"] = float(m.group(1))
        out["hausdorff_limit"] = float(m.group(2))
    m = HAUS_TO_SURF_RE.search(text)
    if m:
        out["hausdorff_to_surface"] = float(m.group(1))
    m = GUARD_RE.search(text)
    if m:
        out["guard"] = m.group(1)
        out["width"] = int(m.group(2))
        out["depth"] = int(m.group(3))
    # a beam that is still searching: the last depth line carries the best
    # block count so far
    best_so_far = BEAM_DEPTH_RE.findall(text)
    if best_so_far:
        out["blocks_best_so_far"] = int(best_so_far[-1])
        out["depths_done"] = len(best_so_far)
    baf = {}
    for m in BAF_RE.finditer(text):
        baf[m.group(1)] = (m.group(2), m.group(3))
    if baf:
        out["before_after"] = baf
    v = VALIDATOR_RE.search(text)
    if v:
        out["validator"] = [l for l in v.group(1).splitlines() if l.strip()]
        for line in out["validator"]:
            if line.strip() in ("VALID", "INVALID"):
                out["valid"] = line.strip() == "VALID"
    return out


def _sh(cmd, log, cwd=HERE):
    """Run one pipeline stage, tee-ing into `log`; return its returncode."""
    print(f"[run] {' '.join(str(c) for c in cmd)}")
    with open(log, "a") as fh:
        fh.write(f"\n===== {time.strftime('%c')} =====\n")
        fh.flush()
        return subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT,
                              env=ENV, cwd=str(cwd), check=False).returncode


def _pid_alive(pidfile):
    """True while the process named in `pidfile` is still running."""
    try:
        pid = int(Path(pidfile).read_text().strip())
    except (OSError, ValueError):
        return False
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def beam_lane(n, row, outdir, force=False):
    """beam collapse, then TFI, then export -- on one replay's hex mesh."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    b = {"n": n, "dir": str(outdir), "hex_ovm": row.get("ovm")}
    blog = outdir / "beam.log"
    blocks = outdir / "blocks.vtk"
    pidfile = outdir / "beam.pid"
    # a lane that another process is running right now is left alone
    running = pidfile.exists() and _pid_alive(pidfile)
    parsed = parse_beam_log(blog) if blog.exists() else {}
    if running and not force:
        print(f"[beam] n={n}: still running in {outdir}, not restarted")
        b["beam_running"] = True
    elif parsed.get("status") is not None and not force:
        if blocks.exists():
            print(f"[beam] n={n}: already done in {outdir}")
        else:
            # a rejected lane has nothing more to give: its log is the result
            print(f"[beam] n={n}: already attempted (rejected) in {outdir}")
    else:
        pidfile.write_text(str(os.getpid()))
        try:
            cmd = [PY, "beam_collapse.py", row["ovm"], "--input-vtk", str(TET),
                   "--out", str(blocks), "--width", str(BEAM_WIDTH),
                   "--depth", str(BEAM_DEPTH), "--guard", BEAM_GUARD,
                   "--hausdorff-tol", str(BEAM_HTOL),
                   "--workers", str(BEAM_WORKERS)]
            rc = _sh(cmd, blog)
            b["beam_returncode"] = rc
        finally:
            pidfile.unlink(missing_ok=True)
        parsed = parse_beam_log(blog)
    b.update(parsed)
    if b.get("beam_running"):
        b["status"] = "RUNNING"
    # the headline numbers, preferring the validator's own verdict; a rejected
    # lane has no start/result dict, so what the log did say is kept
    src = b.get("result") or b.get("start") or {}
    b["blocks"] = src.get("blocks") or b.get("reject_blocks") \
        or b.get("blocks_best_so_far")
    b["cells"] = src.get("cells") or b.get("reject_cells")
    for key, val in (("inverted_cells", src.get("inverted")),
                     ("min_scaled_jacobian", src.get("min_sj")),
                     ("hausdorff", src.get("hausdorff")),
                     ("hausdorff_to_surface", src.get("hausdorff")),
                     ("cuboids", src.get("cuboids"))):
        if val is not None:
            b[key] = val
    if "valid" not in b:
        b["valid"] = None
    if b["status"] == "REJECT":
        b["valid"] = False
    if b.get("beam_running"):
        b["tfi"] = {"skipped": "beam still running"}
        b["export"] = {"skipped": "beam still running"}
        return b
    if not blocks.exists():
        b["tfi"] = {"skipped": "no blocks.vtk (beam rejected)"}
        b["export"] = {"skipped": "no blocks.vtk (beam rejected)"}
        return b
    # TFI
    tfi = outdir / "tfi.vtk"
    tlog = outdir / "tfi.log"
    if tfi.exists() and not force:
        b["tfi"] = {"done": True, "out": str(tfi), "log": str(tlog)}
    else:
        rc = _sh([PY, "tfi.py", str(blocks), "--input-vtk", str(TET),
                  "--target-h", "0.05", "--require-lattices",
                  "--solve-divisions", "--apply-divisions", "--out", str(tfi)],
                 tlog)
        b["tfi"] = {"done": tfi.exists(), "returncode": rc, "out": str(tfi),
                    "log": str(tlog)}
    # export
    sample = outdir / "sample.npz"
    elog = outdir / "export.log"
    params = outdir / "params.json"
    cmd = [PY, "export_sample.py", str(blocks), "--tet", str(TET),
           "--target-h", "0.05", "--n", str(n), "--out", str(sample), "--check"]
    if params.exists():
        cmd += ["--params", str(params)]
    if sample.exists() and not force:
        b["export"] = {"done": True, "out": str(sample), "log": str(elog)}
    else:
        rc = _sh(cmd, elog)
        b["export"] = {"done": sample.exists(), "returncode": rc,
                       "out": str(sample), "log": str(elog)}
    # what the export log says
    if Path(elog).exists():
        text = Path(elog).read_text(errors="replace")
        for pat, key in ((r"\[export_sample\] (.*)", "export_lines"),):
            b["export"][key] = [m for m in re.findall(pat, text)][:12]
        b["export"]["check_ok"] = "OK" in text
    return b


def lane_dirs(sweep, n):
    """Every directory that holds a replay of this n.

    `n20000` is the sweep's own lane; a `n20000_*` sibling (a retry, a
    different timeout) is picked up as well, and the first one that is
    complete wins. Existing lanes are never rewritten.
    """
    sweep = Path(sweep)
    dirs = [sweep / f"n{n}"] + sorted(sweep.glob(f"n{n}_*"))
    return [d for d in dirs if d.is_dir()]


def best_lane(sweep, n, do_replay=False, force=False):
    """The most complete lane for this n, running the replay if none is."""
    rows = [read_run(d, n) for d in lane_dirs(sweep, n)]
    for row in rows:
        if row.get("complete"):
            return row
    if not do_replay:
        return rows[0] if rows else {"n": n}
    return replay(n, lane_dirs(sweep, n)[0], force=force)


# --------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------
def build_report(sweep=SWEEP, ns=NS, beam_k=2, force=False, do_replay=False,
                 beam_ns=None):
    """Collect the sweep into a report dict (and run what is missing)."""
    sweep = Path(sweep)
    runs = [best_lane(sweep, n, do_replay=do_replay, force=force) for n in ns]
    runs.sort(key=rank_key)
    for i, row in enumerate(runs, 1):
        row["rank"] = i
    best = [r["n"] for r in runs if r.get("complete")][:beam_k]
    # the lanes that get the beam stage: the best `beam_k` unless named
    targets = list(beam_ns) if beam_ns else best
    beam = []
    for n in targets:
        row = next((r for r in runs if r["n"] == n), {"n": n, "ovm": None})
        if not row.get("ovm"):
            continue
        b = beam_lane(n, row, sweep / f"beam_n{n}", force=force)
        b["lane_role"] = "best" if n in best else "diagnostic"
        beam.append(b)
    return {
        "schema": "canada-sweep/1",
        "sweep_dir": str(sweep),
        "state": str(STATE),
        "tet_vtk": str(TET),
        "ns": list(ns),
        "ranking": "bad valency (n_invalid_valencies), then invalid tets "
                   "(n_invalid_param_tets)",
        "created": time.strftime("%Y-%m-%d %H:%M:%S"),
        "runs": runs,
        "best": best,
        "beam_targets": targets,
        "beam": beam,
        "beam_incomplete": beam_incomplete(sweep, targets),
    }


def beam_incomplete(sweep, targets):
    """Beam lanes that were started but never produced a result.

    A depth-20 / width-32 beam search on a 100k-cell mesh costs ~18 min per
    depth level, so a lane can still be running when the report is written.
    """
    out = []
    for n in targets:
        d = Path(sweep) / f"beam_n{n}"
        log = d / "beam.log"
        if not log.exists():
            continue
        b = parse_beam_log(log)
        if b.get("status") in (None, "RUNNING"):
            out.append({"n": n, "dir": str(d), "log": str(log),
                        "blocks_best_so_far": b.get("blocks_best_so_far"),
                        "log_bytes": b.get("log_bytes")})
    return out


def check(report, verbose=True):
    """The report is usable only if it holds a full sweep and real beams."""
    problems = []
    runs = report.get("runs") or []
    beam = report.get("beam") or []
    with_metrics = [r for r in runs
                    if r.get("metrics") and r.get("valid_volume") is not None]
    with_beam = [b for b in beam if b.get("blocks") is not None]
    if len(with_metrics) < 4:
        problems.append(f"{len(with_metrics)}/4 n rows carry metrics "
                        f"(need 4: " +
                        ", ".join(f"n={r.get('n')}" for r in runs) + ")")
    if len(with_beam) < 2:
        problems.append(f"{len(with_beam)}/2 beam rows carry a block count")
    if verbose:
        print(f"canada sweep: {len(with_metrics)}/4 n rows with metrics, "
              f"{len(with_beam)}/2 beam rows")
        for r in runs:
            print(f"  n={str(r.get('n')):>7}  rank={r.get('rank')}  "
                  f"valid_volume={r.get('valid_volume')}  "
                  f"bad_valency={r.get('bad_valency')}  "
                  f"invalid_tets={r.get('invalid_tets')}  "
                  f"hexes={r.get('hexes')}  complete={r.get('complete')}")
        for b in beam:
            print(f"  beam n={str(b.get('n')):>7}  blocks={b.get('blocks')}  "
                  f"cells={b.get('cells')}  "
                  f"inverted={b.get('inverted_cells')}  "
                  f"min_sj={b.get('min_scaled_jacobian')}  "
                  f"Hausdorff={b.get('hausdorff')}  "
                  f"valid={b.get('valid')}  status={b.get('status')}")
    for p in problems:
        print(f"CHECK FAIL: {p}")
    return not problems


# --------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", default=None, metavar="REPORT",
                    help="validate REPORT and exit: it fails unless the report "
                         "holds >=4 n rows with metrics and >=2 beam rows")
    ap.add_argument("mode", nargs="?", default="run",
                    choices=["run", "check"],
                    help="run the sweep (default) or check a report")
    ap.add_argument("path", nargs="?", default=None,
                    help="report.json: written in run mode, read in check mode")
    ap.add_argument("--sweep", default=str(SWEEP))
    ap.add_argument("--ns", type=int, nargs="+", default=list(NS))
    ap.add_argument("--beam", type=int, default=2, dest="beam_k")
    ap.add_argument("--force", action="store_true",
                    help="re-run lanes that are already done")
    ap.add_argument("--no-replay", action="store_true",
                    help="only read what is already on disk; never start an "
                         "AlgoHex replay")
    ap.add_argument("--beam-ns", type=int, nargs="+", default=None,
                    metavar="N",
                    help="which n values get the beam lane (default: the best "
                         "--beam of them)")
    a = ap.parse_args(argv)

    if a.check is not None or a.mode == "check":
        path = Path(a.check or a.path or str(SWEEP / "report.json"))
        if not path.exists():
            print(f"CHECK FAIL: {path} does not exist")
            return 1
        try:
            report = json.loads(path.read_text())
        except ValueError as e:
            print(f"CHECK FAIL: {path} is not valid json: {e}")
            return 1
        return 0 if check(report) else 1

    report = build_report(sweep=Path(a.sweep), ns=tuple(a.ns),
                          beam_k=a.beam_k, force=a.force,
                          do_replay=not a.no_replay, beam_ns=a.beam_ns)
    out = Path(a.path or str(SWEEP / "report.json"))
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=1, sort_keys=False) + "\n")
    print(f"[report] wrote {out}")
    return 0 if check(report) else 1


if __name__ == "__main__":
    sys.exit(main())
