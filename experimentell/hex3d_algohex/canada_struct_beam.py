#!/usr/bin/env python3
"""T20 -- canada n=100000: the struct guard (inverted cells allowed) + TFI.

The question this lane answers: the abs guard of the T1 sweep refuses any
collapse that leaves one inverted cell, and on canada n=100000 that pins the
beam at the raw mesh (10 inverted cells after untangle, Hausdorff 0.1005) --
see /work/canada/sweep/beam_n100000. Here the structural guard (--guard
struct) is used instead: the Hausdorff to the input surface still binds
(0.15), but inverted cells no longer gate, so the beam may trade a few of
them for block structure. The lane then asks whether TFI can still build
lattices on the result and whether the deliverable exports.

The width-4/depth-3 beam of the task is measured too slow: one beam state
costs ~1095 s (T1 lane, same input), so depth 2 and 3 add 4 states each and
a full depth-3 lane would run ~10 ks -- past the 7200 s timeout and, with
the host restarting under long lanes, past what a background lane survives.
The lane therefore runs width 4 depth 1 (one state), records the curve point
and the ETA for the deeper levels, and spends the remaining time on the TFI
and export questions.

Usage:
    canada_struct_beam.py run [--depth N] [--timeout S] [--force]
    canada_struct_beam.py --check /work/canada/sweep/struct_beam/report.json

Everything lands in /work/canada/sweep/struct_beam/; the running T1 lane
(/work/canada/sweep/beam_n100000) is never touched -- this lane takes 4
cores from 0-19, T1 owns 20-23 and T5 owns 48-71.
"""
import argparse
import ast
import json
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

PY = "/opt/venv-sci/bin/python"
HERE = Path(__file__).resolve().parent
SWEEP = Path("/work/canada/sweep")
OUT = SWEEP / "struct_beam"
OVM = SWEEP / "n100000" / "hex_n100000.ovm"
TET = Path("/work/canada/dt/cad1/tet.vtk")
BLOCKS = OUT / "blocks.vtk"
N = 100000
CORES = "0,1,2,3"                     # 4 of the 0-19 cores this task may use
WIDTH = 4
GUARD = "struct"
SJ_FLOOR = 0.05
HAUS_TOL = 0.15
TIMEOUT = 7200

BASE_RE = re.compile(r"^\[beam\] start (\{.*\})$", re.M)
GUARD_RE = re.compile(r"^\[beam\] guard=(\w+): inverted <= (\d+), min sJ >= "
                      r"([-\d.]+), Hausdorff <= ([-\d.]+); width (\d+), "
                      r"depth (\d+)$", re.M)
DEPTH_RE = re.compile(r"^\[beam\] depth (\d+): (\d+) admissible, beam "
                      r"\[(.*?)\]  best (\d+) blocks \((\d+)s\)$", re.M)
STATE_RE = re.compile(r"(\d+)bl/ex(\d+)/sJ(-?\d+(?:\.\d+)?)")
RESULT_RE = re.compile(r"^\[beam\] result: (\{.*\}) via sheets \[(.*?)\]$", re.M)
BA_RE = re.compile(r"^  (\w+) +(.*?)  ->  (.*?)$", re.M)
HD_RE = re.compile(r"boundary Hausdorff to input surface: ([-\d.]+)")
EDGES_RE = re.compile(r"^\[clean_blocks\] block edges: (\d+) segments in "
                      r"(\d+) curves$", re.M)
UNTANGLE_RE = re.compile(r"^\[clean_blocks\] untangle: (\d+) vertices moved "
                         r"over (\d+) cells; inverted (\d+) -> (\d+), "
                         r"min scaled Jacobian ([-\d.]+) -> ([-\d.]+), "
                         r"boundary drift ([-\d.]+)$", re.M)
TFI_MISSING_RE = re.compile(r"^\[tfi\] (\d+) blocks are not lattices", re.M)
TFI_REJECT_RE = re.compile(r"^\[tfi\] REJECTED: (\d+) block\(s\) have no "
                           r"lattice", re.M)
TFI_HEAD_RE = re.compile(r"^\[tfi\] .*: \d+ cells, \d+ blocks$", re.M)
TFI_REFILL_RE = re.compile(r"^\[tfi\] refilled mesh: scaled Jacobian min "
                           r"([-\d.]+), mean ([-\d.]+), (\d+) inverted$", re.M)
TFI_CELLS_RE = re.compile(r"^\[tfi\] cells in the lattice blocks: (\d+) -> "
                          r"(\d+)$", re.M)
TFI_REBUILT_RE = re.compile(r"^\[tfi\] (\d+) blocks rebuilt, (\d+) not a "
                            r"lattice$", re.M)
TFI_MESH_RE = re.compile(r"^\[tfi\] whole mesh scaled Jacobian: min "
                         r"([-\d.]+) -> ([-\d.]+), (\d+) -> (\d+) inverted$",
                         re.M)
TFI_CONFORM_RE = re.compile(r"^\[tfi\] conformity of the existing divisions: "
                            r"(\d+) classes disagree$", re.M)
TFI_ASSERT_RE = re.compile(r"^AssertionError: (.+)$", re.M)
TFI_CLASSES_RE = re.compile(r"^\[tfi\] (\d+) direction classes over (\d+) block "
                            r"axes \((\d+) links\)", re.M)
EXPORT_WROTE_RE = re.compile(r"^\[export_sample\] wrote (\S+)", re.M)
EXPORT_ROW_RE = re.compile(r"^  vertices (\d+), blocks (\d+), shell (\d+), "
                           r"directed edges (\d+)$", re.M)


def _num(s):
    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            return s                      # e.g. the validator's status dict


def _sh(cmd, log, timeout=None):
    """Run cmd, tee its output into log, return (returncode, wall seconds)."""
    t0 = time.time()
    with open(log, "w") as fh:
        fh.write("+ " + " ".join(cmd) + "\n\n")
        fh.flush()
        try:
            rc = subprocess.call(cmd, stdout=fh, stderr=subprocess.STDOUT,
                                 timeout=timeout)
        except subprocess.TimeoutExpired:
            rc = "timeout"
    return rc, time.time() - t0


def run_beam(depth, timeout, force=False):
    """The width-4 beam with --guard struct, pinned to 4 of the 0-19 cores."""
    log = OUT / f"beam_d{depth}.log"
    if BLOCKS.exists() and log.exists() and not force:
        print(f"[run] beam depth {depth}: reusing {log.name} and blocks.vtk")
        return "reused", None, log
    cmd = ["taskset", "-c", CORES, PY, "-u", str(HERE / "beam_collapse.py"),
           str(OVM), "--input-vtk", str(TET), "--out", str(BLOCKS),
           "--width", str(WIDTH), "--depth", str(depth), "--guard", GUARD,
           "--hausdorff-tol", str(HAUS_TOL), "--workers", "4"]
    print(f"[run] beam: {' '.join(cmd)}")
    rc, wall = _sh(cmd, log, timeout=timeout)
    print(f"[run] beam rc={rc} wall={wall:.0f}s log={log}")
    return rc, wall, log


def parse_beam_log(path):
    """Everything the beam log says, as plain data."""
    text = Path(path).read_text()
    out = {"log": str(path)}
    m = BASE_RE.search(text)
    out["base"] = ast.literal_eval(m.group(1)) if m else None
    m = GUARD_RE.search(text)
    if m:
        out["guard"] = {"kind": m.group(1), "inverted_max": int(m.group(2)),
                        "sj_floor": float(m.group(3)),
                        "hausdorff_tol": float(m.group(4)),
                        "width": int(m.group(5)), "depth": int(m.group(6))}
    depths = []
    for m in DEPTH_RE.finditer(text):
        depths.append({
            "depth": int(m.group(1)), "admissible": int(m.group(2)),
            "beam": [{"blocks": int(a), "excess": int(b), "sj": float(c)}
                     for a, b, c in STATE_RE.findall(m.group(3))],
            "best_blocks": int(m.group(4)), "wall_s": int(m.group(5))})
    out["depths"] = depths
    m = RESULT_RE.search(text)
    if m:
        out["result"] = ast.literal_eval(m.group(1))
        out["sheets"] = [int(x) for x in
                         re.findall(r"\d+", m.group(2))]
    ba = {}
    for m in BA_RE.finditer(text):
        k, a, b = m.group(1), m.group(2), m.group(3)
        if k in ("blocks", "cuboids", "inverted_cells",
                 "min_scaled_jacobian", "mean_scaled_jacobian"):
            ba[k] = [_num(a), _num(b)]
    out["before_after"] = ba
    # untangle rounds: the raw mesh starts with 148 inverted cells
    rounds = [{"vertices": int(m.group(1)), "cells": int(m.group(2)),
               "inverted": [int(m.group(3)), int(m.group(4))],
               "sj": [float(m.group(5)), float(m.group(6))],
               "drift": float(m.group(7))}
              for m in UNTANGLE_RE.finditer(text)]
    out["untangle"] = rounds
    if rounds:
        out["pre_untangle_inverted"] = rounds[0]["inverted"][0]
    v = {"status": None, "failures": [], "quality": {}, "faces": {}}
    i = text.find("----- HexBlockValidator -----")
    if i >= 0:
        block = text[i:].split("\n")
        block = block[1:]
        for ln in block:
            if ln.startswith("  ! "):
                v["failures"].append(ln[4:])
            elif ln.startswith("  block face counts: "):
                v["faces"] = ast.literal_eval(ln.split(": ", 1)[1])
            elif re.match(r"^  [a-z_]+: ", ln):
                k, val = ln[2:].split(": ", 1)
                v["quality"][k] = _num(val)
            elif ln in ("VALID", "INVALID"):
                v["status"] = ln
            elif ln.strip() and not ln.startswith("["):
                break
    out["validator"] = v
    # "is it a valid subdivision?" -- the inverted cells are a cell-quality
    # failure, not a structural one, so they do not count here.
    out["structurally_valid"] = not [f for f in v.get("failures", [])
                                     if "inverted cells" not in f]
    m = HD_RE.search(text)
    out["hausdorff_final"] = float(m.group(1)) if m else None
    m = EDGES_RE.search(text)
    out["block_edges"] = ({"segments": int(m.group(1)),
                           "curves": int(m.group(2))} if m else None)
    out["status"] = "OK" if (out.get("result") and out["validator"]["status"]
                             ) else ("REJECT" if "[beam] REJECT" in text
                                     else "incomplete")
    return out


def run_tfi(force=False):
    """Can TFI build lattices on the struct-guard structure? And refill it?"""
    log = OUT / "tfi.log"
    outv = OUT / "tfi_blocks.vtk"
    if log.exists() and not force:
        print(f"[run] tfi: reusing {log.name}")
        return parse_tfi_log(log)
    if not BLOCKS.exists():
        return {"status": "skipped", "why": "no blocks.vtk (beam did not "
                                            "finish)", "log": str(log)}
    cmd = [PY, "-u", str(HERE / "tfi.py"), str(BLOCKS), "--input-vtk",
           str(TET), "--target-h", "0.05", "--solve-divisions",
           "--apply-divisions", "--require-lattices", "--out", str(outv)]
    print(f"[run] tfi: {' '.join(cmd)}")
    rc, wall = _sh(cmd, log, timeout=TIMEOUT)
    print(f"[run] tfi rc={rc} wall={wall:.0f}s log={log}")
    return parse_tfi_log(log, rc, wall)


def parse_tfi_log(log, rc=None, wall=None):
    text = Path(log).read_text() if Path(log).exists() else ""
    out = {"status": "skipped", "log": str(log)}
    if rc is not None:
        out["returncode"] = rc
    if wall is not None:
        out["wall_s"] = round(wall, 1)
    m = TFI_MISSING_RE.search(text)
    if m:
        out["missing"] = int(m.group(1))
    else:
        r = TFI_REJECT_RE.search(text)
        if r:
            out["missing"] = int(r.group(1))
        elif TFI_HEAD_RE.search(text):
            out["missing"] = 0             # the stage ran and said nothing
        else:
            out["missing"] = None          # never got that far
    out["lattices_built"] = (out["missing"] == 0)
    m = TFI_CELLS_RE.search(text)
    if m:
        out["cells_now"] = int(m.group(1))
        out["cells_after"] = int(m.group(2))
    m = TFI_REFILL_RE.search(text)
    if m:
        out["refill"] = {"sj_min": float(m.group(1)),
                         "sj_mean": float(m.group(2)),
                         "inverted": int(m.group(3))}
    # The division solve can die before the --require-lattices gate; say so.
    m = TFI_ASSERT_RE.search(text)
    if m:
        out["error"] = m.group(1)
        out["failed_stage"] = "division_solve"
    m = TFI_CONFORM_RE.search(text)
    if m:
        out["conformity_classes_disagree"] = int(m.group(1))
    m = TFI_CLASSES_RE.search(text)
    if m:
        out["direction_classes"] = int(m.group(1))
        out["block_axes"] = int(m.group(2))
        out["links"] = int(m.group(3))
    if rc is None:
        # A reused log carries no exit code, but a crash in it is a crash.
        if out.get("failed_stage"):
            out["status"] = "failed"
        return out
    if rc == 0:
        out["status"] = "ok"
    elif rc == 3:
        out["status"] = "rejected"        # no lattice: TFI cannot refill
    elif rc == 2:
        out["status"] = "not_watertight"
    else:
        out["status"] = "failed"
    if "REFUSING to write" in text:
        out["status"] = "not_watertight"
    return out


def run_tfi_rebuild(force=False):
    """Refill the structure with TFI at its *existing* divisions.

    The --solve-divisions path can die before the --require-lattices gate,
    so this second run answers the other half of the question: for the
    blocks that do carry a lattice, does the refill hold up?
    """
    log = OUT / "tfi_rebuild.log"
    if log.exists() and not force:
        print(f"[run] tfi_rebuild: reusing {log.name}")
        return parse_tfi_rebuild_log(log)
    if not BLOCKS.exists():
        return {"status": "skipped", "why": "no blocks.vtk"}
    cmd = [PY, "-u", str(HERE / "tfi.py"), str(BLOCKS), "--input-vtk",
           str(TET), "--target-h", "0.05"]
    print(f"[run] tfi_rebuild: {' '.join(cmd)}")
    rc, wall = _sh(cmd, log, timeout=TIMEOUT)
    print(f"[run] tfi_rebuild rc={rc} wall={wall:.0f}s log={log}")
    return parse_tfi_rebuild_log(log, rc, wall)


def parse_tfi_rebuild_log(log, rc=None, wall=None):
    """Read the plain-TFI rebuild log: how many blocks rebuilt, what happened
    to the mesh quality."""
    text = Path(log).read_text() if Path(log).exists() else ""
    out = {"status": "skipped", "log": str(log)}
    if rc is not None:
        out["returncode"] = rc
    if wall is not None:
        out["wall_s"] = round(wall, 1)
    m = TFI_REBUILT_RE.search(text)
    if m:
        out["rebuilt"] = int(m.group(1))
        out["not_a_lattice"] = int(m.group(2))
    m = TFI_MESH_RE.search(text)
    if m:
        out["mesh"] = {"sj_min_before": float(m.group(1)),
                       "sj_min_after": float(m.group(2)),
                       "inverted_before": int(m.group(3)),
                       "inverted_after": int(m.group(4))}
    # the seven-metric table: "name  before (n bad)  after (n bad)"
    metrics = {}
    for line in text.splitlines():
        f = line.split()
        if len(f) == 5 and f[0] != "metric" and f[2].startswith("(") \
                and f[4].startswith("("):
            try:
                metrics[f[0]] = {"before": float(f[1]),
                                 "before_bad": int(f[2][1:].split()[0]),
                                 "after": float(f[3]),
                                 "after_bad": int(f[4][1:].split()[0])}
            except ValueError:
                pass
    if metrics:
        out["metrics"] = metrics
    if rc == 0:
        out["status"] = "ok"
    elif rc is not None:
        out["status"] = "failed"
    elif out.get("rebuilt") is not None:
        # a reused log carries no exit code, but it did run
        out["status"] = "ok"
    return out


def run_export(force=False):
    """Export the deliverable with export_sample.py --check, as the sweep does."""
    log = OUT / "export.log"
    out = OUT / f"sample_{N}.npz"
    if log.exists() and out.exists() and not force:
        print(f"[run] export: reusing {log.name}")
        return parse_export_log(log, None, None)
    if not BLOCKS.exists():
        return {"status": "skipped", "why": "no blocks.vtk", "log": str(log)}
    cmd = [PY, "-u", str(HERE / "export_sample.py"), str(BLOCKS), "--tet",
           str(TET), "--target-h", "0.05", "--n", str(N), "--out", str(out),
           "--check"]
    print(f"[run] export: {' '.join(cmd)}")
    rc, wall = _sh(cmd, log, timeout=TIMEOUT)
    print(f"[run] export rc={rc} wall={wall:.0f}s log={log}")
    return parse_export_log(log, rc, wall)


def parse_export_log(log, rc=None, wall=None):
    text = Path(log).read_text() if Path(log).exists() else ""
    out = {"status": "skipped", "log": str(log)}
    if rc is not None:
        out["returncode"] = rc
    if wall is not None:
        out["wall_s"] = round(wall, 1)
    m = EXPORT_WROTE_RE.search(text)
    if m:
        out["out"] = m.group(1)
        out["exists"] = Path(m.group(1)).exists()
    m = EXPORT_ROW_RE.search(text)
    if m:
        out["sample"] = {"vertices": int(m.group(1)),
                         "blocks": int(m.group(2)),
                         "shell": int(m.group(3)),
                         "edges": int(m.group(4))}
    if rc is None:
        out["status"] = "done" if out.get("exists") else "skipped"
        return out
    out["status"] = "done" if (rc == 0 and out.get("exists")) else "failed"
    return out


def measure_blocks_vtk():
    """Independent measurement of the written deliverable, for the check."""
    import numpy as np
    import meshio
    if not BLOCKS.exists():
        return None
    m = meshio.read(BLOCKS)
    hexes = [c for c in m.cells if c.type == "hexahedron"]
    bid = [d for b, d in zip(m.cells, m.cell_data["block_id"])
           if b.type == "hexahedron"]
    H = np.vstack([c.data for c in hexes])
    B = np.concatenate([np.asarray(d).ravel() for d in bid])
    q = Path(str(BLOCKS).replace(".vtk", "_quality.vtk"))
    meas = {"cells": int(len(H)), "blocks": int(B.max()) + 1,
            "points": int(len(m.points))}
    if q.exists():
        mq = meshio.read(q)
        sq = [d for b, d in zip(mq.cells,
                                mq.cell_data["scaled_jacobian_x1000"])
              if b.type == "hexahedron"]
        sj = np.concatenate([np.asarray(d).ravel() for d in sq]) / 1000.0
        meas["inverted"] = int((sj <= 0).sum())
        meas["min_sj"] = round(float(sj.min()), 4)
    return meas


def build_report(beam_rc, beam_wall, beam_log, depth, tfi, export,
                 eta, note, rebuild=None):
    beam = parse_beam_log(beam_log)
    beam.update({"guard_kind": GUARD, "width": WIDTH, "depth": depth,
                 "workers": 4, "cores": CORES, "hausdorff_tol": HAUS_TOL,
                 "sj_floor": SJ_FLOOR, "timeout_s": TIMEOUT,
                 "returncode": beam_rc,
                 "wall_s": round(beam_wall, 1) if beam_wall else None})
    meas = measure_blocks_vtk()
    res = beam.get("result") or {}
    if meas:
        beam["measured_blocks_vtk"] = meas
        # the report's numbers must agree with what is on disk
        res = dict(res)
        res["cells"] = meas["cells"]
        res["blocks"] = meas["blocks"]
        if "inverted" in meas:
            res["inverted"] = meas["inverted"]
        if "min_sj" in meas:
            res["min_sj"] = meas["min_sj"]
        beam["result"] = res
    base = beam.get("base") or {}
    curve = [{"depth": d["depth"], "admissible": d["admissible"],
              "best_blocks": d["best_blocks"], "wall_s": d["wall_s"],
              "states": d["beam"]} for d in beam.get("depths", [])]
    report = {
        "task": "T20",
        "continued_by": "T40",
        "title": "canada n=100000: struct guard (inverted cells allowed) "
                 "+ TFI attempt",
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "input": {"ovm": str(OVM), "tet": str(TET),
                  "cells": base.get("cells"),
                  "blocks_base": base.get("blocks"),
                  "inverted_base": base.get("inverted"),
                  "min_sj_base": base.get("min_sj"),
                  "hausdorff_base": base.get("hausdorff"),
                  "pre_untangle_inverted": beam.get("pre_untangle_inverted")},
        "beam": beam,
        "curve": curve,
        "tfi": tfi,
        "tfi_rebuild": rebuild or {"status": "skipped"},
        "export": export,
        "eta": eta,
        "note": note,
    }
    return report


def verdict(report):
    """The answer the lane exists for, as one sentence."""
    beam = report["beam"]
    res = beam.get("result") or {}
    v = beam.get("validator") or {}
    fails = [f for f in v.get("failures", []) if "inverted cells" not in f]
    inv = res.get("inverted")
    tfi = report["tfi"]
    tfi_bits = []
    if tfi.get("failed_stage") == "division_solve":
        tfi_bits.append(
            f"TFI cannot build the lattices: {tfi.get('missing')} of "
            f"{res.get('blocks')} blocks have no lattice and the division "
            f"solve dies first ({tfi.get('error')})")
    elif tfi.get("status") == "ok":
        tfi_bits.append("TFI builds lattices on every block and refills "
                        "watertight")
    elif tfi.get("status") == "rejected":
        tfi_bits.append(f"TFI cannot refill: {tfi.get('missing')} blocks "
                        "have no lattice")
    elif tfi.get("status") == "skipped":
        tfi_bits.append(f"TFI skipped ({tfi.get('why')})")
    else:
        tfi_bits.append(f"TFI {tfi.get('status')}")
    exp = report["export"]
    if exp.get("status") == "done":
        tfi_bits.append("the deliverable exports and round-trips")
    elif exp.get("status") != "skipped":
        tfi_bits.append("the export fails")
    rb = report.get("tfi_rebuild") or {}
    if rb.get("rebuilt") is not None:
        m = rb.get("mesh") or {}
        tfi_bits.append(
            f"refilling the {rb.get('rebuilt')} blocks that do carry a "
            f"lattice takes the inverted cells from "
            f"{m.get('inverted_before')} to {m.get('inverted_after')}")
    return (f"depth {beam.get('depth')} keeps {res.get('blocks')} blocks "
            f"with {inv} inverted cells (Hausdorff {beam.get('hausdorff_final')}, "
            f"min sJ {res.get('min_sj')}); structurally "
            f"{'valid' if not fails else 'INVALID (' + '; '.join(fails) + ')'} "
            f"as a subdivision" +
            (f" but the validator says {v.get('status')} because inverted "
             "cells are a hard failure" if v.get("status") == "INVALID" and
             not fails else "") +
            "; " + ", ".join(tfi_bits) + ".")


# ------------------------------------------------------------------ run mode
def run(depth, timeout, force):
    OUT.mkdir(parents=True, exist_ok=True)
    note = []
    if depth < 3:
        note.append(
            f"width {WIDTH}/depth 3 was not run: one beam state costs "
            f"~1095 s on this input (T1 lane, /work/canada/sweep/"
            f"beam_n100000/beam.log) and the 7200 s timeout covers only "
            f"depth 1-2; the lane stopped at depth {depth} and reports the "
            f"ETA below. The first depth-3 attempt (beam.log) died in the "
            f"pre-collapse untangle when the host restarted, not in the "
            f"beam itself.")
    beam_rc, beam_wall, beam_log = run_beam(depth, timeout, force)
    tfi = run_tfi(force)
    rebuild = run_tfi_rebuild(force)
    export = (run_export(force) if tfi.get("lattices_built")
              else {"status": "skipped",
                    "why": "TFI cannot build the lattices",
                    "log": str(OUT / "export.log")})
    # ETA for the depths the task asked for, from the measured state cost
    d1 = next((d for d in parse_beam_log(beam_log)["depths"]
               if d["depth"] == 1), None)
    per_state = d1["wall_s"] if d1 else None
    eta = {"per_state_s": per_state, "states_per_level": WIDTH,
           "note": "one _children() call per beam state, states expanded "
                   "one level at a time; a collapsed state has fewer "
                   "sheets than the 330-block base, so this is an upper "
                   "bound"}
    if per_state:
        eta["depth2_s"] = per_state * (1 + WIDTH)
        eta["depth3_s"] = per_state * (1 + WIDTH + WIDTH)
        eta["depth2_within_timeout"] = eta["depth2_s"] + 600 < TIMEOUT
        eta["depth3_within_timeout"] = eta["depth3_s"] + 600 < TIMEOUT
    report = build_report(beam_rc, beam_wall, beam_log, depth, tfi, export,
                          eta, " ".join(note), rebuild)
    report["verdict"] = verdict(report)
    rp = OUT / "report.json"
    rp.write_text(json.dumps(report, indent=1) + "\n")
    print(f"\n[T20] wrote {rp}")
    print(json.dumps({k: report[k] for k in ("input", "curve", "tfi",
                                             "export", "eta", "verdict")},
                     indent=1))
    return 0 if report["beam"].get("result") else 1


# ---------------------------------------------------------------- check mode
def check(report_path):
    """Verify report.json against the logs and the files on disk."""
    bad = []
    rp = Path(report_path)
    if not rp.exists():
        print(f"BAD: {rp} does not exist")
        return 1
    try:
        rep = json.loads(rp.read_text())
    except Exception as e:
        print(f"BAD: {rp} is not JSON: {e}")
        return 1
    if rep.get("task") != "T20":
        bad.append("report.task is not T20")

    beam = rep.get("beam") or {}
    if beam.get("guard_kind") != GUARD:
        bad.append(f"beam.guard_kind is {beam.get('guard_kind')}, want struct")
    if beam.get("width") != WIDTH:
        bad.append(f"beam.width is {beam.get('width')}, want {WIDTH}")
    if beam.get("hausdorff_tol") != HAUS_TOL:
        bad.append(f"beam.hausdorff_tol is {beam.get('hausdorff_tol')}, "
                   f"want {HAUS_TOL}")
    if beam.get("workers") != 4:
        bad.append("beam.workers is not 4")
    if not re.fullmatch(r"[0-9]+(,[0-9]+)*", str(beam.get("cores", ""))):
        bad.append(f"beam.cores {beam.get('cores')} is not a core list")
    else:
        cpus = [int(c) for c in str(beam["cores"]).split(",")]
        if len(cpus) != 4 or min(cpus) < 0 or max(cpus) > 19:
            bad.append(f"beam.cores {beam['cores']} is not 4 cores from 0-19")
    depth = beam.get("depth")
    if not isinstance(depth, int) or depth < 1 or depth > 3:
        bad.append(f"beam.depth {depth} is not 1..3")
    if beam.get("returncode") not in (0, "reused", "timeout", None):
        bad.append(f"beam.returncode is {beam.get('returncode')}")

    inp = rep.get("input") or {}
    if inp.get("cells") != 105837:
        bad.append(f"input.cells is {inp.get('cells')}, want 105837")
    if not Path(inp.get("ovm", "")).exists():
        bad.append(f"input.ovm {inp.get('ovm')} does not exist")
    if not Path(inp.get("tet", "")).exists():
        bad.append(f"input.tet {inp.get('tet')} does not exist")

    # every number in the report must match the beam log it came from
    log = beam.get("log")
    if not log or not Path(log).exists():
        bad.append(f"beam.log {log} does not exist")
    else:
        parsed = parse_beam_log(log)
        for k in ("cells", "blocks", "inverted", "min_sj", "hausdorff"):
            if (beam.get("base") or {}).get(k) != (parsed.get("base")
                                                   or {}).get(k):
                bad.append(f"base.{k} disagrees with {log}")
        if beam.get("hausdorff_final") != parsed.get("hausdorff_final"):
            bad.append("hausdorff_final disagrees with the beam log")
        if (beam.get("validator") or {}).get("status") != \
                (parsed.get("validator") or {}).get("status"):
            bad.append("validator.status disagrees with the beam log")
        if beam.get("structurally_valid") != parsed.get("structurally_valid"):
            bad.append("structurally_valid disagrees with the beam log")
        for got, want in zip(beam.get("depths", []),
                             parsed.get("depths", [])):
            for k in ("depth", "admissible", "best_blocks", "wall_s"):
                if got.get(k) != want.get(k):
                    bad.append(f"depth line {got.get('depth')}.{k} "
                               f"disagrees with the log")
        if len(beam.get("depths", [])) != len(parsed.get("depths", [])):
            bad.append("the report's depth curve is not the log's")
        res = beam.get("result") or {}
        if res.get("blocks") is None or res.get("blocks") <= 0:
            bad.append("beam.result.blocks is missing")
        if beam.get("depth") == 1 and res.get("blocks") is not None:
            base_blocks = (beam.get("base") or {}).get("blocks")
            if res["blocks"] >= base_blocks:
                bad.append(f"depth 1 did not collapse: {res['blocks']} >= "
                           f"{base_blocks} blocks")

    # the deliverable itself
    if not BLOCKS.exists():
        bad.append(f"{BLOCKS} was not written")
    else:
        meas = measure_blocks_vtk()
        res = beam.get("result") or {}
        for k in ("cells", "blocks"):
            if meas.get(k) != res.get(k):
                bad.append(f"blocks.vtk has {meas.get(k)} {k}, report says "
                           f"{res.get(k)}")
        if "inverted" in meas and meas["inverted"] != res.get("inverted"):
            bad.append(f"blocks.vtk has {meas['inverted']} inverted cells, "
                       f"report says {res.get('inverted')}")
        if "min_sj" in meas and abs(meas["min_sj"] - res.get("min_sj",
                                                             1e9)) > 0.002:
            bad.append(f"blocks.vtk min sJ {meas['min_sj']} != report "
                       f"{res.get('min_sj')}")
        if (beam.get("result") or {}).get("hausdorff") is not None and \
                beam.get("hausdorff_final") is not None:
            if beam["result"]["hausdorff"] > HAUS_TOL + 1e-9:
                bad.append(f"result hausdorff {beam['result']['hausdorff']} "
                           f"exceeds the {HAUS_TOL} guard")
        if (beam.get("validator") or {}).get("faces"):
            n_faces = sum((beam["validator"]["faces"] or {}).values())
            if n_faces != meas["blocks"]:
                bad.append(f"validator sees {n_faces} blocks, blocks.vtk "
                           f"has {meas['blocks']}")

    v = beam.get("validator") or {}
    if v.get("status") not in ("VALID", "INVALID"):
        bad.append("no validator status in the report")
    fails = [f for f in v.get("failures", []) if "inverted cells" not in f]
    if fails:
        bad.append(f"structural hard failures: {fails}")
    rep.setdefault("verdict", "")
    if not rep["verdict"] or "blocks" not in rep["verdict"]:
        bad.append("verdict is empty or does not state the block count")

    tfi = rep.get("tfi") or {}
    tlog = tfi.get("log")
    if tfi.get("status") == "skipped":
        if not tfi.get("why"):
            bad.append("tfi skipped without a why")
    else:
        if not tlog or not Path(tlog).exists():
            bad.append(f"tfi.log {tlog} does not exist")
        else:
            parsed = parse_tfi_log(tlog)
            if parsed.get("missing") != tfi.get("missing"):
                bad.append(f"tfi.missing {tfi.get('missing')} disagrees with "
                           f"the tfi log ({parsed.get('missing')})")
            if parsed.get("lattices_built") != tfi.get("lattices_built"):
                bad.append("tfi.lattices_built disagrees with the tfi log")
            if tfi.get("lattices_built"):
                for k in ("cells_now", "cells_after"):
                    if parsed.get(k) is not None and parsed[k] != tfi.get(k):
                        bad.append(f"tfi.{k} {tfi.get(k)} != log {parsed[k]}")
                if parsed.get("refill") and parsed["refill"] != tfi.get("refill"):
                    bad.append("tfi.refill disagrees with the tfi log")
        if tfi.get("status") == "ok":
            if tfi.get("missing") != 0:
                bad.append("tfi ok but blocks have no lattice")
            if not (OUT / "tfi_blocks.vtk").exists():
                bad.append("tfi ok but tfi_blocks.vtk is missing")

    rb = rep.get("tfi_rebuild") or {}
    rlog = rb.get("log")
    if rb.get("status") == "skipped":
        pass
    elif not rlog or not Path(rlog).exists():
        bad.append(f"tfi_rebuild.log {rlog} does not exist")
    else:
        prb = parse_tfi_rebuild_log(rlog)
        for k in ("rebuilt", "not_a_lattice"):
            if prb.get(k) is not None and prb[k] != rb.get(k):
                bad.append(f"tfi_rebuild.{k} {rb.get(k)} != log {prb[k]}")
        if prb.get("mesh") and prb["mesh"] != rb.get("mesh"):
            bad.append("tfi_rebuild.mesh disagrees with the log")

    exp = rep.get("export") or {}
    if tfi.get("lattices_built"):
        if exp.get("status") != "done":
            bad.append(f"export is {exp.get('status')}, want done")
        else:
            if not Path(exp.get("out", "")).exists():
                bad.append(f"export out {exp.get('out')} does not exist")
            if exp.get("returncode") != 0:
                bad.append(f"export returncode {exp.get('returncode')}")
            elog = exp.get("log")
            if not elog or not Path(elog).exists():
                bad.append(f"export.log {elog} does not exist")
    else:
        if exp.get("status") != "skipped":
            bad.append(f"export is {exp.get('status')} while TFI cannot "
                       f"build lattices")

    eta = rep.get("eta") or {}
    for k in ("per_state_s", "depth2_s", "depth3_s"):
        if not isinstance(eta.get(k), (int, float)) or eta[k] <= 0:
            bad.append(f"eta.{k} is {eta.get(k)}")
    if eta.get("depth3_s", 0) <= eta.get("depth2_s", 0):
        bad.append("eta depth3_s is not above depth2_s")

    curve = rep.get("curve") or []
    if not curve:
        bad.append("no blocks-vs-depth curve in the report")
    else:
        for c in curve:
            for k in ("depth", "admissible", "best_blocks", "wall_s"):
                if k not in c:
                    bad.append(f"curve point missing {k}")
        if any(c.get("best_blocks") is None for c in curve):
            bad.append("a curve point has no block count")

    for f in bad:
        print(f"BAD: {f}")
    print(f"\n{len(bad)} problem(s)" if bad else "\nreport OK")
    return 1 if bad else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("command", nargs="?", default=None, choices=("run",))
    ap.add_argument("--check", metavar="REPORT", default=None,
                    help="verify a report.json against the logs and the files")
    ap.add_argument("--depth", type=int, default=1,
                    help="beam depth to run (1: what fits, see the docstring)")
    ap.add_argument("--timeout", type=int, default=TIMEOUT)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    if a.command == "run":
        raise SystemExit(run(a.depth, a.timeout, a.force))
    if a.check:
        raise SystemExit(check(a.check))
    ap.error("give 'run', or --check REPORT.json")


if __name__ == "__main__":
    main()
