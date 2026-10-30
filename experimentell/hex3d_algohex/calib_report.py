"""Aggregate the hexme calibration sweep (calib_sweep.py runs) into one report
and check the report against the tables in CALIBRATION.md.

    python calib_report.py --build /work/calib/hexme_sweep \\
        --out /work/calib/report_hexme.json
    python calib_report.py --tables /work/calib/report_hexme.json
    python calib_report.py --check /work/calib/report_hexme.json --rows 18

--build walks every <root>/<run_dir>/ of the sweep and gives each row the
fields CALIBRATION.md must show: status, valid_volume, hex count, inverted
cells (min corner scaled Jacobian <= 0, ovm_io.scaled_jacobian), n, setting, wall
time, and -- for every row that is not "works" -- the AlgoHex failure mode
verbatim (exit code + last non-empty line of hexmeshing.log). A run directory
without a metrics.json yet is a row with status "running", its exact taskset
command (taken from the calib_sweep launcher log) and the wall time elapsed so
far.

--tables prints the per-geometry markdown tables for CALIBRATION.md.

--check fails (exit 1) when any row of the report is missing a status, when the
status is not one of STATUSES, when the row count is not exactly --rows, or
when CALIBRATION.md has no matching table row (the table is absent or does not
carry that row's status). Exit 0 otherwise.

Status classification (deterministic, in classify()):
  timeout   the harness killed the run after --timeout (metrics stopped=
            "timeout"); calib_sweep.py writes exit 124 or -9 for these.
  segfault  HexMeshing died on SIGSEGV (exit -11 or 139).
  works     exit 0 with a hex mesh (hex_cells > 0) and valid_volume >= 0.9,
            the same validity bar as calib_sweep --smoke-check.
  diverged  everything else: the parametrization failed or covered less than
            90 % of the domain, the solver aborted, the process died on any
            other signal.
  running   no metrics.json yet: the run is still executing (or was killed
            outside the harness); the row carries its exact command.
"""

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CALIBRATION_MD = HERE / "CALIBRATION.md"

STATUSES = ("works", "diverged", "segfault", "timeout", "running")
# Geometry order for the tables: the sweep matrix order of the task.
GEOMETRY_ORDER = ("s01c_cube", "s16u_torus", "i29c_bracket", "i31u_dlr_f6",
                  "n04u_transition_prism", "i26c_gc_turbo_vol1")
SETTING_FLAGS = {"A": "default", "B": "-d 70",
                 "C": "-d 70 --force-feature-threshold"}


def classify(m):
    """One of STATUSES for one calib_sweep metrics dict (see module doc)."""
    if m.get("stopped") == "timeout":
        return "timeout"
    if m.get("exit") in (-11, 139):
        return "segfault"
    good = (m.get("exit") == 0 and (m.get("hex_cells") or 0) > 0
            and m.get("valid_volume") is not None
            and m["valid_volume"] >= 0.9)
    return "works" if good else "diverged"


def failure_mode(m, run_dir):
    """AlgoHex failure mode verbatim: exit code + last log line, or None."""
    if classify(m) == "works":
        return None
    log = run_dir / "hexmeshing.log"
    last = ""
    if log.is_file():
        lines = [l for l in log.read_text(errors="replace").splitlines()
                 if l.strip()]
        if lines:
            last = lines[-1].strip()
    return f"exit {m.get('exit')}: {last}"


def param_note(run_dir):
    """One diagnostic line from HexMeshing's own hex_hex_metrics.json."""
    p = run_dir / "hex_hex_metrics.json"
    try:
        m = json.load(open(p))
        par = m.get("Parametrization") or {}
        pv = par.get("parametric_volume_seamless",
                     par.get("parametric_volume"))
        ni = par.get("n_invalid_param_tets")
        if pv is None and ni is None:
            return None
        return (f"parametric_volume {pv:.3g}"
                f", n_invalid_param_tets {ni}")
    except (OSError, ValueError):
        return None


def parse_run_name(name):
    """(geometry, setting, n) from a <geo>_<setting>_<n> run dir name."""
    for geo in GEOMETRY_ORDER:
        if name.startswith(geo + "_"):
            rest = name[len(geo) + 1:]
            setting, _, n = rest.rpartition("_")
            try:
                return geo, setting, int(n)
            except ValueError:
                return geo, setting, None
    return None, None, None


def run_command(root, run_dir):
    """The exact taskset command of a run, from its calib_sweep launcher log."""
    name = re.sub(r"_\d+$", "", run_dir.name) + ".log"
    cand = root / "logs" / name
    logs = [cand] if cand.is_file() else sorted((root / "logs").glob("*.log"))
    for path in logs:
        try:
            for line in path.read_text(errors="replace").splitlines():
                if line.startswith("[calib] taskset") and str(run_dir) in line:
                    return line[len("[calib] "):].strip()
        except OSError:
            continue
    return None


def elapsed_s(run_dir):
    """Seconds the run's wrapper has been alive (for a running row)."""
    try:
        p = subprocess.run(["pgrep", "-f", f"out-dir {run_dir}"],
                           capture_output=True, text=True)
        pids = [x for x in p.stdout.split() if x]
        if not pids:
            return None
        q = subprocess.run(["ps", "-o", "etimes=", "-p", ",".join(pids)],
                           capture_output=True, text=True)
        vals = [int(v) for v in q.stdout.split() if v.strip().isdigit()]
        return float(max(vals)) if vals else None
    except (OSError, ValueError):
        return None


def inverted_cells(run_dir):
    """Hexes with min corner scaled Jacobian <= 0 in hex_hex.ovm, or None when
    there is no (parseable) output mesh."""
    import numpy as np                     # noqa: F401 - ovm_io needs it
    sys.path.insert(0, str(HERE))
    import ovm_io
    ovm = run_dir / "hex_hex.ovm"
    if not ovm.is_file():
        return None
    try:
        pts, edges, faces, polys = ovm_io.read_ovm(ovm)
        hexes, _ = ovm_io.ovm_to_cells(pts, edges, faces, polys)
        if len(hexes) == 0:
            return None
        return int(sum(1 for c in hexes
                       if ovm_io.scaled_jacobian(pts, np.asarray(c)) <= 0))
    except Exception as exc:                # noqa: BLE001 - report what broke
        print(f"[calib_report] inversion count failed for {run_dir.name}: "
              f"{exc}", file=sys.stderr)
        return None


def build(root, out_path):
    """Collect every <root>/<run_dir>/ into one report json."""
    root = Path(root)
    rows = []
    for run_dir in sorted(p for p in root.iterdir()
                          if p.is_dir() and p.name != "logs"):
        mpath = run_dir / "metrics.json"
        if mpath.is_file():
            m = json.load(open(mpath))
            rows.append({
                "geometry": m.get("geometry")
                or parse_run_name(run_dir.name)[0],
                "setting": m.get("setting"),
                "setting_flags": SETTING_FLAGS.get(m.get("setting"), ""),
                "n": m.get("n"),
                "status": classify(m),
                "valid_volume": m.get("valid_volume"),
                "hex_cells": m.get("hex_cells"),
                "inverted_cells": inverted_cells(run_dir),
                "wall_s": m.get("wall_s"),
                "wall_elapsed": False,
                "exit": m.get("exit"),
                "cores": m.get("cores"),
                "failure_mode": failure_mode(m, run_dir),
                "note": param_note(run_dir),
                "command": run_command(root, run_dir),
                "run_dir": str(run_dir),
            })
        else:
            geo, setting, n = parse_run_name(run_dir.name)
            rows.append({
                "geometry": geo,
                "setting": setting,
                "setting_flags": SETTING_FLAGS.get(setting, ""),
                "n": n,
                "status": "running",
                "valid_volume": None,
                "hex_cells": None,
                "inverted_cells": None,
                "wall_s": elapsed_s(run_dir),
                "wall_elapsed": True,
                "exit": None,
                "cores": None,
                "failure_mode": None,
                "note": None,
                "command": run_command(root, run_dir),
                "run_dir": str(run_dir),
            })
    rows.sort(key=lambda r: (str(r["geometry"]), str(r["setting"])))
    report = {"sweep": "hexme", "n_runs": len(rows), "rows": rows}
    json.dump(report, open(out_path, "w"), indent=1)
    print(f"[calib_report] wrote {out_path}: {len(rows)} rows")
    return report


def _fmt(v, spec="{}"):
    return "—" if v is None else spec.format(v)


def tables(report):
    """The per-geometry markdown tables for CALIBRATION.md, one string."""
    out = []
    rows = report["rows"]
    geos = [g for g in GEOMETRY_ORDER
            if any(r["geometry"] == g for r in rows)]
    geos += sorted({r["geometry"] for r in rows} - set(geos))
    for geo in geos:
        out.append(f"### {geo}\n")
        out.append("| setting | flags | n | status | valid_volume | hexes | "
                   "inverted | wall_s | failure mode (exit + last log line) | "
                   "note |")
        out.append("|---|---|---|---|---|---|---|---|---|---|")
        for r in sorted((r for r in rows if r["geometry"] == geo),
                        key=lambda r: str(r["setting"])):
            if r["status"] == "running":
                fail = "running: " + (r["command"] or "command unknown")
                wall = ("&gt;%d (running)" % r["wall_s"]
                        if r.get("wall_s") is not None else "—")
            else:
                fail = r["failure_mode"] or "—"
                wall = _fmt(r["wall_s"], "{:.1f}")
            out.append("| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |"
                       .format(r["setting"], r["setting_flags"], r["n"],
                               r["status"],
                               _fmt(r["valid_volume"], "{:.3f}"),
                               _fmt(r["hex_cells"]),
                               _fmt(r["inverted_cells"]), wall,
                               fail.replace("|", "\\|"),
                               (r.get("note") or "—").replace("|", "\\|")))
        out.append("")
    return "\n".join(out)


def check(report_path, n_rows):
    """0 iff the report has exactly n_rows statused rows and CALIBRATION.md
    has a table row carrying each of them."""
    report = json.load(open(report_path))
    rows = report.get("rows") or []
    problems = []
    if len(rows) != n_rows:
        problems.append(f"report has {len(rows)} rows, expected {n_rows}")
    seen = set()
    for i, r in enumerate(rows):
        tag = f"row {i} ({r.get('geometry')}/{r.get('setting')}/" \
              f"{r.get('n')})"
        if not r.get("status"):
            problems.append(f"{tag}: missing status")
        elif r["status"] not in STATUSES:
            problems.append(f"{tag}: status {r['status']!r} not one of "
                            f"{STATUSES}")
        key = (str(r.get("geometry")), str(r.get("setting")),
               str(r.get("n")))
        if key in seen:
            problems.append(f"{tag}: duplicate row")
        seen.add(key)

    md = ""
    if not CALIBRATION_MD.is_file():
        problems.append(f"table absent: no {CALIBRATION_MD}")
    else:
        md = CALIBRATION_MD.read_text()
    if md:
        # (geometry, setting, n) -> status, from the markdown tables. The
        # geometry is the first word of the nearest preceding heading.
        in_table = {}
        geo = None
        for line in md.splitlines():
            if line.startswith("#"):
                geo = line.lstrip("#").strip().split()[0] if \
                    line.lstrip("#").strip() else None
                continue
            if not (geo and line.startswith("|")):
                continue
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            if len(cells) < 4:
                continue
            if set(cells[0]) <= {"-", ":"}:      # separator row
                continue
            if cells[0] == "setting":            # header row
                continue
            if cells[0] not in SETTING_FLAGS:    # not a sweep table row
                continue
            in_table[(geo, cells[0], cells[2])] = cells[3]
        for r in rows:
            key = (str(r.get("geometry")), str(r.get("setting")),
                   str(r.get("n")))
            if key not in in_table:
                problems.append(f"table row absent in CALIBRATION.md for "
                                f"{key[0]}/{key[1]}/n={key[2]}")
            elif in_table[key] != r.get("status"):
                problems.append(
                    f"table status {in_table[key]!r} != report status "
                    f"{r.get('status')!r} for {key[0]}/{key[1]}/n={key[2]}")

    if problems:
        print(f"[calib_report] CHECK FAILED ({len(problems)} problem(s)):")
        for p in problems:
            print(f"  - {p}")
        return 1
    print(f"[calib_report] check passed: {len(rows)} rows, all statused, "
          f"all present in {CALIBRATION_MD.name}")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build", metavar="RUN_ROOT",
                    help="collect <RUN_ROOT>/*/ into --out")
    ap.add_argument("--out", default="/work/calib/report_hexme.json",
                    help="report json to write with --build")
    ap.add_argument("--tables", metavar="REPORT",
                    help="print the per-geometry markdown tables of REPORT")
    ap.add_argument("--check", metavar="REPORT",
                    help="check REPORT's rows against CALIBRATION.md, "
                         "exit 1 on any missing status or missing table")
    ap.add_argument("--rows", type=int, default=18,
                    help="row count --check expects (default 18)")
    a = ap.parse_args()
    if a.build:
        rep = build(a.build, a.out)
        print(tables(rep))
    elif a.tables:
        print(tables(json.load(open(a.tables))))
    elif a.check:
        sys.exit(check(a.check, a.rows))
    else:
        ap.print_help()
