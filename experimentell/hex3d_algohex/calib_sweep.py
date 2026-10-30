"""P4a calibration harness: drive AlgoHex (HexMeshing) over the calibration
geometries with fixed settings, collect wall time + hex count + valid_volume
per run, and check the smoke runs before any sweep is trusted.

    python calib_sweep.py run --geometry <name|path> --setting A|B|C \\
        --n 2000 --out /work/calib/sweep
    python calib_sweep.py --smoke-check /work/calib/smoke

Geometry names are resolved in this order: an existing path;
``/work/calib/hexme[/meshes]/<name>.vtk`` (the hexme pack, extracted from
data/hexme/hexme.zip); else ``data/charite/test_geometries_LV_aneurysm/<name>``
(Charité LV_all.stl case, which is tet-prepped first via charite_prep into
``/work/calib/prep/<name>/tet.vtk``, reused when it exists).

Settings (the calibration matrix for T4/T5):
  A  default                -- no extra HexMeshing flags
  B  -d 70                  -- analytic-free dihedral feature detection at
                               70 deg (AlgoHex default, PROGRESS.md)
  C  -d 70 --force-feature-threshold -- recompute the feature edges from the
       dihedral angle even when the input vtk already carries
       `AlgoHex::FeatureEdges`/`edge_colors`. NOTE: the flag is a boolean
       (main.cc `add_flag`); the threshold VALUE is -d, and 70 deg is the
       justified AlgoHex default -- a lower value tags partial creases as
       features and AlgoHex then hard-aligns the octahedral field to them
       (PROGRESS.md, classifySurfaces root cause).

Each run gets its own directory <out>/<geometry>_<setting>_<n>/ with
run_algohex.py called as a subprocess, pinned to 2 cores with taskset
(--backend native refuses --cpus, so CPU is limited outside the process);
the cores are --cores (default '0,1', override per lane with $CALIB_CORES
so several sweep lanes can run side by side on disjoint core pairs),
killed after --timeout seconds (default 7200). metrics.json (this script) and
hexmeshing.log (run_algohex.py) land in the run directory; metrics.json holds
`valid_volume` and `Parametrization` values from HexMeshing's own -j file plus
the hex count parsed from the output .ovm with ovm_io.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from collections import OrderedDict
from datetime import datetime, timezone
from pathlib import Path

REPO_DIR = Path(__file__).resolve().parent.parent.parent
HERE = Path(__file__).resolve().parent

HEXME_ZIP = Path("/opt/stack/domain_partition_3D/data/hexme/hexme.zip")
HEXME_DIR = Path("/work/calib/hexme")
CHARITE_ROOT = Path("/opt/stack/domain_partition_3D/data/charite/"
                    "test_geometries_LV_aneurysm")
PREP_DIR = Path("/work/calib/prep")
CHARITE_CASES = ("Case09_pre", "Case09_post", "Case199_pre", "Case199_post")

SETTINGS = OrderedDict((
    ("A", {"label": "default",
           "args": []}),
    ("B", {"label": "-d 70",
           "args": ["-d", "70"]}),
    # --force-feature-threshold is a boolean flag on this HexMeshing build
    # (main.cc add_flag, "Refresh feature edges even if the feature edge
    # property exists"); the threshold it refreshes with is -d. 70 deg is the
    # default PROGRESS.md justifies, so no separate value flag is needed.
    ("C", {"label": "-d 70 --force-feature-threshold",
           "args": ["-d", "70", "--force-feature-threshold"]}),
))


def resolve_hexme(name):
    """The extracted pack keeps the zip structure meshes/<name>.vtk."""
    for cand in (HEXME_DIR / f"{name}.vtk", HEXME_DIR / "meshes" / f"{name}.vtk",
                 HEXME_DIR / name, HEXME_DIR / "meshes" / name):
        if cand.is_file():
            return cand
    return None


def resolve_geometry(name):
    """-> ("hexme", path) | ("charite", case_dir) or None."""
    p = Path(name)
    if p.is_file():
        return "hexme", p
    if p.is_dir() and (p / "LV_all.stl").is_file():
        return "charite", p
    h = resolve_hexme(name)
    if h:
        return "hexme", h
    c = CHARITE_ROOT / name
    if (c / "LV_all.stl").is_file():
        return "charite", c
    return None


def extract_hexme(dest=HEXME_DIR, zip_path=HEXME_ZIP):
    """Unpack data/hexme/hexme.zip with python's zipfile (no unzip needed)."""
    dest.mkdir(parents=True, exist_ok=True)
    import zipfile
    with zipfile.ZipFile(zip_path) as z:
        for info in z.infolist():
            if info.is_dir():
                continue
            target = dest / info.filename
            target.parent.mkdir(parents=True, exist_ok=True)
            if target.exists() and target.stat().st_size == info.file_size:
                continue
            with z.open(info) as src, open(target, "wb") as dst:
                shutil.copyfileobj(src, dst)
    n = len(list(dest.glob("meshes/*.vtk")))
    print(f"[calib] hexme pack extracted to {dest}: {n} vtk meshes")


# Bump when charite_prep or its tet-meshing parameters change meaningfully:
# prep_charite redoes the prep when the cached marker doesn't match.
PREP_VERSION = "v1"


def prep_charite(case_dir):
    """STL case -> cached AlgoHex input tet.vtk (charite_prep, 2-core)."""
    name = Path(case_dir).name
    out = PREP_DIR / name
    tet = out / "tet.vtk"
    marker = out / "prep_marker.json"
    cur = {"prep_version": PREP_VERSION}
    if tet.is_file() and marker.is_file():
        try:
            if json.load(open(marker)) == cur:
                print(f"[calib] charite prep cache hit: {tet}")
                return tet
        except (json.JSONDecodeError, OSError):
            pass
        tet.replace(out / f"tet_prev_{datetime.now(timezone.utc):%Y%m%d}.vtk")
        print(f"[calib] charite prep cache stale ({marker.read_text()} "
              f"!= {cur}), redoing")
    out.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(HERE))
    import charite_prep
    t0 = time.time()
    charite_prep.main(case_dir, out)
    json.dump(cur, open(marker, "w"))
    print(f"[calib] charite prep {name}: {time.time() - t0:.1f}s -> {tet}")
    return tet


def default_cores():
    try:
        return ",".join(str(c) for c in sorted(os.sched_getaffinity(0))[:2])
    except AttributeError:
        return "0,1"


def valid_volume(param):
    """HexMeshing's own validity numbers. The plain 'valid_volume' key is
    written in some builds; the final one is 'valid_volume_seamless'
    (fraction of the domain covered by a valid seamless map)."""
    for k in ("valid_volume", "valid_volume_seamless"):
        v = param.get(k)
        if isinstance(v, (int, float)) and not isinstance(v, bool) and v > 0:
            return float(v)
    return None


def run_one(geometry, setting, n, out_dir, cores=None, timeout=7200,
            hexmeshing_bin=None):
    """One calibrated HexMeshing run. Returns the path to its metrics.json."""
    resolved = resolve_geometry(geometry)
    if resolved is None:
        sys.exit(f"[calib] unknown geometry {geometry!r}: no such path and "
                 f"nothing under {HEXME_DIR} or {CHARITE_ROOT}")
    kind, src = resolved
    if kind == "charite":
        in_vtk = prep_charite(src)
    else:
        in_vtk = src

    geo_name = Path(geometry).name
    run_dir = Path(out_dir) / f"{geo_name}_{setting}_{n}"
    run_dir.mkdir(parents=True, exist_ok=True)
    cores = cores or default_cores()

    cmd = ["taskset", "-c", cores, sys.executable,
           str(HERE / "run_algohex.py"),
           "--backend", "native", "--prefix", "hex",
           "--in-vtk", str(in_vtk), "--out-dir", str(run_dir)]
    if hexmeshing_bin:
        cmd += ["--hexmeshing-bin", hexmeshing_bin]
    # "--" always: argparse REMAINDER as run_algohex parses it fails on bare
    # option-like extras after the option block (documented usage is
    # "run_algohex ... -- <HexMeshing args>"); run_algohex filters it out.
    cmd += ["--"] + SETTINGS[setting]["args"] + ["-n", str(n)]
    print(f"[calib] {run_dir.name}: {kind} {src}")
    print(f"[calib] {' '.join(cmd)}")

    import ovm_io
    metrics = {"geometry": geo_name, "kind": kind, "source": str(src),
               "setting": setting, "extra_args": SETTINGS[setting]["args"],
               "n": n, "cores": cores, "timeout_s": timeout,
               "in_vtk": str(in_vtk), "run_dir": str(run_dir),
               "exit": -1, "wall_s": None, "valid_volume": None,
               "hex_cells": None, "skipped_cells": None,
               "parametrization": None,
               "stopped": None,
               "timestamp": datetime.now(timezone.utc).isoformat()}
    # The HexMeshing log itself is written by run_algohex (hexmeshing.log for
    # the empty tag); the wrapper log captures the python wrapper's stderr.
    run_log = open(run_dir / "calib_wrapper.log", "ab")
    t0 = time.time()
    # Popen + process-group kill: run_algohex spawns HexMeshing as its own
    # child, so on timeout the whole group must die (subprocess.run would kill
    # only the wrapper and orphan the meshing run).
    exit_code = 124                              # like GNU timeout's exit
    try:
        proc = subprocess.Popen(cmd, stdout=run_log,
                                stderr=subprocess.STDOUT,
                                start_new_session=True)
        exit_code = proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        print(f"[calib] TIMEOUT after {timeout}s, killing process group")
        try:
            os.killpg(os.getpgid(proc.pid), 9)
            exit_code = proc.wait()
        except (ProcessLookupError, NameError):
            pass
        metrics["stopped"] = "timeout"
    except Exception as exc:                      # noqa: BLE001 - logged below
        metrics["stopped"] = f"error: {exc}"
    metrics["wall_s"] = time.time() - t0
    metrics["exit"] = exit_code

    jm = run_dir / "hex_hex_metrics.json"
    ovm = run_dir / "hex_hex.ovm"
    logf = run_dir / "hexmeshing.log"
    if jm.is_file():
        try:
            hj = json.load(open(jm))
            param = hj.get("Parametrization") or {}
            metrics["parametrization"] = param
            metrics["valid_volume"] = valid_volume(param)
        except json.JSONDecodeError as exc:
            print(f"[calib] metrics json unreadable: {exc}")
            metrics["stopped"] = metrics.get("stopped") or \
                "hexmeshing metrics json corrupt"
    if ovm.is_file():
        try:
            pts, edges, faces, polys = ovm_io.read_ovm(ovm)
            hexes, skipped = ovm_io.ovm_to_cells(pts, edges, faces, polys)
            metrics["hex_cells"] = len(hexes)
            metrics["skipped_cells"] = skipped
        except Exception as exc:                  # noqa: BLE001
            print(f"[calib] ovm parse failed ({exc})")
    if logf.is_file():
        metrics["log"] = str(logf)

    json.dump(metrics, open(run_dir / "metrics.json", "w"), indent=1)
    print(f"[calib] {run_dir.name}: exit={metrics['exit']} "
          f"wall={metrics['wall_s']}s "
          f"valid_volume={metrics['valid_volume']} "
          f"hexes={metrics['hex_cells']}")
    return run_dir / "metrics.json"


def smoke_check(dirpath):
    """Two finished smokes with valid_volume > 0.9 and hexes > 0 -> exit 0."""
    root = Path(dirpath)
    # layout <root>/<geometry>/<run_dir>/metrics.json (also accept any depth)
    runs = sorted(root.rglob("metrics.json"))
    runs = [Path(m) for m in runs]
    print(f"[calib] --smoke-check {root}: {len(runs)} run(s) found")
    ok = True
    for mpath in runs:
        m = json.load(open(mpath))
        good = (m.get("exit") == 0 and m.get("valid_volume") is not None
                and m["valid_volume"] > 0.9 and (m.get("hex_cells") or 0) > 0)
        ok &= good
        run_name = Path(mpath).parent.name
        print(f"  {'OK ' if good else 'BAD'} {run_name} "
              f"exit={m.get('exit')} wall={m.get('wall_s')}s "
              f"valid_volume={m.get('valid_volume')} "
              f"hexes={m.get('hex_cells')}")
    if len(runs) < 2:
        print("[calib] FAIL: need 2 runs (s01c_cube, Case09_pre), "
              f"found {len(runs)}")
        return 1
    if not ok:
        print("[calib] FAIL: a run is missing, timed out, has "
              "valid_volume <= 0.9 or 0 hexes")
        return 1
    print("[calib] smoke check PASSED: both runs valid")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", nargs="?", help="run one (geometry, setting, n)")
    ap.add_argument("--geometry", help="hexme mesh name, Charité case name "
                                       "or a direct .vtk path")
    ap.add_argument("--setting", default="A", choices=list(SETTINGS),
                    help="HexMeshing flag set (default A)")
    ap.add_argument("--n", type=int, default=10000,
                    help="target hex cell count for HexMeshing -n")
    ap.add_argument("--out", default="/work/calib/sweep",
                    help="sweep root; the run goes to "
                         "<out>/<geometry>_<setting>_<n>/")
    ap.add_argument("--cores", default=os.environ.get("CALIB_CORES", "0,1"),
                    help="taskset core list for the run (default: $CALIB_CORES, "
                         "else '0,1')")
    ap.add_argument("--timeout", type=int, default=7200,
                    help="kill HexMeshing after this many seconds")
    ap.add_argument("--hexmeshing-bin", default=None,
                    help="HexMeshing binary (default: run_algohex's pick)")
    ap.add_argument("--extract-hexme", action="store_true",
                    help="just unpack data/hexme/hexme.zip to "
                         "/work/calib/hexme and exit")
    ap.add_argument("--smoke-check", metavar="DIR",
                    help="parse the smoke metrics under DIR, exit 0 iff both "
                         "runs finished with valid_volume > 0.9 and hexes "
                         "> 0")
    a = ap.parse_args()
    if a.smoke_check:
        sys.exit(smoke_check(a.smoke_check))
    if a.extract_hexme:
        extract_hexme()
    elif a.run == "run":
        if not a.geometry:
            ap.error("run needs --geometry <name|path>")
        run_one(a.geometry, a.setting, a.n, a.out, a.cores, a.timeout,
                a.hexmeshing_bin)
    else:
        ap.print_help()
        if a.run and a.run != "run":
            ap.error(f"unknown subcommand {a.run!r} (expected 'run')")
