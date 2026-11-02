#!/usr/bin/env python3
"""T24 -- does the sheet collapse need an O-grid / BL structure?

`beam_collapse.py` was developed and measured on the keep-BL samples (tistos
v9/v11, canada n=100000): every mesh there carries a boundary layer and an
O-grid, so a "sheet" was never tested on a mesh that has neither. This probe
asks the question directly: run the same beam collapse on the plain AlgoHex
outputs of the hexme sweep (T4), which are all-hex, carry no boundary layer and
no O-grid at all, and see whether it still produces a valid block structure.

hexme ships TETS. `/work/calib/hexme/meshes/*.vtk` are the tet INPUT meshes of
the sweep -- that is what `--input-vtk` wants (triangles + a per-triangle
colour, read by `clean_blocks.read_input_surface`). The all-hex meshes are the
AlgoHex OUTPUTS, `hex_hex.ovm`, under
`/work/calib/hexme_sweep/<geometry>_<setting>_<n>/`. What is there:

  * `s01c_cube_{A,B,C}_8000` -- 8000 hexes, 0 singular edges, 1 base block: a
    fully structured 20^3 cube. CONTROL (nothing to collapse).
  * `s16u_torus_{A,B}_8000` -- `hex_hex.ovm` is EMPTY (the parametrization
    diverged, final energy 12592 / seamless 3.6e8), so there is no hex mesh.
  * `s16u_torus_C_8000` -- 8742 hexes, 529 singular edges, 3 base blocks: the
    only torus mesh. It also holds 14 internal voids, 11 of which
    `fill_cavities` cannot close, so its boundary quad centres sit up to 3.70
    from the input surface and the beam's raw gate (0.15) refuses the run
    before any collapse. It is still run twice: once as specified (that is the
    answer for this mesh) and once with the two Hausdorff gates raised above
    the mesh's own 3.70 and everything else identical, so the collapse itself
    is measured too.
  * `n04u_transition_prism_{A,B}_8000` -- 8100 hexes, 200 singular edges, 13
    base blocks, no cavities, boundary Hausdorff 0.006. All-hex, no BL, no
    O-grid. This is the mesh that answers the question as specified.
  * `n04u_transition_prism_C_8000` -- 11434 hexes, 245 singular edges, 44 base
    blocks, Hausdorff 0.040, same family: a second data point.
  * `i26c_*`, `i29c_*`, `i31u_*` -- logs only, no hex output at all.

Each case runs `beam_collapse.py --guard struct --hausdorff-tol 0.15 --width 4
--depth 4 --workers 2` (timeout 1800, pinned to $CALIB_CORES, default 92,93),
then `tfi.py --require-lattices --solve-divisions` on the result -- whether
every block of the collapsed mesh carries a lattice is the real question, since
that is what stage 4 needs. Everything is recorded per case: the mesh it came
from, the base state the beam starts from (blocks, cells, singular edges,
inverted, min sJ, Hausdorff), the blocks after each depth, the final state, the
validator verdict, the TFI outcome and the wall time.

The keep-BL reference is conditional: the same parameters are run on the tistos
keep-BL mesh if T1 has left a `blocks.vtk` (searched under /work/calib). It had
not when this probe ran, which is recorded instead of guessed.

usage:
    python collapse_probe.py                 # run every case, write report.json
    python collapse_probe.py --check /work/calib/collapse_probe/report.json
    python collapse_probe.py --only n04u_transition_prism_A   # one case
"""

import argparse
import ast
import glob
import json
import os
import re
import subprocess
import sys
import time
from collections import defaultdict, deque
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import base_complex as bc                                            # noqa: E402
import clean_blocks as cb                                            # noqa: E402
import ovm_io                                                        # noqa: E402
import tet_prep_v5 as v5                                             # noqa: E402

PROBE_DIR = Path("/work/calib/collapse_probe")
SWEEP = Path("/work/calib/hexme_sweep")
HEXME = Path("/work/calib/hexme/meshes")

SCHEMA = 1
PY = sys.executable

# the task's parameters; only the torus' documented relaxation deviates
PARAMS = {"guard": "struct", "hausdorff_tol": 0.15, "width": 4, "depth": 4,
          "workers": 2, "timeout_s": 1800, "tfi_timeout_s": 900,
          "cores": os.environ.get("CALIB_CORES", "92,93")}

# variant -> the two Hausdorff gates. "as-specified" is the task's line;
# "relaxed" only lifts the gates above the mesh's own boundary deviation,
# because the torus' unclosable voids put its raw gate at 3.70.
GATES = {"as-specified": (None, None),
         "relaxed": (5.0, 5.0)}

CASES = [
    dict(name="s01c_cube_A", role="control", variant="as-specified",
         ovm=SWEEP / "s01c_cube_A_8000" / "hex_hex.ovm",
         input_vtk=HEXME / "s01c_cube.vtk",
         note="structured 20^3 cube, 0 singular edges, 1 base block"),
    dict(name="s16u_torus_C", role="probe", variant="as-specified",
         ovm=SWEEP / "s16u_torus_C_8000" / "hex_hex.ovm",
         input_vtk=HEXME / "s16u_torus.vtk",
         note="all-hex torus, 14 internal voids (11 unclosable)"),
]

# the same parameters on a keep-BL mesh, if T1 left one
TISTOS_GLOBS = ["/work/calib/tistos*/**/blocks*.vtk",
                "/work/calib/**/tistos*/**/blocks*.vtk"]

START_RE = re.compile(r"^\[beam\] start (\{.*\})$", re.M)
DEPTH_RE = re.compile(r"^\[beam\] depth (\d+): (\d+) admissible.*?"
                      r"best (\d+) blocks \((\d+)s\)", re.M)
RESULT_RE = re.compile(r"^\[beam\] result: (\{.*\}) via sheets", re.M)
VALID_RE = re.compile(r"^----- HexBlockValidator -----\n(VALID|INVALID)", re.M)
HD_RE = re.compile(r"^  boundary Hausdorff to input surface: "
                   r"([0-9.eE+-]+)$", re.M)
AFTER_RE = re.compile(r"^  (\w+) +(\S+)  ->  (\S+)$", re.M)
REJECT_RE = re.compile(r"^\[beam\] REJECT: (.*)$", re.M)
FILLED_RE = re.compile(r"^\[clean_blocks\] (\d+) cells filled", re.M)
TFI_NOT_LATTICE_RE = re.compile(r"^\[tfi\] (\d+) blocks are not lattices", re.M)
TFI_CELLS_RE = re.compile(r"^\[tfi\] (\S+): (\d+) cells, (\d+) blocks$", re.M)


# --------------------------------------------------------------------------
# the meshes
# --------------------------------------------------------------------------
def hex_mesh(ovm):
    """(points, hexes) of an all-hex mesh; raises if the file is empty.

    `.ovm` is an AlgoHex output, `.vtk` a written blocks.vtk -- the beam reads
    both (beam_collapse.py has the same branch), and both give one hex mesh."""
    if str(ovm).endswith(".vtk"):
        import tfi
        P, H, _B, _f2h = tfi.load_blocks(ovm)
    else:
        pts, e, f, poly = ovm_io.read_ovm(ovm)
        P, H = pts, ovm_io.ovm_to_cells(pts, e, f, poly)[0]
    if not len(H):
        raise RuntimeError(f"{ovm}: 0 hexes (the run diverged)")
    return P, H


def labeller_for(input_vtk):
    P, tri, tid = cb.read_input_surface(input_vtk)
    return cb.SurfaceLabeller(P, tri, tid, v5.NAMES)


def base_state(ovm, labeller):
    """The state the beam starts from, including its own preprocessing.

    `BlockStructure(fill=True, untangle_mesh=True)` fills cavities and untangles
    BEFORE the collapse, and the beam's raw gate measures that mesh, not the
    file. Both numbers are recorded: the raw file and what the beam sees."""
    P, H = hex_mesh(ovm)
    f2h, e2h = bc.build_topology(H)
    add = cb.fill_cavities(P, H, f2h, e2h)
    raw = {"hexes": int(len(H)), "cavities_filled": int(len(add))}
    if len(add):
        H = np.vstack([H, add])
    st = cb._structure_stats(P, H, labeller)
    return raw, {k: (round(v, 5) if isinstance(v, float) else int(v))
                 for k, v in st.items()}


def boundary_anatomy(P, H, labeller):
    """Where the boundary quads are: closed patches = cavities, open = surface.

    A closed quad patch that is not on the input surface is an internal void,
    and it is what puts a mesh's raw Hausdorff above the beam gate."""
    f2h, _ = bc.build_topology(H)
    loops = {}
    for h in range(len(H)):
        c = H[h]
        for fc in cb.HF:
            lp = [int(c[i]) for i in fc]
            fk = frozenset(lp)
            if len(f2h[fk]) == 1:
                loops[fk] = lp
    faces = list(loops)
    if not faces:
        return {"faces": 0, "components": 0, "closed": 0, "voids": 0,
                "dist_max": 0.0, "dist_median": 0.0}
    adj = defaultdict(list)
    for i, fk in enumerate(faces):
        for a in fk:
            adj[int(a)].append(i)
    seen, comps = set(), []
    for s in range(len(faces)):
        if s in seen:
            continue
        q, comp = deque([s]), []
        seen.add(s)
        while q:
            i = q.popleft()
            comp.append(i)
            for a in faces[i]:
                for j in adj[int(a)]:
                    if j not in seen:
                        seen.add(j)
                        q.append(j)
        comps.append(comp)
    ctrs = np.asarray([P[np.asarray(loops[faces[i]], int)].mean(0)
                       for i in range(len(faces))])
    T = labeller.tris
    d2, _ = cb._closest_point_on_tris(ctrs, labeller.P[T[:, 0]],
                                      labeller.P[T[:, 1]], labeller.P[T[:, 2]])
    d = np.sqrt(d2.min(1))
    n_void = 0
    for comp in comps:
        edges = defaultdict(int)
        for i in comp:
            lp = loops[faces[i]]
            for a, b in [(lp[k], lp[(k + 1) % 4]) for k in range(4)]:
                edges[(min(a, b), max(a, b))] += 1
        closed = all(n == 2 for n in edges.values())
        if closed and d[comp].min() > 0.5:
            n_void += 1
    return {"faces": len(faces), "components": len(comps), "closed": n_void,
            "voids": n_void, "dist_max": round(float(d.max()), 4),
            "dist_median": round(float(np.median(d)), 4)}


# --------------------------------------------------------------------------
# running the two stages
# --------------------------------------------------------------------------
def run_beam(case, outdir, params):
    """beam_collapse.py on one all-hex mesh; everything goes to <outdir>/beam.log."""
    outdir.mkdir(parents=True, exist_ok=True)
    raw_gate, ceiling = GATES[case["variant"]]
    cmd = ["taskset", "-c", params["cores"], PY,
           str(HERE / "beam_collapse.py"), str(case["ovm"]),
           "--input-vtk", str(case["input_vtk"]),
           "--out", str(outdir / "blocks.vtk"),
           "--guard", params["guard"],
           "--hausdorff-tol", str(ceiling if ceiling is not None
                                  else params["hausdorff_tol"]),
           "--width", str(params["width"]),
           "--depth", str(params["depth"]),
           "--workers", str(params["workers"])]
    if raw_gate is not None:
        cmd += ["--max-raw-hausdorff", str(raw_gate)]
    print(f"[probe] {case['name']}: {' '.join(cmd)}")
    log = outdir / "beam.log"
    t0 = time.time()
    with log.open("w") as fh:
        try:
            rc = subprocess.call(cmd, stdout=fh, stderr=subprocess.STDOUT,
                                 timeout=params["timeout_s"])
        except subprocess.TimeoutExpired:
            fh.write(f"\n[probe] TIMEOUT after {params['timeout_s']}s\n")
            rc = -1
    wall = time.time() - t0
    return parse_beam(log.read_text(), rc, wall, outdir / "blocks.vtk")


def _num(s):
    try:
        return int(s)
    except ValueError:
        return float(s)


def parse_beam(text, rc, wall, blocks):
    """Pull the beam's own numbers out of its log, keyed as the task asks."""
    out = {"exit": rc, "wall_s": round(wall, 1),
           "log": str(blocks.parent / "beam.log"), "blocks_vtk": str(blocks)}
    m = START_RE.search(text)
    out["initial"] = ast.literal_eval(m.group(1)) if m else None
    out["depth_blocks"] = [int(m.group(3)) for m in DEPTH_RE.finditer(text)]
    out["depth_seconds"] = [int(m.group(4)) for m in DEPTH_RE.finditer(text)]
    out["depth_admissible"] = [int(m.group(2)) for m in DEPTH_RE.finditer(text)]
    m = RESULT_RE.search(text)
    out["result"] = ast.literal_eval(m.group(1)) if m else None
    m = VALID_RE.search(text)
    out["valid"] = m.group(1) if m else None
    m = HD_RE.search(text)
    out["hausdorff"] = float(m.group(1)) if m else None
    after = {}
    for k, a, b in AFTER_RE.findall(text):
        after[k] = (_num(a), _num(b))
    out["after"] = after
    m = REJECT_RE.search(text)
    out["reject"] = m.group(1) if m else None
    m = FILLED_RE.search(text)
    out["cells_filled"] = int(m.group(1)) if m else 0
    out["blocks_written"] = bool(blocks.exists())
    # a run that produced nothing is a result, not a hole in the report
    if out["result"] is None:
        out["verdict"] = "no collapse: " + (out["reject"] or
                                            f"exit {rc}, see the log")
    else:
        out["verdict"] = "collapsed"
    return out


def run_tfi(blocks, input_vtk, outdir, params, target_h):
    """tfi.py --require-lattices: can stage 4 build a lattice per block?

    `--apply-divisions` (it implies `--solve-divisions`) is the full stage 4:
    a lattice per block, the conforming division counts, the refill. It exits 3
    when a block has no lattice, which is exactly the answer wanted."""
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = ["taskset", "-c", params["cores"], PY, str(HERE / "tfi.py"),
           str(blocks), "--input-vtk", str(input_vtk),
           "--apply-divisions", "--target-h", f"{target_h:.6g}",
           "--require-lattices"]
    print(f"[probe] tfi: {' '.join(cmd)}")
    log = outdir / "tfi.log"
    t0 = time.time()
    with log.open("w") as fh:
        try:
            rc = subprocess.call(cmd, stdout=fh, stderr=subprocess.STDOUT,
                                 timeout=params["tfi_timeout_s"])
        except subprocess.TimeoutExpired:
            fh.write(f"\n[probe] TIMEOUT after {params['tfi_timeout_s']}s\n")
            rc = -1
    text = log.read_text()
    out = {"exit": rc, "wall_s": round(time.time() - t0, 1),
           "log": str(log), "target_h": round(target_h, 6),
           "lattices_built": rc == 0, "missing": None}
    m = TFI_NOT_LATTICE_RE.search(text)
    if m:
        out["missing"] = int(m.group(1))
        out["lattices_built"] = False
    m = TFI_CELLS_RE.search(text)
    if m:
        out["cells"], out["blocks"] = int(m.group(2)), int(m.group(3))
    if rc not in (0, 3):
        # the CLI failed for a reason of its own: ask tfi directly, so that
        # "no lattice" and "the driver broke" stay distinguishable
        out["lattices_built"] = tfi_lattices_in_process(blocks)
        out["note"] = f"tfi.py exited {rc}; lattices() called in process"
    out["verdict"] = ("every block has a lattice" if out["lattices_built"]
                      else f"{out.get('missing')} blocks without a lattice")
    return out


def tfi_lattices_in_process(blocks):
    import tfi
    P, H, B, f2h = tfi.load_blocks(blocks)
    _lat, missing = tfi.lattices(P, H, B, f2h, verbose=False)
    return not len(missing)


def target_h(blocks):
    """Median edge length of the collapsed mesh: what TFI refills towards."""
    import tfi
    P, H, _B, _f2h = tfi.load_blocks(blocks)
    return float(np.median(np.linalg.norm(P[H[:, 1]] - P[H[:, 0]], axis=1)))


# --------------------------------------------------------------------------
# the keep-BL reference (only if T1 left one)
# --------------------------------------------------------------------------
def find_tistos_blocks():
    """T1's keep-BL deliverable, if it exists: (blocks.vtk, surface.vtk).

    The beam needs an `--input-vtk` to label the boundary, so a blocks.vtk only
    counts when its run directory also holds the surface it was built from."""
    for g in TISTOS_GLOBS:
        for blocks in sorted(Path(p) for p in
                             glob.glob(g, recursive=True)):
            surface = blocks.parent / "surface.vtk"
            if surface.exists():
                return blocks, surface
    return None, None


def tistos_case(blocks, surface):
    """The same beam parameters on a keep-BL mesh, for comparison."""
    return dict(name=f"keep_bl_{blocks.parent.name}", role="keep_bl",
                variant="as-specified", ovm=blocks, input_vtk=surface,
                note="tistos keep-BL mesh (T1), the mesh is already collapsed")


# --------------------------------------------------------------------------
# the report
# --------------------------------------------------------------------------
def run(only=None, out=None, timeout=None):
    out = Path(out or PROBE_DIR / "report.json")
    params = dict(PARAMS)
    if timeout is not None:                    # re-run past the 1800 s budget
        params["timeout_s"] = int(timeout)
    # a re-run of one case merges into the report that is already there, so
    # the other cases are not lost and their budget stays what it was
    old = {}
    if only and out.exists():
        try:
            prev = json.loads(out.read_text())
        except Exception:
            prev = {}
        if prev.get("schema") == SCHEMA:
            old = {c["name"]: c for c in prev.get("cases", [])}
    t0 = time.time()
    rep = {"schema": SCHEMA, "probe": "T24 sheet collapse on all-hex meshes "
           "without O-grid / BL", "params": PARAMS, "cases": []}
    cases = [c for c in CASES if not only or c["name"] in only]
    tistos, tistos_surface = find_tistos_blocks()
    if tistos is not None and not only:
        cases.append(tistos_case(tistos, tistos_surface))
    recs = dict(old)

    for case in cases:
        print(f"\n[probe] === {case['name']} ({case['role']}, "
              f"{case['variant']}) -- {case['note']}")
        rec = {"name": case["name"], "role": case["role"],
               "variant": case["variant"], "note": case["note"],
               "ovm": str(case["ovm"]), "input_vtk": str(case["input_vtk"]),
               "timeout_s": params["timeout_s"]}
        if params["timeout_s"] != PARAMS["timeout_s"]:
            rec["timeout_note"] = (
                f"re-run past the {PARAMS['timeout_s']}s budget: the "
                "post-collapse untangle did not fit in it")
        try:
            lab = labeller_for(case["input_vtk"])
        except Exception as e:
            rec["error"] = f"cannot read the input surface: {e}"
            rep["cases"].append(rec)
            continue
        try:
            raw, base = base_state(case["ovm"], lab)
        except Exception as e:
            rec["error"] = f"cannot read the hex mesh: {e}"
            rep["cases"].append(rec)
            continue
        P, H = hex_mesh(case["ovm"])
        rec["raw_mesh"] = raw
        rec["base"] = base
        rec["boundary"] = boundary_anatomy(P, H, lab)
        rec["beam"] = run_beam(case, PROBE_DIR / case["name"], params)
        # a run that is rejected by the raw surface gate never reaches the
        # beam's own "[beam] start" line, so recover the start state from the
        # base the beam measured (identical keys).
        if rec["beam"].get("initial") is None and isinstance(rec.get("base"), dict):
            rec["beam"]["initial"] = dict(rec["base"])
        b = rec["beam"]
        if b["blocks_written"]:
            try:
                h = target_h(PROBE_DIR / case["name"] / "blocks.vtk")
                rec["tfi"] = run_tfi(PROBE_DIR / case["name"] / "blocks.vtk",
                                     case["input_vtk"],
                                     PROBE_DIR / case["name"], params, h)
            except Exception as e:
                rec["tfi"] = {"lattices_built": False, "missing": None,
                              "error": str(e)}
        else:
            rec["tfi"] = {"lattices_built": False, "missing": None,
                          "skipped": "the beam wrote no blocks.vtk"}
        recs[case["name"]] = rec
        print(f"[probe] {case['name']}: exit {b['exit']}, "
              f"blocks {b['initial'] and b['initial']['blocks']} -> "
              f"{b['result'] and b['result']['blocks']}, "
              f"validator {b['valid']}, tfi {rec['tfi'].get('verdict')}, "
              f"{b['wall_s']}s")

    rep["tistos_keep_bl"] = (
        {"available": True, "blocks_vtk": str(tistos),
         "note": "the same parameters on T1's keep-BL mesh"} if tistos
        else {"available": False,
              "reason": "no keep-BL blocks.vtk under /work/calib when the "
                        "probe ran (T1's reference run was still going); "
                        "the comparison is left open rather than guessed",
              "searched": TISTOS_GLOBS})
    names = [c["name"] for c in CASES] + [n for n in recs
                                          if n not in {c["name"] for c in CASES}]
    rep["cases"] = [recs[n] for n in names if n in recs]
    probes = [c for c in rep["cases"] if c["role"] == "probe"]
    ran = [c for c in probes if c.get("beam", {}).get("blocks_written")]
    ok = [c for c in ran if c["beam"].get("valid") == "VALID"
          and c.get("tfi", {}).get("lattices_built")]
    rep["summary"] = {
        "probe_cases": len(probes), "probe_cases_with_blocks": len(ran),
        "probe_cases_valid_with_lattices": len(ok),
        "depends_on_ogrid_bl": not ok,
        "initial_blocks_means": "connected regions of the mesh cut along "
                                "every sheet (mesh_sheets + blocks_from_cut), "
                                "not one block per hex",
        "rejected": [c["name"] for c in rep["cases"]
                     if c.get("beam", {}).get("exit") not in (0, None)],
    }
    if not probes:
        rep["answer"] = ("No probe case selected (--only ran the control "
                         "alone); nothing to conclude here.")
    elif ok:
        rep["answer"] = (
            f"No: {len(ok)} of {len(probes)} O-grid-free/BL-free all-hex meshes "
            f"({', '.join(c['name'] for c in ok)}) went through the same beam "
            f"collapse (guard struct, tol 0.15, width 4, depth 4) and came out "
            f"VALID with a TFI lattice in every block, so the method needs "
            f"all-hex + sheets, not an O-grid or a boundary layer.")
    elif ran:
        rep["answer"] = (
            f"Inconclusive: {len(ran)} of {len(probes)} O-grid-free meshes "
            f"collapsed but none came out VALID with lattices; see the cases.")
    else:
        rep["answer"] = (
            f"Blocked: none of the {len(probes)} O-grid-free meshes reached the "
            f"collapse -- the raw gate or the mesh itself stopped them "
            f"({rep['summary']['rejected']}).")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rep, indent=1) + "\n")
    print(f"\n[probe] wrote {out} ({time.time() - t0:.0f}s)")
    print(f"[probe] {rep['answer']}")
    return rep


# --------------------------------------------------------------------------
# check: every claim against the saved geometry
# --------------------------------------------------------------------------
class Bad(Exception):
    pass


def _need(cond, msg):
    if not cond:
        raise Bad(msg)


def _recompute(blocks_vtk, input_vtk):
    """The beam's own stat keys, recomputed straight from the mesh."""
    import tfi
    P, H, B, _f2h = tfi.load_blocks(blocks_vtk)
    sj = cb.scaled_jacobians(P, H)
    f2h, _e2h = bc.build_topology(H)
    loops = []
    for h in range(len(H)):
        c = H[h]
        for fc in cb.HF:
            lp = [int(c[i]) for i in fc]
            if len(f2h[frozenset(lp)]) == 1:
                loops.append(lp)
                break
    lab = labeller_for(input_vtk)
    lab.label(P[np.asarray(loops)].mean(1))
    return {"cells": int(len(H)), "blocks": int(B.max()) + 1,
            "inverted": int((sj <= 0).sum()), "min_sj": float(sj.min()),
            "hausdorff": float(lab.max_dist)}


def check(path):
    """Verify the probe report against itself and against the saved meshes."""
    path = Path(path)
    _need(path.exists(), f"{path} does not exist")
    rep = json.loads(path.read_text())
    _need(rep.get("schema") == SCHEMA,
          f"schema is {rep.get('schema')!r}, expected {SCHEMA!r}")
    for key in ("params", "cases", "tistos_keep_bl", "summary", "answer"):
        _need(key in rep, f"missing section {key!r}")
    p = rep["params"]
    for key, val in (("guard", "struct"), ("hausdorff_tol", 0.15),
                     ("width", 4), ("depth", 4), ("workers", 2)):
        _need(p.get(key) == val, f"params.{key} is {p.get(key)!r}, the task "
                                 f"specifies {val!r}")
    cases = rep["cases"]
    _need(len(cases) >= 2, f"{len(cases)} cases, the probe needs at least a "
                           f"control and a probe")
    roles = {c["role"] for c in cases}
    _need("control" in roles, "no control case (the structured cube)")
    _need("probe" in roles, "no probe case (an O-grid-free all-hex mesh)")
    _need("keep_bl" in roles or not rep["tistos_keep_bl"]["available"],
          "the keep-BL run is claimed but its case is missing")
    names = [c["name"] for c in cases]
    _need(len(set(names)) == len(names), "duplicate case names")

    for c in cases:
        for key in ("name", "role", "variant", "ovm", "input_vtk", "note",
                    "base", "raw_mesh", "boundary", "beam", "tfi"):
            _need(key in c, f"{c.get('name')}: missing {key!r}")
        _need(Path(c["ovm"]).exists(), f"{c['name']}: {c['ovm']} does not exist")
        _need(Path(c["input_vtk"]).exists(),
              f"{c['name']}: {c['input_vtk']} does not exist")
        _need(c["variant"] in GATES, f"{c['name']}: unknown variant "
                                     f"{c['variant']!r}")
        _need(c.get("timeout_s", p["timeout_s"]) > 0,
              f"{c['name']}: no beam timeout recorded")
        base, b = c["base"], c["beam"]
        for key in ("cells", "blocks", "sing", "inverted", "min_sj",
                    "hausdorff"):
            _need(key in base, f"{c['name']}: base.{key} missing")
        _need(base["cells"] > 0, f"{c['name']}: the hex mesh has no cells")
        _need(base["blocks"] >= 1, f"{c['name']}: no base blocks")
        for key in ("exit", "wall_s", "initial", "depth_blocks", "result",
                    "valid", "hausdorff", "blocks_written"):
            _need(key in b, f"{c['name']}: beam.{key} missing")
        _need(isinstance(b["initial"], dict) and
              b["initial"]["cells"] == base["cells"],
              f"{c['name']}: the beam started from "
              f"{b['initial'] and b['initial'].get('cells')} cells, the mesh "
              f"has {base['cells']}")
        _need(len(b["depth_blocks"]) <= p["depth"],
              f"{c['name']}: {len(b['depth_blocks'])} depth lines for a "
              f"depth-{p['depth']} beam")
        _need(b["valid"] in ("VALID", "INVALID", None),
              f"{c['name']}: validator says {b['valid']!r}")
        if b["blocks_written"]:
            _need(b["exit"] == 0,
                  f"{c['name']}: exit {b['exit']} but a blocks.vtk exists")
            _need(b["valid"] in ("VALID", "INVALID"),
                  f"{c['name']}: blocks.vtk written but no validator verdict")
            _need(isinstance(b["result"], dict),
                  f"{c['name']}: no beam result although it wrote blocks.vtk")
            blocks = Path(b["blocks_vtk"])
            _need(blocks.exists(), f"{c['name']}: {blocks} is missing")
            got = _recompute(blocks, c["input_vtk"])
            for key in ("cells", "blocks", "inverted"):
                _need(got[key] == b["result"][key],
                      f"{c['name']}: blocks.vtk holds {got[key]} {key}, the "
                      f"report says {b['result'][key]}")
            _need(abs(got["min_sj"] - b["result"]["min_sj"]) < 1e-3,
                  f"{c['name']}: min sJ is {got['min_sj']:.4f} in the mesh, "
                  f"{b['result']['min_sj']} in the report")
            _need(abs(got["hausdorff"] - b["hausdorff"]) < 1e-3,
                  f"{c['name']}: Hausdorff is {got['hausdorff']:.4f} in the "
                  f"mesh, {b['hausdorff']} in the report")
            t = c["tfi"]
            _need(isinstance(t.get("lattices_built"), bool),
                  f"{c['name']}: tfi.lattices_built missing")
            _need(bool(t.get("verdict")), f"{c['name']}: tfi.verdict missing")
            _need(Path(t["log"]).exists(), f"{c['name']}: {t['log']} missing")
            if t["lattices_built"]:
                _need(t.get("exit") == 0, f"{c['name']}: tfi claims lattices "
                                          f"but exited {t.get('exit')}")
        else:
            _need(b["exit"] != 0, f"{c['name']}: no blocks.vtk but the beam "
                                  f"exited 0")
            _need(bool(b.get("reject") or b.get("verdict")),
                  f"{c['name']}: a run that wrote nothing has no reason")
            _need(c["tfi"].get("skipped"),
                  f"{c['name']}: TFI ran on a mesh that was never written")

    s = rep["summary"]
    probes = [c for c in cases if c["role"] == "probe"]
    _need(s["probe_cases"] == len(probes),
          f"summary.probe_cases is {s['probe_cases']}, there are {len(probes)}")
    ran = [c for c in probes if c["beam"]["blocks_written"]]
    ok = [c for c in ran if c["beam"]["valid"] == "VALID"
          and c["tfi"]["lattices_built"]]
    _need(s["probe_cases_with_blocks"] == len(ran),
          f"summary.probe_cases_with_blocks is {s['probe_cases_with_blocks']}, "
          f"{len(ran)} cases wrote blocks.vtk")
    _need(s["probe_cases_valid_with_lattices"] == len(ok),
          f"summary.probe_cases_valid_with_lattices is "
          f"{s['probe_cases_valid_with_lattices']}, {len(ok)} cases qualify")
    _need(s["depends_on_ogrid_bl"] == (not ok),
          "summary.depends_on_ogrid_bl contradicts the cases")
    _need(bool(rep["answer"].strip()),
          "the answer is empty -- the probe must state its conclusion")

    kb = rep["tistos_keep_bl"]
    _need(kb.get("available") in (True, False),
          "tistos_keep_bl.available is not a boolean")
    if kb["available"]:
        kc = [c for c in cases if c["role"] == "keep_bl"]
        _need(len(kc) == 1, "the keep-BL case is claimed but not present")
        _need(kc and kc[0]["beam"]["blocks_written"],
              "the keep-BL run wrote no blocks.vtk")
    else:
        _need(bool(kb.get("reason")), "no keep-BL run and no reason recorded")

    print(f"OK {path}: {len(cases)} cases "
          f"({', '.join(f'{c['name']} {c['beam']['valid']}' for c in cases)}), "
          f"{len(ok)}/{len(probes)} O-grid-free cases VALID with lattices; "
          f"depends_on_ogrid_bl={s['depends_on_ogrid_bl']}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mode", nargs="?", default="run", choices=("run",),
                    help="run the probe (the default; --check verifies a report)")
    ap.add_argument("--check", metavar="PATH",
                    help="verify a probe report and its saved meshes")
    ap.add_argument("--out", default=str(PROBE_DIR / "report.json"))
    ap.add_argument("--only", metavar="NAME[,NAME...]",
                    help="run only these cases (default: all)")
    ap.add_argument("--timeout", type=int, default=None, metavar="S",
                    help="beam timeout in seconds (default 1800, the probe's "
                         "budget); a longer one is only for re-running a case "
                         "whose post-processing did not fit in the budget")
    a = ap.parse_args(argv)
    if a.check:
        try:
            check(a.check)
        except Bad as e:
            print(f"FAIL {a.check}: {e}")
            return 1
        return 0
    only = [s for s in (a.only or "").split(",") if s] or None
    run(only=only, out=a.out, timeout=a.timeout)
    return 0


if __name__ == "__main__":
    sys.exit(main())
