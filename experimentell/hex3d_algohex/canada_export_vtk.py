#!/usr/bin/env python3
"""T19 -- make the inverted cells of canada n = 100000 visible in ParaView.

P3b localized the 10 inverted cells of the canada n = 100000 hex mesh
and repaired 9 of them (`canada_repair_probe.py`, section in HANDOFF.md), but
a report is hard to argue about: nobody can see WHERE a cell with sJ -0.9906
sits, or whether the 1-ring around it is what blocks the repair. This module
turns the two saved mesh stages (`repair_probe_stages.npz`) into legacy VTK
files that carry the defect as cell data, so ParaView's threshold filter
shows it.

Three files in /work/canada/sweep/inspect/:

    hex_n100000.vtk      the pre-repair mesh -- the stage the beam measures,
                         i.e. after `fill_cavities` + `untangle`, which is
                         where the 10 inverted cells are. All hexes, cell
                         data `scaled_jacobian` (float) and `inverted`
                         (int, 1 for the 10).
    inverted_cells.vtk   just those 10 hexes plus their 1-ring neighbourhood,
                         the same two scalars plus `distance_to_inverted`
                         (int): 0 = inverted, 1 = shares a vertex with an
                         inverted cell. That is the whole patch a repair has
                         to move, a few hundred cells instead of 105k.
    repaired_n100000.vtk the post-repair mesh, same two scalars: 1 cell is
                         still inverted (the topologically locked half of
                         the doublet).

`scaled_jacobian` is the minimum corner scaled Jacobian
(`ovm_io.scaled_jacobian`, vectorised in `clean_blocks.scaled_jacobians`).
A hex can have a positive signed volume and still inverted corners, so this
-- not the volume -- is the validity test, and it is what the beam's
`min sJ >= 0.05` guard measures.

Nothing under /work/canada is read but the three inputs listed below, and
nothing is written outside the output directory.

usage:
  canada_export_vtk.py run [--out DIR] [--stages PATH] [--report PATH]
  canada_export_vtk.py --check DIR [--no-ovm]
"""

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import clean_blocks as cb  # noqa: E402
import ovm_io  # noqa: E402

# the canada lane, fixed by the task (same paths as canada_repair_probe.py)
SWEEP = Path("/work/canada/sweep")
OVM = SWEEP / "n100000" / "hex_n100000.ovm"
STAGES = SWEEP / "repair_probe_stages.npz"
REPORT = SWEEP / "repair_probe.json"
OUT = SWEEP / "inspect"

RAW_VTK = "hex_n100000.vtk"
SUBSET_VTK = "inverted_cells.vtk"
REPAIRED_VTK = "repaired_n100000.vtk"


# --------------------------------------------------------------------------
# the mesh stages
# --------------------------------------------------------------------------
def load_stages(path):
    """The two mesh stages the probe cached: (points, hexes) per stage.

    ``P0``/``H0`` is the mesh the beam measures (after `fill_cavities` and
    `untangle`, the stage with the 10 inverted cells), ``P1``/``H1`` is it
    after the repair.
    """
    with np.load(path) as d:
        return {"P0": d["P0"], "H0": d["H0"], "P1": d["P1"], "H1": d["H1"]}


# --------------------------------------------------------------------------
# writing
# --------------------------------------------------------------------------
def write_vtk(path, points, hexes, arrays, title):
    """Legacy VTK UNSTRUCTURED_GRID of hexahedra with per-cell scalars.

    ``arrays`` is a sequence of (name, kind, values) with kind in
    {"float", "int"}; one SCALARS block is written per entry, which is what
    ParaView's threshold filter reads. Same writer style as
    `ovm_io.write_hex_vtk`, extended by the cell data.
    """
    with open(path, "w") as f:
        f.write("# vtk DataFile Version 2.0\n")
        f.write(f"{title}\n")
        f.write("ASCII\nDATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {len(points)} float\n")
        for p in points:
            f.write(f"{p[0]} {p[1]} {p[2]}\n")
        f.write(f"CELLS {len(hexes)} {9 * len(hexes)}\n")
        for c in hexes:
            f.write("8 " + " ".join(str(int(v)) for v in c) + "\n")
        f.write(f"CELL_TYPES {len(hexes)}\n")
        for _ in hexes:
            f.write("12\n")          # VTK_HEXAHEDRON
        f.write(f"CELL_DATA {len(hexes)}\n")
        for name, kind, values in arrays:
            f.write(f"SCALARS {name} {kind} 1\nLOOKUP_TABLE default\n")
            if kind == "float":
                f.writelines(f"{float(v):.9g}\n" for v in values)
            else:
                f.writelines(f"{int(v)}\n" for v in values)
    print(f"[vtk] wrote {path}  ({len(points)} pts, {len(hexes)} cells, "
          f"{len(arrays)} scalars)")


def cell_rings(hexes, seeds, rings=1):
    """Cells within ``rings`` vertex-adjacency steps of ``seeds``.

    Returns (selected, distance): the cell ids of the patch and, for each of
    them, the number of steps to the nearest seed (0 for a seed itself).
    """
    v2c = defaultdict(list)
    for i, c in enumerate(hexes):
        for v in c:
            v2c[int(v)].append(i)
    dist = {int(s): 0 for s in seeds}
    frontier = list(dist)
    for step in range(1, rings + 1):
        nxt = []
        for c in frontier:
            for v in hexes[c]:
                for d in v2c[int(v)]:
                    if d not in dist:
                        dist[d] = step
                        nxt.append(d)
        frontier = nxt
    return sorted(dist), dist


def export(outdir, stages, verbose=True):
    """Write the three VTKs. Returns a summary dict."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    d = load_stages(stages)
    P0, H0, P1, H1 = d["P0"], d["H0"], d["P1"], d["H1"]

    sj0 = cb.scaled_jacobians(P0, H0)
    sj1 = cb.scaled_jacobians(P1, H1)
    inv0 = (sj0 <= 0).astype(int)
    inv1 = (sj1 <= 0).astype(int)
    bad0 = np.where(inv0)[0]
    sel, dist = cell_rings(H0, bad0, rings=1)

    title = {RAW_VTK: "canada n=100000 hex mesh, pre-repair "
                      "(cell data: scaled_jacobian, inverted)",
             SUBSET_VTK: "canada n=100000: the inverted cells and their "
                         "1-ring (distance_to_inverted)",
             REPAIRED_VTK: "canada n=100000 hex mesh, post-repair "
                           "(cell data: scaled_jacobian, inverted)"}
    write_vtk(outdir / RAW_VTK, P0, H0,
              [("scaled_jacobian", "float", sj0), ("inverted", "int", inv0)],
              title[RAW_VTK])
    write_vtk(outdir / REPAIRED_VTK, P1, H1,
              [("scaled_jacobian", "float", sj1), ("inverted", "int", inv1)],
              title[REPAIRED_VTK])
    write_vtk(outdir / SUBSET_VTK, P0, H0[sel],
              [("scaled_jacobian", "float", sj0[sel]),
               ("inverted", "int", inv0[sel]),
               ("distance_to_inverted", "int", [dist[c] for c in sel])],
              title[SUBSET_VTK])
    if verbose:
        print(f"[vtk] pre-repair : {len(H0)} cells, {int(inv0.sum())} "
              f"inverted, min sJ {sj0.min():.5f}")
        print(f"[vtk] post-repair: {len(H1)} cells, {int(inv1.sum())} "
              f"inverted, min sJ {sj1.min():.5f}")
        print(f"[vtk] 1-ring around the {len(bad0)} inverted cells: "
              f"{len(sel)} cells ({sum(1 for c in sel if dist[c] == 1)} at "
              f"distance 1)")
    return {"pre": (len(H0), int(inv0.sum()), float(sj0.min())),
            "post": (len(H1), int(inv1.sum()), float(sj1.min())),
            "subset": (len(sel), sorted(int(b) for b in bad0))}


# --------------------------------------------------------------------------
# reading back (plot_stages.read_vtk carries one int array only, so it cannot
# check a float one)
# --------------------------------------------------------------------------
def read_vtk(path):
    """Read a legacy UNSTRUCTURED_GRID written by `write_vtk`.

    Returns (points, cells, cell_types, arrays) with arrays as
    {name: np.ndarray}."""
    lines = Path(path).read_text().split("\n")

    def find(tag, start=0):
        for j in range(start, len(lines)):
            if lines[j].startswith(tag):
                return j
        raise ValueError(f"{path}: no {tag} section")

    j = find("POINTS")
    n = int(lines[j].split()[1])
    P = np.array([[float(x) for x in lines[j + 1 + k].split()]
                  for k in range(n)])
    j = find("CELLS")
    nc = int(lines[j].split()[1])
    cells, k = [], j + 1
    for _ in range(nc):
        pr = [int(x) for x in lines[k].split()]
        cells.append(pr[1:1 + pr[0]])
        k += 1
    j = find("CELL_TYPES")
    ty = np.array([int(lines[j + 1 + c]) for c in range(nc)])
    j = find("CELL_DATA")
    arrays, k = {}, j + 1
    while k < len(lines) and lines[k].startswith("SCALARS"):
        name, kind = lines[k].split()[1], lines[k].split()[2]
        if not lines[k + 1].startswith("LOOKUP_TABLE"):
            raise ValueError(f"{path}: {name} has no lookup table")
        arrays[name] = np.array([float(lines[k + 2 + c]) for c in range(nc)]) \
            if kind == "float" else \
            np.array([int(lines[k + 2 + c]) for c in range(nc)])
        k += 2 + nc
    return P, cells, ty, arrays


# --------------------------------------------------------------------------
# check: every claim of the deliverable against the files on disk
# --------------------------------------------------------------------------
class Bad(Exception):
    pass


def _need(cond, msg):
    if not cond:
        raise Bad(msg)


def check(directory, stages=STAGES, report=REPORT, ovm=OVM, verbose=True):
    """Verify the exported VTKs against the saved meshes and the report."""
    directory = Path(directory)
    _need(directory.is_dir(), f"{directory} does not exist")
    for name in (RAW_VTK, SUBSET_VTK, REPAIRED_VTK):
        _need((directory / name).is_file() and
              (directory / name).stat().st_size > 0,
              f"{directory / name} is missing or empty")

    rep = json.loads(Path(report).read_text())
    cells = rep["cells"]
    ids = sorted(int(c["id"]) for c in cells)
    _need(len(cells) == 10, f"{len(cells)} localized cells in the report, "
                            f"the task's lane has 10")

    d = load_stages(stages)
    P0, H0, P1, H1 = d["P0"], d["H0"], d["P1"], d["H1"]
    sj0 = cb.scaled_jacobians(P0, H0)
    sj1 = cb.scaled_jacobians(P1, H1)
    inv0 = sj0 <= 0
    inv1 = sj1 <= 0
    inv0_ids = [int(i) for i in np.where(inv0)[0]]
    inv1_ids = [int(i) for i in np.where(inv1)[0]]
    _need(sorted(inv0_ids) == ids,
          "the stage mesh's inverted cells are not the report's cells")
    _need(len(inv1_ids) == rep["repair"]["inverted_after"],
          "the repaired stage has a different inverted count than the report")
    sel, dist = cell_rings(H0, inv0_ids, rings=1)

    # --- the three files -------------------------------------------------
    got = {}
    for name in (RAW_VTK, SUBSET_VTK, REPAIRED_VTK):
        try:
            P, C, ty, A = read_vtk(directory / name)
        except (KeyError, IndexError, ValueError) as e:
            raise Bad(f"{name}: not a readable legacy VTK unstructured "
                      f"grid ({e})")
        _need(len(C) > 0, f"{name}: no cells")
        _need(int((ty == 12).sum()) == len(C),
              f"{name}: {len(C) - int((ty == 12).sum())} cells are not "
              f"VTK_HEXAHEDRON")
        for c in C:
            _need(len(c) == 8, f"{name}: a cell with {len(c)} vertices")
            _need(min(c) >= 0 and max(c) < len(P),
                  f"{name}: a cell indexes a point that is not there")
        for key in ("scaled_jacobian", "inverted"):
            _need(key in A, f"{name}: no cell data {key!r}")
        got[name] = (P, C, ty, A)

    # cell counts: the two full meshes are the two stages, the subset is the
    # 1-ring of the inverted cells of the pre-repair stage
    for name, (Pk, Hk) in ((RAW_VTK, (P0, H0)), (REPAIRED_VTK, (P1, H1))):
        P, C, _ty, A = got[name]
        _need(len(C) == len(Hk),
              f"{name}: {len(C)} cells, the mesh stage has {len(Hk)}")
        _need(len(C) == rep["input"]["cells"],
              f"{name}: {len(C)} cells, the report says "
              f"{rep['input']['cells']}")
        _need(len(P) == len(Pk),
              f"{name}: {len(P)} points, stage has {len(Pk)}")
        _need(np.allclose(P, Pk, atol=1e-12), f"{name}: points moved")
    P, C, _ty, A = got[SUBSET_VTK]
    _need(len(C) == len(sel),
          f"{SUBSET_VTK}: {len(C)} cells, the 1-ring of the inverted cells "
          f"has {len(sel)}")
    # the subset keeps the pre-repair point array, so a cell id in the subset
    # is the same cell id in the full mesh
    _need(len(P) == len(P0) and np.allclose(P, P0, atol=1e-12),
          f"{SUBSET_VTK}: points are not the pre-repair stage's")
    _need("distance_to_inverted" in A,
          f"{SUBSET_VTK}: no distance_to_inverted")

    # --- the inverted flag ----------------------------------------------
    for name, want in ((RAW_VTK, 10), (REPAIRED_VTK, 1), (SUBSET_VTK, 10)):
        n = int(got[name][3]["inverted"].sum())
        _need(n == want,
              f"{name}: {n} cells flagged inverted, expected {want}")
    flagged = [int(i) for i in np.where(got[RAW_VTK][3]["inverted"])[0]]
    _need(sorted(flagged) == ids,
          f"{RAW_VTK}: flags cells {sorted(flagged)}, "
          f"repair_probe.json has {ids}")
    left = [int(i) for i in np.where(got[REPAIRED_VTK][3]["inverted"])[0]]
    _need(sorted(left) == sorted(int(r["id"])
                                 for r in rep["repair"]["remaining"]),
          f"{REPAIRED_VTK}: flags {sorted(left)}, the report leaves "
          f"{[r['id'] for r in rep['repair']['remaining']]}")

    # --- the scaled Jacobians -------------------------------------------
    for name, sjk in ((RAW_VTK, sj0), (REPAIRED_VTK, sj1)):
        A = got[name][3]
        _need(np.allclose(A["scaled_jacobian"], sjk, atol=1e-6),
              f"{name}: scaled_jacobian does not match the mesh stage")
        _need(np.array_equal(A["inverted"], (sjk <= 0).astype(int)),
              f"{name}: inverted is not (scaled_jacobian <= 0)")
    A = got[SUBSET_VTK][3]
    _need(np.allclose(A["scaled_jacobian"], sj0[sel], atol=1e-6),
          f"{SUBSET_VTK}: scaled_jacobian does not match the mesh stage")
    _need(np.array_equal(A["inverted"], inv0[sel].astype(int)),
          f"{SUBSET_VTK}: inverted is not the stage's inverted cells")

    # --- the subset: right cells, right distances ------------------------
    D = A["distance_to_inverted"]
    _need(set(np.unique(D)) <= {0, 1},
          f"{SUBSET_VTK}: distance_to_inverted holds "
          f"{sorted(set(np.unique(D)) - {0, 1})}")
    _need(np.array_equal(D, [dist[c] for c in sel]),
          f"{SUBSET_VTK}: distance_to_inverted is not the 1-ring distance")
    _need(int((D == 0).sum()) == 10,
          f"{SUBSET_VTK}: {int((D == 0).sum())} cells at distance 0")
    _need(int((A["inverted"] == 1).sum()) == int((D == 0).sum()),
          f"{SUBSET_VTK}: the distance-0 cells are not the inverted ones")
    for c in cells:
        Q = P0[H0[c["id"]]]
        hit = [k for k in range(len(C))
               if np.allclose(P[C[k]].mean(0), Q.mean(0), atol=1e-9)]
        _need(len(hit) == 1, f"{SUBSET_VTK}: cell {c['id']} not found once")
        _need(abs(float(A["scaled_jacobian"][hit[0]]) - c["min_sj"]) < 1e-6,
              f"{SUBSET_VTK}: cell {c['id']}: sJ does not match the report")

    # --- the mesh really is the canada .ovm ------------------------------
    if ovm is not None:
        ovm = Path(ovm)
        if not ovm.exists():
            print(f"[check] {ovm} is not here, skipping the .ovm cross-check")
        else:
            P_ovm, e, f, poly = ovm_io.read_ovm(ovm)
            H_ovm, skipped = ovm_io.ovm_to_cells(P_ovm, e, f, poly)
            _need(skipped == 0, f"{ovm}: {skipped} polyhedra are not hexes")
            _need(len(H_ovm) == len(H0),
                  f"{ovm}: {len(H_ovm)} hexes, the stage has {len(H0)}")
            _need(len(P_ovm) == len(P0),
                  f"{ovm}: {len(P_ovm)} vertices, the stage has {len(P0)}")
            _need(Counter(frozenset(int(v) for v in c) for c in H_ovm) ==
                  Counter(frozenset(int(v) for v in c) for c in H0),
                  f"{ovm}: its hexes are not the exported stage's hexes")

    if verbose:
        print(f"OK {directory}: {RAW_VTK} {len(got[RAW_VTK][1])} cells, "
              f"{int(got[RAW_VTK][3]['inverted'].sum())} inverted; "
              f"{SUBSET_VTK} {len(got[SUBSET_VTK][1])} cells "
              f"({int((D == 1).sum())} at distance 1); "
              f"{REPAIRED_VTK} {len(got[REPAIRED_VTK][1])} cells, "
              f"{int(got[REPAIRED_VTK][3]['inverted'].sum())} inverted; "
              f"flagged ids {ids}")


# --------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mode", nargs="?", default="run", choices=("run",),
                    help="write the three VTKs (the default; --check "
                         "verifies them)")
    ap.add_argument("--check", metavar="DIR",
                    help="verify the exported VTKs in DIR")
    ap.add_argument("--out", default=str(OUT),
                    help="output directory of run")
    ap.add_argument("--stages", default=str(STAGES),
                    help="npz holding both mesh stages")
    ap.add_argument("--report", default=str(REPORT),
                    help="repair_probe.json the cell ids come from")
    ap.add_argument("--ovm", default=str(OVM),
                    help="the canada .ovm, cross-checked in --check")
    ap.add_argument("--no-ovm", action="store_true",
                    help="skip the .ovm cross-check (it costs ~30 s)")
    a = ap.parse_args(argv)
    if a.check:
        try:
            check(a.check, stages=Path(a.stages), report=Path(a.report),
                  ovm=None if a.no_ovm else Path(a.ovm))
        except Bad as e:
            print(f"FAIL {a.check}: {e}")
            return 1
        return 0
    export(a.out, a.stages)
    return 0


if __name__ == "__main__":
    sys.exit(main())
