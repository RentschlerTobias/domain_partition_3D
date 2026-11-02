"""P1 keep-BL tistos: export the finished P1 result for viewing (ParaView).

Reads the P1 stages of the keep-BL tistos run

    /work/datagen_keepbl/runs/tistos_s1_0131/blocks.vtk   12-block structure
    /work/datagen_keepbl/runs/tistos_s1_0131/tfi.vtk      TFI refill
    /work/datagen_keepbl/runs/tistos_s1_0131/sample.npz   the sample + its
                                                          quality report
    /work/datagen_keepbl/inspect/keepbl_features.vtk       feature curves
                                                          (T33 export)

and writes three legacy ASCII VTK files to /work/datagen_keepbl/inspect/p1:

    blocks.vtk    the 12-block structure: CELL_DATA block_id (int) and
                  scaled_jacobian (float)
    tfi.vtk       the TFI refill (77372 cells) with scaled_jacobian
    overview.vtk  the blocks plus the feature curves as VTK_LINE cells
                  (CELL_DATA feature_curve_id), so one file shows the
                  structure and the curves together

No new dependencies: the VTK writing is the repo's own legacy-ASCII writer
style (as in `ovm_io.write_hex_vtk` / `export_vtk.write_vtk`) and the scaled
Jacobian comes from `ovm_io.scaled_jacobian` (vectorised copy in
`clean_blocks.scaled_jacobians`).  The `vtk` python module is NOT used.

    # write the three VTKs
    /opt/venv-sci/bin/python experimentell/hex3d_algohex/p1_export_vtk.py

    # verify them
    /opt/venv-sci/bin/python experimentell/hex3d_algohex/p1_export_vtk.py \
        --check /work/datagen_keepbl/inspect/p1
"""

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import clean_blocks as cb                                          # noqa: E402
import export_features_vtk as ef                                    # noqa: E402

VTK_VERTEX = 1
VTK_LINE = 3
VTK_HEXAHEDRON = 12

RUN = Path("/work/datagen_keepbl/runs/tistos_s1_0131")
INSPECT = Path("/work/datagen_keepbl/inspect")
OUT = INSPECT / "p1"

BLOCKS_VTK = RUN / "blocks.vtk"
TFI_VTK = RUN / "tfi.vtk"
SAMPLE_NPZ = RUN / "sample.npz"
TET_VTK = RUN / "tet.vtk"
FEATURES_VTK = INSPECT / "keepbl_features.vtk"

EXPECTED_BLOCKS = 12
EXPECTED_TFI_CELLS = 77372
EXPECTED_FEATURE_EDGES = 580


# --------------------------------------------------------------------------
# reading
# --------------------------------------------------------------------------
def read_unstructured(path):
    """Read a legacy ASCII VTK UNSTRUCTURED_GRID.

    Returns (points, cells, types, cell_data) where cell_data is a dict
    name -> (kind, values) for every SCALARS block under CELL_DATA.
    """
    with open(path) as f:
        lines = f.read().split("\n")
    n_lines = len(lines)

    def next_line(i):
        while i < n_lines and (not lines[i].strip()
                                or lines[i].startswith("#")):
            i += 1
        return i

    i = next_line(0)
    while i < n_lines and not lines[i].startswith("POINTS"):
        i = next_line(i + 1)
    if i >= n_lines or not lines[i].startswith("POINTS"):
        raise ValueError("%s: no POINTS section" % path)
    npoints = int(lines[i].split()[1])
    pts = np.empty((npoints, 3), float)
    for k in range(npoints):
        pts[k] = [float(x) for x in lines[i + 1 + k].split()]

    i = next_line(i + 1 + npoints)
    if not lines[i].startswith("CELLS"):
        raise ValueError("%s: expected CELLS at line %d" % (path, i + 1))
    ncells = int(lines[i].split()[1])
    cells = []
    j = i + 1
    for _ in range(ncells):
        parts = lines[j].split()
        cells.append([int(x) for x in parts[1:1 + int(parts[0])]])
        j += 1

    j = next_line(j)
    if not lines[j].startswith("CELL_TYPES"):
        raise ValueError("%s: expected CELL_TYPES at line %d" % (path, j + 1))
    types = [int(lines[j + 1 + k].split()[0]) for k in range(ncells)]

    # CELL_DATA may follow; each SCALARS block holds one value per cell
    data = {}
    j = next_line(j + 1 + ncells)
    while j < n_lines and lines[j].startswith("CELL_DATA"):
        j = next_line(j + 1)
        while j < n_lines and lines[j].startswith("SCALARS"):
            tok = lines[j].split()
            name, kind = tok[1], tok[2]
            j = next_line(j + 1)
            if j < n_lines and lines[j].startswith("LOOKUP_TABLE"):
                j = next_line(j + 1)      # skip LOOKUP_TABLE default
            vals = [lines[j + k].split() for k in range(ncells)]
            if kind == "float":
                arr = np.array([float(v[0]) for v in vals])
            else:
                arr = np.array([int(v[0]) for v in vals])
            data[name] = (kind, arr)
            j = next_line(j + ncells)

    return pts, cells, types, data


# --------------------------------------------------------------------------
# writing
# --------------------------------------------------------------------------
def write_vtk(path, points, cells, cell_types, arrays, title):
    """Legacy VTK UNSTRUCTURED_GRID with per-cell scalars.

    ``arrays`` is a sequence of (name, kind, values) with kind in
    {"float", "int"}; one SCALARS block is written per entry, which is what
    ParaView's colour-by / threshold filter reads.  Same writer style as
    `canada_export_vtk.write_vtk`, minus the hex-only assumption so that
    the overview file can hold hexes and lines in one dataset.
    """
    with open(path, "w") as f:
        f.write("# vtk DataFile Version 2.0\n")
        f.write(f"{title}\n")
        f.write("ASCII\nDATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {len(points)} float\n")
        for p in points:
            f.write("%r %r %r\n" % (float(p[0]), float(p[1]), float(p[2])))
        total = sum(1 + len(c) for c in cells)
        f.write(f"CELLS {len(cells)} {total}\n")
        for c in cells:
            f.write(f"{len(c)} " + " ".join(str(int(v)) for v in c) + "\n")
        f.write(f"CELL_TYPES {len(cells)}\n")
        for t in cell_types:
            f.write(f"{int(t)}\n")
        f.write(f"CELL_DATA {len(cells)}\n")
        for name, kind, values in arrays:
            if len(values) != len(cells):
                raise ValueError("%s: array %s has %d values for %d cells"
                                 % (path, name, len(values), len(cells)))
            f.write(f"SCALARS {name} {kind} 1\nLOOKUP_TABLE default\n")
            if kind == "float":
                f.writelines(f"{float(v):.9g}\n" for v in values)
            else:
                f.writelines(f"{int(v)}\n" for v in values)
    size = os.path.getsize(path)
    print("[vtk] wrote %s  (%d pts, %d cells, %d bytes)"
          % (path, len(points), len(cells), size))
    return size


def hexes_of(cells, types):
    """The hexahedral cells of a parsed VTK (index list into `cells`)."""
    return [c for c, t in zip(cells, types) if t == VTK_HEXAHEDRON]


# --------------------------------------------------------------------------
# export
# --------------------------------------------------------------------------
def export(outdir=OUT, verbose=True):
    """Write blocks.vtk, tfi.vtk and overview.vtk.  Returns a summary dict."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    summary = {}

    # ---- the 12-block structure -----------------------------------------
    bp, bcells, btypes, bdata = read_unstructured(BLOCKS_VTK)
    bhex = hexes_of(bcells, btypes)
    if len(bhex) != len(bcells):
        raise ValueError("blocks.vtk holds non-hex cells")
    H = np.asarray(bhex, np.int64)
    if "block_id" not in bdata:
        raise ValueError("blocks.vtk has no CELL_DATA block_id")
    block_id = bdata["block_id"][1]
    sj = cb.scaled_jacobians(bp, H)
    write_vtk(outdir / "blocks.vtk", bp, bhex, [VTK_HEXAHEDRON] * len(bhex),
              [("block_id", "int", block_id),
               ("scaled_jacobian", "float", sj)],
              "P1 keep-BL tistos: 12-block structure "
              "(CELL_DATA block_id, scaled_jacobian)")
    summary["blocks"] = {
        "npoints": len(bp), "ncells": len(bhex),
        "block_ids": sorted(set(int(b) for b in block_id)),
        "sj_min": float(sj.min()), "sj_mean": float(sj.mean()),
        "n_inverted": int((sj <= 0).sum()),
    }

    # ---- the TFI refill --------------------------------------------------
    tp, tcells, ttypes, tdata = read_unstructured(TFI_VTK)
    thex = hexes_of(tcells, ttypes)
    if len(thex) != len(tcells):
        raise ValueError("tfi.vtk holds non-hex cells")
    TH = np.asarray(thex, np.int64)
    tsj = cb.scaled_jacobians(tp, TH)
    tblock = (tdata["block_id"][1] if "block_id" in tdata
              else np.zeros(len(thex), int))
    write_vtk(outdir / "tfi.vtk", tp, thex, [VTK_HEXAHEDRON] * len(thex),
              [("block_id", "int", tblock),
               ("scaled_jacobian", "float", tsj)],
              "P1 keep-BL tistos: TFI refill (CELL_DATA block_id, "
              "scaled_jacobian)")
    summary["tfi"] = {
        "npoints": len(tp), "ncells": len(thex),
        "block_ids": sorted(set(int(b) for b in tblock)),
        "sj_min": float(tsj.min()), "sj_mean": float(tsj.mean()),
        "n_inverted": int((tsj <= 0).sum()),
    }

    # ---- overview: blocks + the feature curves ---------------------------
    fp = ef.parse_vtk(FEATURES_VTK)
    fpoints = np.asarray(ef.read_points(str(FEATURES_VTK)), float)
    flines = [(c, cid) for c, t, cid in
              zip(fp["cells"], fp["types"], fp["colors"])
              if t == VTK_LINE]
    fverts = [c for c, t in zip(fp["cells"], fp["types"]) if t == VTK_VERTEX]

    # one point set: the blocks' points first, then the feature points that
    # are not already there (the feature VTK carries the tet mesh's vertices)
    offset = len(bp)
    extra = []
    used = set()
    for c, _cid in flines:
        for v in c:
            used.add(int(v))
    for c in fverts:
        for v in c:
            used.add(int(v))
    # keep the feature points that the blocks do not already contain, by
    # position: a block point is "the same" as a feature point if the
    # coordinates coincide within a tight tolerance
    bkey = {}
    for k, p in enumerate(bp):
        bkey.setdefault(_pkey(p), k)
    remap = {}
    new_pts = []
    for v in sorted(used):
        p = fpoints[v]
        k = _find_point(bp, bkey, p)
        if k is None:
            k = offset + len(new_pts)
            new_pts.append(p)
        remap[v] = k
    all_pts = np.vstack([bp, np.asarray(new_pts, float).reshape(-1, 3)]) \
        if new_pts else bp

    cells = list(bhex)
    types = [VTK_HEXAHEDRON] * len(bhex)
    curve_id = np.full(len(bhex), -1, int)
    for c, cid in flines:
        cells.append([remap[int(v)] for v in c])
        types.append(VTK_LINE)
        curve_id = np.append(curve_id, int(cid))
    for c in fverts:
        cells.append([remap[int(v)] for v in c])
        types.append(VTK_VERTEX)
        curve_id = np.append(curve_id, -1)

    # the per-cell arrays must have one entry per cell of the dataset: the
    # hexes carry block_id / scaled_jacobian, the curve cells carry -1
    n_cells = len(cells)
    ov_block_id = np.full(n_cells, -1, int)
    ov_block_id[:len(bhex)] = block_id.astype(int)
    ov_sj = np.full(n_cells, 0.0, float)
    ov_sj[:len(bhex)] = sj

    write_vtk(outdir / "overview.vtk", all_pts, cells, types,
              [("block_id", "int", ov_block_id),
               ("scaled_jacobian", "float", ov_sj),
               ("feature_curve_id", "int", curve_id)],
              "P1 keep-BL tistos: 12-block structure + feature curves "
              "(VTK_LINE cells, CELL_DATA feature_curve_id)")
    summary["overview"] = {
        "npoints": len(all_pts), "ncells": len(cells),
        "n_hex": len(bhex), "n_line": len(flines), "n_vertex": len(fverts),
        "curve_ids": sorted(set(int(c) for _c, c in flines)),
    }
    summary["features"] = {
        "npoints": len(fpoints), "n_edges": len(flines),
        "curve_ids": sorted(set(int(c) for _c, c in flines)),
    }

    if verbose:
        print_summary(summary, outdir)
    return summary


def _pkey(p, ndigits=9):
    """Rounded coordinate key, for matching two point sets."""
    return tuple(round(float(x), ndigits) for x in p)


def read_sample():
    """The sample's quality / provenance records (sample.npz), as dicts.

    Returns (quality, provenance); an empty dict if the file is missing or
    unreadable -- the sample is a cross-check, not a hard requirement.
    """
    if not SAMPLE_NPZ.exists():
        return {}, {}
    try:
        d = np.load(str(SAMPLE_NPZ), allow_pickle=True)
    except Exception:
        return {}, {}
    out = {}
    for key, name in (("quality", "quality"), ("provenance", "provenance")):
        if key not in d.files:
            continue
        try:
            out[name] = json.loads(str(d[key]))
        except Exception:
            out[name] = {}
    return out.get("quality", {}), out.get("provenance", {})


def _find_point(bp, bkey, p):
    """Index of `p` in `bp` (via the rounded-key table), or None."""
    k = bkey.get(_pkey(p))
    if k is None:
        return None
    if np.allclose(bp[k], p, rtol=0.0, atol=1e-9):
        return k
    # tolerate the tiny rounding of the ASCII round trip
    d = np.abs(bp - p).max(1)
    j = int(np.argmin(d))
    return j if d[j] < 1e-6 else None


def print_summary(summary, outdir=OUT):
    b, t, o = summary["blocks"], summary["tfi"], summary["overview"]
    print("\n[P1] exported to %s" % outdir)
    print("  blocks.vtk   %5d cells, %d block ids %s"
          % (b["ncells"], len(b["block_ids"]), b["block_ids"]))
    print("               scaled_jacobian min %.4f mean %.4f, %d inverted"
          % (b["sj_min"], b["sj_mean"], b["n_inverted"]))
    print("  tfi.vtk      %5d cells, %d block ids"
          % (t["ncells"], len(t["block_ids"])))
    print("               scaled_jacobian min %.4f mean %.4f, %d inverted"
          % (t["sj_min"], t["sj_mean"], t["n_inverted"]))
    print("  overview.vtk %5d cells (%d hex + %d line + %d vertex), "
          "%d points"
          % (o["ncells"], o["n_hex"], o["n_line"], o["n_vertex"], o["npoints"]))
    print("               feature curves %s, %d edges"
          % (summary["features"]["curve_ids"], summary["features"]["n_edges"]))


# --------------------------------------------------------------------------
# check
# --------------------------------------------------------------------------
def check(outdir=OUT):
    """Verify the three exported files.  Exit 0 only if all is well."""
    outdir = Path(outdir)
    ok = True

    def fail(msg):
        nonlocal ok
        ok = False
        print("MISMATCH: %s" % msg)

    names = ["blocks.vtk", "tfi.vtk", "overview.vtk"]
    for n in names:
        p = outdir / n
        if not p.exists():
            fail("%s does not exist" % p)
        elif p.stat().st_size == 0:
            fail("%s is empty" % p)
        else:
            print("%-13s %10d bytes" % (n, p.stat().st_size))
    if not ok:
        print("CHECK FAILED: %s" % outdir)
        return 1

    # ---- blocks ----------------------------------------------------------
    bp, bcells, btypes, bdata = read_unstructured(outdir / "blocks.vtk")
    bhex = hexes_of(bcells, btypes)
    if len(bhex) != len(bcells):
        fail("blocks.vtk holds non-hex cells")
    ids = sorted(set(int(v) for v in bdata["block_id"][1]))
    if len(ids) != EXPECTED_BLOCKS:
        fail("blocks.vtk has %d block ids %s, expected %d"
             % (len(ids), ids, EXPECTED_BLOCKS))
    else:
        print("blocks.vtk: %d cells, exactly %d block ids %s"
              % (len(bhex), len(ids), ids))
    for name in ("block_id", "scaled_jacobian"):
        if name not in bdata:
            fail("blocks.vtk has no CELL_DATA %s" % name)
    if "scaled_jacobian" in bdata:
        s = bdata["scaled_jacobian"][1]
        if len(s) != len(bhex):
            fail("blocks.vtk scaled_jacobian has %d values for %d cells"
                 % (len(s), len(bhex)))
        else:
            print("blocks.vtk: scaled_jacobian min %.4f mean %.4f"
                  % (s.min(), s.mean()))
    # the scaled Jacobian must be the one of the exported geometry
    sj = cb.scaled_jacobians(bp, np.asarray(bhex, np.int64))
    if "scaled_jacobian" in bdata:
        ref = bdata["scaled_jacobian"][1]
        if len(ref) == len(sj) and not np.allclose(ref, sj, atol=1e-6):
            fail("blocks.vtk scaled_jacobian does not match the mesh")

    # ---- cross-check against the sample's own quality record --------------
    quality, provenance = read_sample()
    q_blocks = quality.get("blocks")
    if q_blocks is not None and int(q_blocks) != len(ids):
        fail("sample.npz reports %d blocks, blocks.vtk has %d"
             % (int(q_blocks), len(ids)))
    elif q_blocks is not None:
        print("sample.npz: %d blocks (matches blocks.vtk)" % int(q_blocks))
    prov_blocks = provenance.get("blocks_vtk")
    if prov_blocks and Path(prov_blocks) != BLOCKS_VTK:
        print("NOTE: sample.npz provenance points at %s, not %s"
              % (prov_blocks, BLOCKS_VTK))

    # ---- tfi -------------------------------------------------------------
    tp, tcells, ttypes, tdata = read_unstructured(outdir / "tfi.vtk")
    thex = hexes_of(tcells, ttypes)
    if len(thex) != EXPECTED_TFI_CELLS:
        fail("tfi.vtk has %d cells, expected %d"
             % (len(thex), EXPECTED_TFI_CELLS))
    else:
        print("tfi.vtk: %d cells (expected %d)" % (len(thex), EXPECTED_TFI_CELLS))
    if "scaled_jacobian" not in tdata:
        fail("tfi.vtk has no CELL_DATA scaled_jacobian")
    else:
        s = tdata["scaled_jacobian"][1]
        if len(s) != len(thex):
            fail("tfi.vtk scaled_jacobian has %d values for %d cells"
                 % (len(s), len(thex)))
        else:
            print("tfi.vtk: scaled_jacobian min %.4f mean %.4f"
                  % (s.min(), s.mean()))

    # ---- overview: the feature lines must match keepbl_features.vtk ------
    op, ocells, otypes, odata = read_unstructured(outdir / "overview.vtk")
    olines = [c for c, t in zip(ocells, otypes) if t == VTK_LINE]
    if "feature_curve_id" not in odata:
        fail("overview.vtk has no CELL_DATA feature_curve_id")
        curve_of = {}
    else:
        curve_of = {}
        for c, t, cid in zip(ocells, otypes, odata["feature_curve_id"][1]):
            if t == VTK_LINE:
                curve_of[(int(c[0]), int(c[1]))] = int(cid)

    fp = ef.parse_vtk(FEATURES_VTK)
    fpoints = np.asarray(ef.read_points(str(FEATURES_VTK)), float)
    src = {}
    for c, t, cid in zip(fp["cells"], fp["types"], fp["colors"]):
        if t == VTK_LINE:
            src[(int(c[0]), int(c[1]))] = int(cid)
    if len(src) != EXPECTED_FEATURE_EDGES:
        print("NOTE: %s has %d feature edges, expected %d"
              % (FEATURES_VTK, len(src), EXPECTED_FEATURE_EDGES))
    if len(olines) != len(src):
        fail("overview.vtk has %d feature lines, %s has %d"
             % (len(olines), FEATURES_VTK.name, len(src)))
    else:
        print("overview.vtk: %d feature lines, %s has %d"
              % (len(olines), FEATURES_VTK.name, len(src)))

    # the curves must be the same edges: compare by geometry (the overview
    # renumbers the points, so the edge identity is checked on the geometry)
    op_key = {}
    for k, p in enumerate(op):
        op_key.setdefault(_pkey(p), k)

    def _key_pair(p0, p1):
        a, b = sorted((_pkey(p0), _pkey(p1)))
        return (a, b)

    src_pairs = {}
    for (a, b), cid in src.items():
        src_pairs[_key_pair(fpoints[a], fpoints[b])] = cid
    exp_pairs = {}
    for (a, b), cid in curve_of.items():
        exp_pairs[_key_pair(op[a], op[b])] = cid

    if set(exp_pairs) != set(src_pairs):
        only_e = sorted(set(exp_pairs) - set(src_pairs))[:5]
        only_s = sorted(set(src_pairs) - set(exp_pairs))[:5]
        fail("overview.vtk feature edges differ from %s "
             "(only in overview: %s, only in source: %s)"
             % (FEATURES_VTK.name, only_e, only_s))
    elif {k: v for k, v in exp_pairs.items()} != src_pairs:
        fail("overview.vtk feature_curve_id partition differs from %s"
             % FEATURES_VTK.name)
    else:
        print("overview.vtk: feature edges and feature_curve_id match %s "
              "exactly (%d edges, curves %s)"
              % (FEATURES_VTK.name, len(exp_pairs),
                 sorted(set(exp_pairs.values()))))

    if not ok:
        print("CHECK FAILED: %s" % outdir)
        return 1
    print("CHECK OK: blocks.vtk (%d blocks), tfi.vtk (%d cells), "
          "overview.vtk (%d cells, %d feature lines)"
          % (len(ids), len(thex), len(ocells), len(olines)))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default=str(OUT),
                    help="where the three VTKs are written (default %(default)s)")
    ap.add_argument("--check", metavar="DIR",
                    help="verify the exported VTKs in DIR instead of writing")
    args = ap.parse_args(argv)

    if args.check:
        return check(args.check)

    export(args.outdir)
    print("\nParaView:")
    print("  open %s/overview.vtk, %s/blocks.vtk, %s/tfi.vtk"
          % (args.outdir, args.outdir, args.outdir))
    print("  colour by 'block_id' (12 blocks) or 'feature_curve_id' "
          "(3 curves, -1 = hex)")
    print("  threshold 'scaled_jacobian' to see the worst cells")
    return 0


if __name__ == "__main__":
    sys.exit(main())
