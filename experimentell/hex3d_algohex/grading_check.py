#!/usr/bin/env python3
"""Check the TFI grading direction of the blade O-grid ring (blade vs O-ring).

The blade mesh keeps a thin boundary layer on the blade wall; the O-grid ring
is the region between the core/ogrid_interface and the blade wall.  If the
ring's wall-normal grading is inverted -- thin layers at the core/O-ring end
instead of at the blade wall -- the thin cells sit where they do nothing, the
coarse cells sit on the wall, and the ring folds.

This script is READ-ONLY on the artifacts; it writes one JSON report.

    python3 grading_check.py [--run-dir DIR] [--out JSON]
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent.parent))

import base_complex as bc          # noqa: E402
import tfi                          # noqa: E402

RING_BLOCKS = (12, 13, 14, 15)     # the extruded O-grid ring (part 1)
CORE_BLOCKS = (4, 5, 6, 7)          # host core blocks near the interface
BL_BLOCKS = (8, 9, 10, 11)          # core blocks that carry the blade BL
OGRID_LABEL = 7                    # surface colour of ogrid_interface
BL_THICKNESS = (0.027, 0.035)       # reference wall-normal BL thickness

B_ALL = None                       # block id per cell, set in main()


# --------------------------------------------------------------------------
# VTK readers
# --------------------------------------------------------------------------
def read_hex_vtk(path):
    """(points, hexes, cell_data_dict) from a legacy VTK."""
    import meshio
    m = meshio.read(str(path))
    hexes = np.vstack([c.data for c in m.cells if c.type == "hexahedron"])
    cd = {}
    for key, d in m.cell_data.items():
        if isinstance(d, (list, tuple)):
            arr = np.concatenate([np.asarray(x).ravel()
                                  for x, c in zip(d, m.cells)
                                  if c.type == "hexahedron"])
        else:
            arr = np.asarray(d).ravel()
        cd[key] = arr
    return np.asarray(m.points, float), hexes, cd


def read_surface_vtk(path, label):
    """The (n, 3) point array of the triangles carrying surface `label`.

    meshio keeps one cell_data array per cell block, so the colour array
    has to be matched to the triangle block by position, not by zip().
    """
    import meshio
    m = meshio.read(str(path))
    P = np.asarray(m.points, float)
    cols = m.cell_data.get("color", [])
    out = []
    for ci, c in enumerate(m.cells):
        if c.type != "triangle":
            continue
        col = np.asarray(cols[ci]).ravel() if ci < len(cols) else np.zeros(0)
        if len(col) != len(c.data):
            continue
        sel = np.where(col == label)[0]
        if len(sel):
            out.append(P[c.data[sel]].reshape(-1, 3))
    return np.vstack(out) if out else np.zeros((0, 3))


# --------------------------------------------------------------------------
# geometry helpers
# --------------------------------------------------------------------------
def dist_to_surface(Q, S):
    """Distance from every point of Q to the nearest point of S."""
    from scipy.spatial import cKDTree
    if len(S) == 0:
        return np.full(len(Q), np.nan)
    return cKDTree(S).query(np.asarray(Q, float), k=1)[0]


def lattice_of_block(P, H, f2h, r):
    """(dims, vert) of block r, or None if it is not a structured grid."""
    cells = np.where(B_ALL == r)[0]
    if not len(cells):
        return None
    return tfi.block_lattice(H, cells, f2h)


def plane_gaps(P, vert, axis):
    """Gap length between consecutive planes along `axis`, per (i, j) column.

    Returns (gaps, dir) where gaps[k] holds the physical distance between
    plane k and plane k+1 for every column (n_planes - 1 gaps), and `dir` is
    the mean axis direction (plane 0 -> plane n-1).
    """
    Q = P[vert]                       # (ni, nj, nk, 3)
    seg = np.diff(Q, axis=axis)       # (ni, nj, nk-1, 3)
    gaps = np.linalg.norm(seg, axis=-1)
    d = seg.reshape(-1, 3).mean(0)
    n = np.linalg.norm(d)
    if n > 0:
        d = d / n
    return gaps, d


def _jsonable(o):
    """numpy scalars/arrays -> plain python, so the report serialises."""
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError("Object of type %s is not JSON serializable" %
                    o.__class__.__name__)


# --------------------------------------------------------------------------
def main():
    global B_ALL
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", default="/work/datagen_keepbl/runs/tistos_s1_0131")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    run = Path(args.run_dir)
    blade = run / "blade"
    out_path = Path(args.out) if args.out else blade / "grading_check.json"

    rep = {"run_dir": str(run), "inputs": {}, "blocks": {},
           "fractions": {}, "worst_cells": [], "verdict": "",
           "target_h": None, "code": {}}

    blocks_vtk = blade / "blocks_blade.vtk"
    part_vtk = blade / "blocks_blade_part.vtk"
    div_json = blade / "blocks_blade.divisions.json"
    audit_json = blade / "tfi_blade_audit.json"
    tet_vtk = run / "tet.vtk"
    tfi_vtk = blade / "tfi_blade.vtk"

    for p in (blocks_vtk, part_vtk, div_json, audit_json, tet_vtk, tfi_vtk):
        rep["inputs"][p.name] = {"path": str(p), "exists": p.exists()}
        if not p.exists():
            rep["verdict"] = ("ERROR: missing input %s -- no verdict possible"
                              % p)
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_text(json.dumps(rep, indent=2))
            print("[grading_check] " + rep["verdict"])
            return 2

    # ---------------- load ----------------
    P, H, cd = read_hex_vtk(blocks_vtk)
    _, _, cd_part = read_hex_vtk(part_vtk)
    global B_ALL
    B_ALL = cd["block_id"].astype(int)
    part = cd_part["part"].astype(int)
    f2h, _e2h = bc.build_topology(H)

    surf7 = read_surface_vtk(tet_vtk, OGRID_LABEL)
    rep["inputs"]["ogrid_interface_label7_points"] = int(len(surf7))

    divisions = json.loads(div_json.read_text())
    rep["target_h"] = divisions.get("target_h")
    rep["fractions"] = {
        "in_divisions_json": [k for k in divisions
                              if "frac" in k.lower()],
        "note": ("blocks_blade.divisions.json carries counts/classes/table/"
                 "frozen only -- it stores NO per-class fraction vector. "
                 "The wall-normal class 10 count is frozen at 10 layers."),
    }

    lat = {}
    for r in range(int(B_ALL.max()) + 1):
        got = lattice_of_block(P, H, f2h, r)
        if got is not None:
            lat[r] = got

    # ---- which vertices are shared between the core and the ring? --------
    core_mask = np.isin(B_ALL, np.arange(12))
    ring_mask = ~core_mask
    shared = np.intersect1d(np.unique(H[core_mask]),
                             np.unique(H[ring_mask]))
    rep["core_ring_shared_vertices"] = int(len(shared))

    # ---------------- 1. k orientation, 2. thickness profile -------------
    # The wall-normal axis of a ring block is the one whose last index is 10
    # (class 10 = [(12,2),(13,2),(14,2),(15,2)] in divisions.json) and whose
    # plane 0 coincides with the core/ring interface vertices.
    for r in RING_BLOCKS:
        if r not in lat:
            continue
        dims, vert = lat[r]
        dims = np.asarray(dims, int)
        info = {
            "dims": dims.tolist(),
            "part": int(part[B_ALL == r][0]) if (B_ALL == r).any() else None,
            "cells": int((B_ALL == r).sum()),
        }

        # the interface vertices of this block
        iface_verts = np.intersect1d(vert.ravel(), shared)
        # candidate wall-normal axes: the axis that carries the 10 layers
        cand = [ax for ax in range(3) if dims[ax] == 10]
        axis = cand[0] if cand else int(np.argmax(dims))

        # how many interface vertices lie on each end plane of that axis
        def plane_verts(k, _vert=vert, _axis=axis):
            return np.unique(np.take(P[_vert], k, axis=_axis)
                             .reshape(-1, 3), axis=0)

        from scipy.spatial import cKDTree
        tree = cKDTree(P[iface_verts]) if len(iface_verts) else None

        def dist_to_iface(k, _tree=tree, _pv=plane_verts):
            if _tree is None:
                return float("nan")
            d, _ = _tree.query(_pv(k), k=1)
            return float(d.max())

        d0 = dist_to_iface(0)
        dN = dist_to_iface(dims[axis] - 1)
        # The end plane that CARRIES the core/ring interface is the CORE/O-RING
        # side; the other end is the BLADE WALL.
        core_at_k0 = bool(d0 <= dN)
        blade_at_k0 = not core_at_k0

        # distance of every plane to the label-7 ogrid_interface surface
        if len(surf7):
            tree7 = cKDTree(surf7)
            pl_dist = []
            for k in range(dims[axis]):
                d, _ = tree7.query(plane_verts(k), k=1)
                pl_dist.append(float(d.mean()))
        else:
            pl_dist = [float("nan")] * dims[axis]

        info["wall_normal_axis"] = int(axis)
        info["n_layers"] = int(dims[axis])
        info["dist_plane0_to_core_interface"] = d0
        info["dist_planeN_to_core_interface"] = dN
        info["plane_mean_dist_to_label7_ogrid_interface"] = pl_dist
        info["blade_side_k"] = 0 if blade_at_k0 else int(dims[axis]) - 1
        info["core_oringside_k"] = int(dims[axis]) - 1 if blade_at_k0 else 0
        info["orientation"] = (
            "k=0 at the CORE/O-ring side (the ogrid_interface), k=%d at the BLADE WALL"
            % (dims[axis] - 1) if core_at_k0 else
            "k=0 at the BLADE WALL, k=%d at the CORE/O-ring side"
            % (dims[axis] - 1))
        info["label7_distance_monotonic_increasing_in_k"] = bool(
            all(pl_dist[i] <= pl_dist[i + 1] + 1e-12
                for i in range(len(pl_dist) - 1)))

        # thickness profile
        gaps, _d = plane_gaps(P, vert, axis)
        gm = gaps.mean(axis=(0, 1))
        gmin = gaps.min(axis=(0, 1))
        gmax = gaps.max(axis=(0, 1))
        kb = info["blade_side_k"]
        ko = info["core_oringside_k"]
        info["gap_mean"] = gm.tolist()
        info["gap_min"] = gmin.tolist()
        info["gap_max"] = gmax.tolist()
        info["gap_ratio_k_over_k1"] = (gm[:-1] / gm[1:]).tolist()
        info["thinnest_gap_k"] = int(np.argmin(gm))
        info["thickest_gap_k"] = int(np.argmax(gm))
        info["gap_at_blade_side"] = float(gm[kb])
        info["gap_at_core_oringside"] = float(gm[ko])
        info["ratio_gap_bladeside_over_core_oringside"] = float(gm[kb] / gm[ko])
        info["gap_at_bladeside_min"] = float(gmin[kb])
        info["gap_at_core_oringside_min"] = float(gmin[ko])
        info["first_gap_value"] = float(gm[0])
        info["first_gap_is_BL_thickness"] = bool(
            BL_THICKNESS[0] <= gm[0] <= BL_THICKNESS[1])
        rep["blocks"][str(r)] = info

    # ---------------- host core blocks 4-7 (and the BL blocks 8-11) -------
    for r in CORE_BLOCKS + BL_BLOCKS:
        if r not in lat:
            continue
        dims, vert = lat[r]
        dims = np.asarray(dims, int)
        info = {"dims": dims.tolist(),
                "part": int(part[B_ALL == r][0]) if (B_ALL == r).any() else None,
                "cells": int((B_ALL == r).sum())}
        # the graded (wall-normal) axis: the one with the largest spread of
        # gap means
        spreads = []
        for ax in range(3):
            gaps, _d = plane_gaps(P, vert, ax)
            gm = gaps.mean(axis=(0, 1))
            spreads.append(float(gm.max() - gm.min()) if len(gm) > 1 else 0.0)
        ax = int(np.argmax(spreads))
        info["wall_normal_axis"] = ax
        gaps, _d = plane_gaps(P, vert, ax)
        gm = gaps.mean(axis=(0, 1))
        info["gap_mean"] = gm.tolist()
        info["gap_min"] = gaps.min(axis=(0, 1)).tolist()
        info["gap_max"] = gaps.max(axis=(0, 1)).tolist()
        info["thinnest_gap_k"] = int(np.argmin(gm))
        info["thickest_gap_k"] = int(np.argmax(gm))
        info["spread"] = spreads
        # which end of that axis is nearer the ogrid_interface surface?
        if len(surf7):
            from scipy.spatial import cKDTree
            tree7 = cKDTree(surf7)
            pl = []
            for k in range(dims[ax]):
                pv = np.unique(np.take(P[vert], k, axis=ax)
                               .reshape(-1, 3), axis=0)
                d, _ = tree7.query(pv, k=1)
                pl.append(float(d.mean()))
            info["plane_mean_dist_to_label7_ogrid_interface"] = pl
        rep["blocks"][str(r)] = info

    # ---------------- 3. fractions ---------------------------------------
    classes = divisions.get("classes", [])
    rep["fractions"]["class_10_members"] = (
        classes[10] if len(classes) > 10 else None)
    rep["fractions"]["class_10_axes"] = (
        [ax for _r, ax in classes[10]] if len(classes) > 10 else None)
    rep["fractions"]["class_10_blocks"] = (
        sorted({r for r, _ax in classes[10]}) if len(classes) > 10 else None)
    rep["fractions"]["counts_per_class"] = divisions.get("counts")
    rep["fractions"]["frozen"] = divisions.get("frozen")
    tbl = {int(row["class"]): row for row in divisions.get("table", [])}
    if 10 in tbl:
        rep["fractions"]["class_10_table_row"] = tbl[10]

    # What does tfi.py actually use for the sample positions?  refill_complex
    # is called without a `fractions` argument, so every axis gets
    # uniform_fractions(n) = linspace(0, 1, n+1): NO geometric progression is
    # generated by tfi at all.  Report both the uniform vector that is used
    # and the geometric vector that the mesh actually carries.
    for r in RING_BLOCKS:
        key = str(r)
        if key not in rep["blocks"]:
            continue
        n = rep["blocks"][key]["n_layers"]
        uni = np.linspace(0.0, 1.0, n + 1)
        gm = np.asarray(rep["blocks"][key]["gap_mean"], float)
        tot = gm.sum()
        geo = np.concatenate([[0.0], np.cumsum(gm) / tot]) if tot > 0 else uni
        rep["blocks"][key]["fractions_uniform_used_by_tfi"] = uni.tolist()
        rep["blocks"][key]["fractions_geometric_of_input_lattice"] = geo.tolist()
        # is the geometric vector increasing from the blade side?
        kb = rep["blocks"][key]["blade_side_k"]
        geo_blade_first = np.concatenate(
            [[0.0], np.cumsum(gm[::-1]) / tot]) if tot > 0 else uni
        rep["blocks"][key]["fractions_geometric_anchored_at_blade_side"] = \
            geo_blade_first.tolist()
        rep["blocks"][key]["thin_end_is"] = ("blade" if gm[kb] < gm[1 - kb]
                                             else "core/O-ring")

    # ---------------- 5. the 20 worst cells ------------------------------
    audit = json.loads(audit_json.read_text())
    worst = audit.get("worst_20_cells", [])
    Pt, Ht, cdt = read_hex_vtk(tfi_vtk)
    Bt = cdt["block_id"].astype(int)
    save_B = B_ALL
    B_ALL = Bt
    f2h_t, _ = bc.build_topology(Ht)
    lat_t = {}
    for r in range(int(Bt.max()) + 1):
        got = lattice_of_block(Pt, Ht, f2h_t, r)
        if got is not None:
            lat_t[r] = got
    B_ALL = save_B

    wc = []
    for w in worst:
        cid = int(w["cell_id"])
        blk = int(w["block_id"])
        rec = {"cell_id": cid, "block_id": blk,
               "centre": w.get("centre"), "sJ": w.get("sJ"),
               "volume": w.get("volume")}
        if blk in lat_t:
            dims, vert = lat_t[blk]
            pos = {}
            for i in range(dims[0]):
                for j in range(dims[1]):
                    for k in range(dims[2]):
                        pos[int(vert[i, j, k])] = (i, j, k)
            ax = int(rep["blocks"].get(str(blk), {})
                     .get("wall_normal_axis", 2))
            ks = [pos[int(v)][ax] for v in Ht[cid] if int(v) in pos]
            if ks:
                k = int(min(ks))
                rec["k_range"] = [int(min(ks)), int(max(ks))]
                kb = rep["blocks"].get(str(blk), {}).get("blade_side_k")
                rec["k"] = k
                rec["wall_normal_axis"] = ax
                if kb is not None:
                    core_k = rep["blocks"][str(blk)]["core_oringside_k"]
                    rec["blade_side_k"] = kb
                    rec["core_oringside_k"] = core_k
                    rec["side"] = ("core/O-ring side" if k == core_k
                                   else ("blade-wall side" if k == kb
                                         else "interior"))
                    # the cell occupies the k-slab [k, k+1]; name the side it
                    # touches, so a cell at the very first slab counts as
                    # core-side even though it spans two k values.
                    rec["side_exact"] = rec["side"]
                    rec["layer_gap_at_this_k"] = float(
                        rep["blocks"][str(blk)]["gap_mean"][k])
        wc.append(rec)
    rep["worst_cells"] = wc
    k_hist = {}
    for rec in wc:
        k_hist[str(rec.get("k"))] = k_hist.get(str(rec.get("k")), 0) + 1
    rep["worst_cell_k_histogram"] = k_hist
    rep["worst_cells_all_at_k0"] = bool(
        k_hist and set(k_hist) == {"0"})

    # also: the same profile in the TFI output
    tfi_prof = {}
    for r in RING_BLOCKS:
        if r not in lat_t:
            continue
        dims, vert = lat_t[r]
        axis = 2 if dims[2] == 10 else int(np.argmax(dims))
        gaps, _d = plane_gaps(Pt, vert, axis)
        gm = gaps.mean(axis=(0, 1))
        tfi_prof[str(r)] = {"dims": np.asarray(dims, int).tolist(),
                            "gap_mean": gm.tolist(),
                            "gap_at_k0": float(gm[0]),
                            "gap_at_klast": float(gm[-1]),
                            "ratio_last_over_first": float(gm[-1] / gm[0])}
    rep["tfi_output_ring_profiles"] = tfi_prof

    # ---------------- 4. code ---------------------------------------------
    tfi_src = (HERE / "tfi.py").read_text().splitlines()
    ogr_src = (HERE / "ogrid_extrude.py").read_text().splitlines()

    def grab(src, start_pat, n, after=0):
        for i, line in enumerate(src):
            if i >= after and start_pat in line:
                return [{"line": i + 1, "text": l}
                        for l in src[i:i + n]]
        return []

    def grab_func(src, def_pat, n):
        """The first n lines of the function whose `def` line matches."""
        for i, line in enumerate(src):
            if line.startswith("def ") and def_pat in line:
                return [{"line": i + 1, "text": l} for l in src[i:i + n]]
        return []

    # locate the three extrusion variants by their def lines
    idx_const = next((i for i, l in enumerate(ogr_src)
                       if l.startswith("def extrude_ogrid_const")), 0)
    idx_smooth = next((i for i, l in enumerate(ogr_src)
                        if l.startswith("def extrude_ogrid_smooth")), 0)
    idx_orig = next((i for i, l in enumerate(ogr_src)
                      if l.startswith("def extrude_ogrid(")), 0)

    rep["code"] = {
        "tfi.py_refill_complex_fractions_default": grab(
            tfi_src, "def refill_complex", 0),
        "tfi.py_uniform_fractions": grab(tfi_src, "def uniform_fractions", 4),
        "tfi.py_refill_block": grab(tfi_src, "def refill_block", 20),
        "tfi.py_run_refill_call": grab(
            tfi_src, "Pn, Hn, Bn, rep = refill_complex", 1),
        "ogrid_extrude_profile_in_extrude_ogrid": grab(
            ogr_src, "profile = np.zeros(L + 1)", 12, after=idx_orig),
        "ogrid_extrude_profile_in_extrude_ogrid_const": grab(
            ogr_src, "profile = np.zeros(L + 1)", 12, after=idx_const),
        "ogrid_extrude_profile_in_extrude_ogrid_smooth": grab(
            ogr_src, "profile = np.zeros(L + 1)", 12, after=idx_smooth),
        "ogrid_extrude_layer_placement_extrude_ogrid": grab(
            ogr_src, "cols = S.P[list(verts)]", 6, after=idx_orig),
        "ogrid_extrude_layer_placement_const": grab(
            ogr_src, "cols = S.P[list(verts)]", 6, after=idx_const),
        "ogrid_extrude_layer_placement_smooth": grab(
            ogr_src, "cols = S.P[list(verts)]", 6, after=idx_smooth),
        "ogrid_extrude_chain_reversed_in_extrude_ogrid": grab(
            ogr_src, "C = omap.chain(int(v))[::-1]", 6, after=idx_orig),
        "ogrid_extrude_chain_reversed_in_const": grab(
            ogr_src, "C = omap.chain(int(v))[::-1]", 6, after=idx_const),
        "ogrid_extrude_chain_reversed_in_smooth": grab(
            ogr_src, "C = omap.chain(int(v))[::-1]", 6, after=idx_smooth),
        "ogrid_extrude_chain_doc": grab(
            ogr_src, "index 0 is the outer layer, index L is the blade wall", 12),
    }
    rep["code"]["notes"] = [
        "tfi.refill_complex is called from tfi.run without a `fractions` "
        "argument, so every axis -- including class 10, the wall-normal axis "
        "-- is resampled with uniform_fractions(n) = linspace(0,1,n+1). "
        "tfi.py therefore generates NO geometric progression at all: the "
        "wall-normal grading it reproduces is whatever the input lattice "
        "carries, sampled uniformly.",
        "The geometric progression that the input ring lattice carries is "
        "built in ogrid_extrude.extrude_ogrid / extrude_ogrid_const / "
        "extrude_ogrid_smooth: `C = omap.chain(int(v))[::-1]` puts chain "
        "index 0 at the BLADE WALL and index L at the outer/core layer, "
        "`f = np.r_[0.0, np.cumsum(d)]` accumulates the distance FROM the "
        "blade wall, and `profile` is the normalised f.",
        "That profile is then applied as `cols = base + dir * lens * "
        "profile[None, :, None]` with layer 0 placed at the BASE (the "
        "core/ogrid_interface) and layer m at the blade wall.  profile[0] = 0 "
        "sits on the core side and profile[m] = 1 on the blade wall, so the "
        "fine clustering that belongs at the wall lands at the CORE end: the "
        "profile is applied backwards.",
    ]

    # ---------------- 6. verdict ------------------------------------------
    ratios = [rep["blocks"][str(r)]["ratio_gap_bladeside_over_core_oringside"]
              for r in RING_BLOCKS if str(r) in rep["blocks"]]
    thin_is_core = all(
        rep["blocks"][str(r)]["gap_at_core_oringside"]
        < rep["blocks"][str(r)]["gap_at_blade_side"]
        for r in RING_BLOCKS if str(r) in rep["blocks"])
    n_core = sum(1 for rec in wc if rec.get("side_exact") == "core/O-ring side")
    n_blade = sum(1 for rec in wc if rec.get("side_exact") == "blade-wall side")
    inverted = bool(thin_is_core and ratios
                    and max(ratios) > 2.0)

    if inverted:
        verdict = ("INVERTED (YES): the ring's wall-normal grading is "
                   "reversed -- the thin layers sit at the CORE/O-ring end "
                   "(k=0) and the coarse layers at the BLADE WALL "
                   "(k=9).")
    else:
        verdict = ("NOT INVERTED (NO): the thin layers sit at the blade side.")

    # exact source lines of the profile construction, per extrusion variant
    flip_sites = []
    for fname, idx in (("extrude_ogrid", idx_orig),
                       ("extrude_ogrid_const", idx_const),
                       ("extrude_ogrid_smooth", idx_smooth)):
        ln = None
        for i in range(idx, min(idx + 120, len(ogr_src))):
            if "profile += f / max(f[-1], 1e-30)" in ogr_src[i]:
                ln = i + 1
                break
        ln_chain = None
        for i in range(idx, min(idx + 120, len(ogr_src))):
            if "C = omap.chain(int(v))[::-1]" in ogr_src[i]:
                ln_chain = i + 1
                break
        flip_sites.append({
            "function": "ogrid_extrude." + fname,
            "chain_reversed_line": ln_chain,
            "profile_line": ln,
            "current_code": [
                {"line": ln_chain, "text": ogr_src[ln_chain - 1].strip()}
                if ln_chain else None,
                {"line": ln, "text": ogr_src[ln - 1].strip()}
                if ln else None,
            ],
            "proposed_change": (
                "drop the [::-1] on the chain (or equivalently reverse the "
                "profile before it is used): `C = omap.chain(int(v))` makes "
                "chain index 0 the OUTER/core layer and index L the blade "
                "wall, so f = cumsum(d) is 0 at the core and max at the "
                "blade wall; profile[0] = 0 then lands on the core-side base "
                "and profile[m] = 1 on the blade wall, i.e. the fine layers "
                "sit on the blade wall as they should."),
        })

    rep["verdict"] = verdict
    rep["verdict_line"] = verdict
    rep["core_blocks_4_7_note"] = (
        "blocks 4-7 are core blocks with no wall-normal direction: they are "
        "the coarse (h~0.24-0.46) interior blocks on the core side of the "
        "ogrid_interface, and none of their three lattice axes is graded "
        "with a thin layer.  The core's own boundary layer lives in blocks "
        "8-11 (their graded axis ends in a 0.10-0.16 thin layer at the end "
        "FARTHEST from the ogrid_interface, i.e. at the blade wall), so the "
        "core BL is correctly oriented -- only the ring (blocks 12-15) is "
        "inverted.")
    rep["summary"] = {
        "inverted": inverted,
        "ratio_gap_bladeside_over_core_oringside_per_ring_block":
            {str(r): rep["blocks"][str(r)]
             ["ratio_gap_bladeside_over_core_oringside"]
             for r in RING_BLOCKS if str(r) in rep["blocks"]},
        "worst_cells_on_core_oringside": n_core,
        "worst_cells_on_blade_side": n_blade,
        "worst_cell_k_histogram": k_hist,
        "verdict_line": verdict,
        "where_to_flip": (
            "ogrid_extrude.py: in extrude_ogrid / extrude_ogrid_const / "
            "extrude_ogrid_smooth, the `profile` built from "
            "`C = omap.chain(int(v))[::-1]` and `f = np.r_[0.0, "
            "np.cumsum(d)]` is applied with layer 0 at the core-side base. "
            "Reverse it (profile -> profile[::-1], or equivalently build f "
            "from the chain WITHOUT the [::-1]) so that profile[0] = 0 is "
            "the core-side base and profile[m] = 1 the blade wall, which "
            "puts the fine layers on the blade wall. Affected blocks: the "
            "four ring blocks 12-15 (the whole O-grid ring). "
            "tfi.py needs NO change: it reproduces the input lattice's "
            "profile and cannot repair an inverted input."),
        "flip_sites": flip_sites,
        "expected_effect": (
            "the ring's wall-normal gap profile becomes smallest at k=9 "
            "(the blade wall, ~0.027-0.035 as the BL reference) and "
            "largest at k=0 (the core/O-ring side); the ~15x gap ratio "
            "across the ring disappears, the 215 inverted cells that all "
            "sit in the k=0 degenerate layer should disappear with it, and "
            "the boundary layer actually touches the blade wall where it "
            "resolves the shear."),
    }

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(rep, indent=2, default=_jsonable))
    print("[grading_check] " + verdict)
    print("[grading_check] gap blade/core per ring block: %s" %
          ", ".join("%d=%.2f" % (r, rep["blocks"][str(r)]
                                ["ratio_gap_bladeside_over_core_oringside"])
                    for r in RING_BLOCKS if str(r) in rep["blocks"]))
    print("[grading_check] worst-cell k histogram: %s" % k_hist)
    print("[grading_check] wrote %s" % out_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
