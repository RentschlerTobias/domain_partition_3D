#!/usr/bin/env python3
"""Viz of the canada singular-curve MERGE study (T86): one VTK per curve coloured by merge group.

Read-only over the study artefacts (it never re-runs the merge physics):

  curve_analysis.json    634 singular edges -> 15 open curves, parallel curve
                         pairs, the 373 edge-level pair candidates
  curve_merge_test.json  the 3 whole-curve merge attempts (pairs 1-5, 0-2, 0-5:
                         all failed; 1-5 and 0-2 on boundary preservation
                         (Hausdorff 0.892 / 0.612 vs tol 1e-3, +322 / +39
                         boundary vertices), 0-5 on new inverted cells)
  struct_beam/blocks.vtk the 210-block baseline as light geometric context
                         (same domain bbox as the n100000 hex mesh)

The curve polylines are not stored in the JSONs (summary + indices only), so
this script does what the analysers did: singular edges from
`base_complex.singular_edges` and the curves from
`viz_singular.analyze_curves` (the exact calls behind curve_analysis.json), on
`--mesh`. It asserts its curve table against the JSON (curve count, n_edges
per curve) and prints a WARNING on any mismatch; no curve is invented, the
files below only re-label the recovered geometry.

Merge groups = colour. One colour per group: two curves share it iff they
would be merged into one. A curve that appears in several tested attempts or
candidate pairs can carry only one colour, so groups form a partition: the
tested pairs claim their curves in attempt order (1-5, 0-2; 0-5 only if both
its curves are still free), the remaining curves pair greedily along the
analysis parallel-pair list (sorted by tangent deviation, then distance), and
leftovers are singleton groups. Every tested pair is also drawn separately in
the VTM highlight block, so a skipped attempt (0-5, lost to the partition by
0-2 and 1-5) stays visible. Arrays:

  status 0 = candidate only (parallel pair from the analysis, merge not tested)
  status 1 = whole-curve merge tested (all three failed; why is in the legend)
  status 2 = merged (none in this study, kept for future runs)

Outputs in `--outdir` (default /work/canada/sweep/pair_probe3/viz/):

  merge_curves.vtp                    PolyData, one line cell per singular
                                      curve, arrays curve_id, group_id
                                      (colour), merge_partner_curve_id
                                      (-1 = none), status, pair_distance,
                                      n_edges
  merge_curves.block_curves.vtp       same data (the VTM curves block)
  merge_curves.block_context.vtp      blocks.vtk boundary faces (light context)
  merge_curves.block_tested_pairs.vtp only the curves of the 3 tested pairs,
                                      with array pair_tag; render thicker
  merge_curves.vtm                    multi-block: curves + beam context +
                                      the 3 tested pairs
  merge_curves_legend.json            colour -> group -> curve ids -> partner
                                      -> status -> what failed and why, one
                                      entry per group, plus the exact camera
                                      of merge_curves.png
  merge_curves.png                    render of the curves on the beam from
                                      the separable view

Run:  PY experimentell/hex3d_algohex/viz_merge_curves.py
      (paths defaulted; no vtk/paraview needed -- the VTK XML PolyData and
      multi-block files are written by hand, the PNG via matplotlib)
"""

import argparse
import json
import re
from collections import Counter
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent.parent
import sys
sys.path.insert(0, str(REPO / "experimentell" / "hex3d_algohex"))

import base_complex as bc
import curve_merge as cm
import ovm_io
import viz_singular as vs

OUT_VTP = "merge_curves.vtp"
OUT_BLOCKS = {
    "curves": "merge_curves.block_curves.vtp",
    "context": "merge_curves.block_context.vtp",
    "tested_pairs": "merge_curves.block_tested_pairs.vtp",
}
OUT_VTM = "merge_curves.vtm"
OUT_LEGEND = "merge_curves_legend.json"
OUT_PNG = "merge_curves.png"

STATUS_NAME = {0: "candidate-only", 1: "tested-failed", 2: "merged"}

# 15 distinct colours, one per group at most (a study of >= 16 curves still
# wraps around; the legend json carries the exact hex of every group).
PALETTE = [
    "#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4",
    "#33e0e0", "#f032e6", "#bcf60c", "#fabebe", "#008080",
    "#9a6324", "#800000", "#808469", "#000075", "#7f7f7f",
]

# The view where the pairs are separable. The curves run ~ (0.20, 0.10, 0.97)
# (mostly along z; the y span is only 3 of 15 units), so seen edge-on from a
# shallow elevation the near-coincident tested pair 1-5 (min dist 0.05) reads
# as two distinct lines over the beam.
CAMERA = {"elev": 12.0, "azim": -115.0, "roll": 0.0, "proj_type": "ortho"}
BOX_ASPECT = [1.0, 0.28, 1.28]    # x/y/z spans of the domain are 11, 3, 14
DOMAIN = {"x": (0.0, 11.0), "y": (-1.5, 4.5), "z": (-4.0, 10.0)}
SUB_FACES_KEEP = 6000             # context Poly3DCollection stays light


# ---------------------------------------------------------------- mesh
def load_mesh(mesh_path):
    pts, e, f, poly = ovm_io.read_ovm(str(mesh_path))
    P, H = pts.astype(float), ovm_io.ovm_to_cells(pts, e, f, poly)[0]
    return P, H.astype(np.int64)


def curve_polyline(P, curve):
    """Ordered vertex path of a traced curve (list of (u, v) edges) -> [n, 3]."""
    path = cm.curve_vertex_path(curve)
    if len(path) < 2:
        raise ValueError("curve trace has no path")
    return P[np.asarray(path, int)]


# ---------------------------------------------------------------- context
def _legacy_section(txt, header_regex, n):
    """`n` numbers after a legacy-vtk header line; generic helper."""
    m = re.search(header_regex, txt)
    if m is None:
        return None, None
    start = m.end()
    vals, count = [], 0
    for line in txt[start:].split("\n"):
        toks = line.split()
        if not toks or toks[0].isalpha():
            continue
        for t in toks:
            try:
                vals.append(float(t))
            except ValueError:
                break
        if len(vals) >= n:
            break
    return m, (np.array(vals[:n], float) if len(vals) >= n else None)


def beam_boundary_faces(blocks_vtk):
    """Blocks.vtk (legacy ASCII UNSTRUCTURED_GRID of hexes): the boundary
    faces (faced by exactly one hex). Returns (points, faces [F,4],
    block_id [F])."""
    txt = Path(blocks_vtk).read_text()
    m = re.search(r"POINTS (\d+) (?:float|double)\n", txt)
    npts = int(m.group(1))
    _, vals = _legacy_section(txt, r"POINTS (\d+) (?:float|double)\n", npts * 3)
    P = vals.reshape(-1, 3)

    mcell = re.search(r"CELLS (\d+) (\d+)\n", txt)
    ncell, total_size = int(mcell.group(1)), int(mcell.group(2))
    cell_size = total_size // ncell          # legacy: total len of cell list
    assert cell_size == 9 and total_size == ncell * 9, "expected hex cells"
    ktok = mcell.end()  # start of first cell line
    kend = txt.index("CELL_TYPES", ktok)
    cells = np.array([int(x) for x in txt[ktok:kend].split()])
    cells = cells[:ncell * cell_size].reshape(-1, cell_size)
    ids = cells[:, 1:]           # per cell: point count "8" + 8 point ids

    # per-cell block_id: the CELL_DATA SCALARS block_id int section
    di = txt.index("SCALARS block_id")
    after_hdr = txt.index("\n", txt.index("\n", di) + 1)   # skip + lookup line
    look = txt[after_hdr + 1:after_hdr + 6].split("\n")[0]
    data_start = after_hdr + 1 + (len(look) + 1 if look.startswith("LOOKUP") else 0)
    data_end = txt.find("FIELD", data_start)
    if data_end == -1:
        data_end = len(txt)
    block_id = np.array([int(x) for x in txt[data_start:data_end].split()])[:ncell]
    if len(block_id) != ncell:
        print(f"[merge_curves] WARNING: block_id scalars missing "
              f"({len(block_id)} of {ncell}); using cell index", flush=True)
        block_id = np.arange(ncell)

    # boundary faces: each (hex, HEX_FACES slot) quad faced by exactly one hex
    quads = []                                # ((point ids), block_id, first seen)
    first = {}
    cnt = Counter()
    for hi, row in enumerate(ids):
        for f4 in bc.HEX_FACES:
            fv = row[list(f4)]
            key = tuple(sorted(int(x) for x in fv))
            cnt[key] += 1
            if key not in first:
                first[key] = (fv.tolist(), int(block_id[hi]))
    for key, n in cnt.items():
        if n == 1:
            quads.append(first[key])
    faces = np.asarray([q[0] for q in quads], np.int64)
    bid = np.asarray([q[1] for q in quads], np.int64)
    return P, faces, bid


# ---------------------------------------------------------------- vtk xml
def _fmt(x):
    if isinstance(x, str):
        return x                      # String arrays: emit verbatim
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    return f"{float(x):.10g}"


def _xml_array(name, data, comp=1):
    arr = np.asarray(data)
    if arr.dtype.kind == "U":                    # string array (e.g. pair_tag)
        kind = "String"
    elif np.issubdtype(arr.dtype, np.integer):
        kind = "Int64"
    else:
        kind = "Float64"
    out = [f'        <DataArray type="{kind}" Name="{name}" '
           f'NumberOfComponents="{comp}" format="ascii">']
    out.append("          " + " ".join(_fmt(v) for v in arr.ravel()))
    out.append('        </DataArray>')
    return "\n".join(out)


def _polydata_xml(piece_attrs, pts, conn, offs, kind, cell_arrays):
    out = ['<?xml version="1.0"?>',
           '<VTKFile type="PolyData" version="1.0" byte_order="LittleEndian" header_type="UInt64">',
           '  <PolyData>',
           f'    <Piece {piece_attrs}>',
           '      <Points>',
           '        <DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">',
           "          " + " ".join(_fmt(v) for v in pts.ravel()),
           '        </DataArray>',
           '      </Points>',
           f'      <{kind}>',
           '        <DataArray type="Int64" Name="connectivity" format="ascii">',
           "          " + " ".join(str(int(v)) for v in conn),
           '        </DataArray>',
           '        <DataArray type="Int64" Name="offsets" format="ascii">',
           "          " + " ".join(str(int(v)) for v in offs),
           '        </DataArray>',
           f'      </{kind}>']
    if cell_arrays:
        out.append('      <CellData>')
        for name, data in cell_arrays.items():
            out.append(_xml_array(name, data))
        out.append('      </CellData>')
    out += ['    </Piece>', '  </PolyData>', '</VTKFile>']
    return "\n".join(out)


def write_vtp_lines(path, pts, conn, offsets, cell_arrays):
    Path(path).write_text(_polydata_xml(
        f'NumberOfPoints="{len(pts)}" NumberOfLines="{len(offsets)}"',
        pts, conn, offsets, "Lines", cell_arrays))


def write_vtp_polys(path, pts, faces, cell_arrays=None):
    offs = np.arange(1, len(faces) + 1) * faces.shape[1]
    Path(path).write_text(_polydata_xml(
        f'NumberOfPoints="{len(pts)}" NumberOfPolys="{len(faces)}"',
        pts, faces.ravel(), offs, "Polys", cell_arrays))


def write_vtm(path, blocks):
    out = ['<?xml version="1.0"?>',
           '<VTKFile type="vtkMultiBlockDataSet" version="0.1" byte_order="LittleEndian" header_type="UInt64">',
           '  <vtkMultiBlockDataSet version="0.1">',
           '    <BlockDataSet name="merge_curves" nested="1">']
    for k, (bname, fname, vtktype) in enumerate(blocks):
        out.append(f'      <Block name="{bname}" vtktype="{vtktype}" file="{fname}" index="{k}"/>')
    out += ['    </BlockDataSet>', '  </vtkMultiBlockDataSet>', '</VTKFile>']
    Path(path).write_text("\n".join(out))


# ---------------------------------------------------------------- groups
def build_groups(n_curves, tested, parallel_pairs):
    """One colour per merge group, as the partition in the module docstring.
    Returns (groups, loser_pairs): `groups` = colour-partition entries
    (tag, curves, status, pair_distance, why); `loser_pairs` = whole-curve
    tests that lost their pair to an earlier attempt (not their own colour,
    still in the highlight block + legend)."""
    free = set(range(n_curves))
    groups, losers = [], {}

    for att in tested:
        ci, cj = int(att["curve_i"]), int(att["curve_j"])
        why = att.get("reason") or "unknown failure"
        pdist = float(att.get("min_dist") if att.get("min_dist") is not None else -1.0)
        entry = {"tag": att["pair"], "curves": [ci, cj], "status": 1,
                 "pair_distance": pdist, "why": why}
        if ci in free and cj in free:
            free.discard(ci); free.discard(cj)
            groups.append(entry)
        else:
            losers[att["pair"]] = entry

    for pp in parallel_pairs:
        ci, cj = int(pp["i"]), int(pp["j"])
        if ci in free and cj in free:
            free.discard(ci); free.discard(cj)
            groups.append({
                "tag": f"cand-{ci}-{cj}", "curves": [ci, cj], "status": 0,
                "pair_distance": float(pp["min_dist"]),
                "why": "parallel-pair candidate only: whole-curve merge not tested here",
            })

    for c in sorted(free):
        groups.append({
            "tag": f"solo-{c}", "curves": [int(c)], "status": 0,
            "pair_distance": -1.0,
            "why": "no parallel curve partner: single curve, nothing to merge",
        })
    return groups, losers


# ---------------------------------------------------------------- render
def render_png(outdir, polylines, groups, tested, camera, P_ctx, surf_faces):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    faces = P_ctx[surf_faces]

    fig = plt.figure(figsize=(13, 9))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_box_aspect(BOX_ASPECT)

    ax.add_collection3d(Poly3DCollection(
        faces, facecolors=(0.82, 0.84, 0.87, 0.10),
        edgecolors="none"))

    palette = PALETTE[:len(groups)]
    segs, cols = [], []
    for gi, g in enumerate(groups):
        for cidx in g["curves"]:
            pl = polylines[cidx]
            for i in range(len(pl) - 1):
                segs.append((pl[i], pl[i + 1])); cols.append(palette[gi])
    ax.add_collection3d(Line3DCollection(segs, colors=cols, linewidths=1.8))

    thick, tcol = [], []
    for att in tested:
        for cidx in (int(att["curve_i"]), int(att["curve_j"])):
            pl = polylines[cidx]
            for i in range(len(pl) - 1):
                thick.append((pl[i], pl[i + 1])); tcol.append("#222222")
    ax.add_collection3d(Line3DCollection(thick, colors=tcol,
                                         linewidths=4.0,
                                         linestyles=(0, (4, 2)),
                                         alpha=0.5, zorder=10))

    handles = [plt.Line2D([0], [0], color=palette[gi], lw=2,
                          label=f"grp {gi} {g['tag']}: curves {g['curves']}, "
                                f"{g['status']}")
               for gi, g in enumerate(groups)]
    handles.append(plt.Line2D([0], [0], color="#222222", lw=4,
                              linestyle=(0, (4, 2)),
                              label="black dashed: took part in a tested "
                                    "whole-curve merge (3 pairs, all failed)"))
    ax.legend(handles=handles, loc="upper left", fontsize=7, framealpha=0.6)

    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
    ax.set_xticks(np.arange(0, 12, 2.0))
    ax.set_yticks(np.arange(0, 4, 1.0))
    ax.set_zticks(np.arange(-4, 11, 2.0))
    ax.set_xlim(*DOMAIN["x"]); ax.set_ylim(*DOMAIN["y"]); ax.set_zlim(*DOMAIN["z"])
    ax.view_init(elev=camera["elev"], azim=camera["azim"],
                 roll=camera.get("roll", 0.0))
    ax.set_proj_type(camera.get("proj_type", "persp"))
    png = outdir / OUT_PNG
    fig.savefig(png, dpi=180)
    plt.close(fig)
    print(f"[merge_curves] wrote {png}", flush=True)


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", default="/work/canada/sweep/n100000/hex_n100000.ovm")
    ap.add_argument("--analysis", default="/work/canada/sweep/pair_probe3/curve_analysis.json")
    ap.add_argument("--merge-test", default="/work/canada/sweep/pair_probe3/curve_merge_test.json")
    ap.add_argument("--blocks-vtk", default="/work/canada/sweep/struct_beam/blocks.vtk")
    ap.add_argument("--outdir", default="/work/canada/sweep/pair_probe3/viz")
    ap.add_argument("--camera", default=None,
                    help="JSON dict of elev/azim/roll/proj_type to override "
                         "the default separable view")
    ap.add_argument("--no-render", action="store_true", help="skip merge_curves.png")
    ap.add_argument("--no-context", action="store_true",
                    help="skip the block-baseline context surfaces")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    analysis = json.loads(Path(args.analysis).read_text())
    merge_test = json.loads(Path(args.merge_test).read_text())

    # ---- recover the curve geometry exactly as the analysers did -----------
    P, H = load_mesh(args.mesh)
    print(f"[merge_curves] mesh {len(P)} pts, {len(H)} hexes", flush=True)
    f2h, e2h = bc.build_topology(H)
    sing = bc.singular_edges(H, P, f2h, e2h, interior_only=True)
    curves, parallel = vs.analyze_curves(P, H, f2h, e2h, set(sing))
    n_curves = len(curves)
    polylines = [curve_polyline(P, c) for c in curves]
    n_edges = [len(c) for c in curves]
    print(f"[merge_curves] {len(sing)} singular edges, {n_curves} curves, "
          f"{len(parallel)} parallel pairs", flush=True)

    ref_edges = {int(c["idx"]): int(c["n_edges"]) for c in analysis["curves"]}
    if (n_curves != analysis["curve_count"]
            or len(sing) != analysis["singular_edges_total"]
            or len(parallel) != analysis["parallel_curve_pairs_total"]
            or any(ref_edges.get(ci) != n_edges[ci] for ci in range(n_curves))):
        print(f"[merge_curves] WARNING: recomputed curve table disagrees with "
              f"{args.analysis} (curves {n_curves} vs {analysis['curve_count']}, "
              f"sing edges {len(sing)} vs {analysis['singular_edges_total']}, "
              f"parallel pairs {len(parallel)} vs "
              f"{analysis['parallel_curve_pairs_total']}, "
              f"n_edges {n_edges} vs analysis {ref_edges})", flush=True)

    # ---- merge groups ------------------------------------------------------
    tested = merge_test["attempts"]
    if not tested:
        raise SystemExit(f"no whole-curve attempts in {args.merge_test}")
    groups, losers = build_groups(n_curves, tested, parallel)
    curve_group = {c: gi for gi, g in enumerate(groups) for c in g["curves"]}
    partner_of = {}
    for g in groups:
        if len(g["curves"]) == 2:
            a, b = g["curves"]
            partner_of[a] = b; partner_of[b] = a
    covered = len(curve_group)
    if covered != n_curves:
        raise SystemExit(f"partition covers {covered} of {n_curves} curves")

    a_curve = np.arange(n_curves, dtype=np.int64)
    a_group = np.array([curve_group[ci] for ci in range(n_curves)], np.int64)
    a_partner = np.array([partner_of.get(ci, -1) for ci in range(n_curves)], np.int64)
    a_status = np.array([groups[curve_group[ci]]["status"] for ci in range(n_curves)], np.int64)
    a_dist = np.array([groups[curve_group[ci]]["pair_distance"] for ci in range(n_curves)], float)
    a_nedges = np.array(n_edges, np.int64)
    cell_arrays = {"curve_id": a_curve, "group_id": a_group,
                   "merge_partner_curve_id": a_partner, "status": a_status,
                   "pair_distance": a_dist, "n_edges": a_nedges}

    # ---- write files -------------------------------------------------------
    concat_pts = np.concatenate(polylines)
    pts_per_curve = [len(p) for p in polylines]          # n_edges + 1
    conn = np.concatenate([base + np.arange(L) for base, L in
                           zip(np.r_[0, np.cumsum(pts_per_curve[:-1])], pts_per_curve)])
    offs = np.cumsum(pts_per_curve).astype(np.int64)     # offsets are exclusive end
    write_vtp_lines(outdir / OUT_VTP, concat_pts, conn, offs, cell_arrays)
    write_vtp_lines(outdir / OUT_BLOCKS["curves"], concat_pts, conn, offs, cell_arrays)
    print(f"[merge_curves] wrote {outdir / OUT_VTP} and block "
          f"{OUT_BLOCKS['curves']} ({n_curves} line cells, {len(concat_pts)} pts)",
          flush=True)

    P_ctx = surf_faces = None
    if not args.no_context:
        P_ctx, surf_faces, ctx_bid = beam_boundary_faces(args.blocks_vtk)
        write_vtp_polys(outdir / OUT_BLOCKS["context"], P_ctx, surf_faces,
                        {"block_id": ctx_bid})
        print(f"[merge_curves] wrote block {OUT_BLOCKS['context']} "
              f"({len(surf_faces)} boundary quads, "
              f"{len(set(int(x) for x in ctx_bid))} blocks)", flush=True)

    pts_t, conn_t, off_t, curve_t, tag_t = [], [], [], [], []
    base = 0
    for att in tested:
        for cidx in (int(att["curve_i"]), int(att["curve_j"])):
            pl = polylines[cidx]
            pts_t.append(pl)
            curve_t.append(cidx); tag_t.append(att["pair"])
            conn_t.extend(range(base, base + len(pl))); base += len(pl)
            off_t.append(base)
    write_vtp_lines(outdir / OUT_BLOCKS["tested_pairs"], np.concatenate(pts_t),
                    np.asarray(conn_t, np.int64), np.asarray(off_t, np.int64),
                    {"curve_id": np.asarray(curve_t, np.int64),
                     "pair_tag": np.array(tag_t, dtype="U8")})
    print(f"[merge_curves] wrote block {OUT_BLOCKS['tested_pairs']} "
          f"(tested pairs {[t['pair'] for t in tested]})", flush=True)

    blocks = [("curves_colored_by_merge_group", OUT_BLOCKS["curves"], "vtkPolyData")]
    if P_ctx is not None:
        blocks.append(("beam_block_context", OUT_BLOCKS["context"], "vtkPolyData"))
    blocks.append(("tested_pairs_highlight", OUT_BLOCKS["tested_pairs"], "vtkPolyData"))
    write_vtm(outdir / OUT_VTM, blocks)
    print(f"[merge_curves] wrote {outdir / OUT_VTM}", flush=True)

    # ---- legend ------------------------------------------------------------
    camera = dict(CAMERA)
    if args.camera:
        camera.update(json.loads(args.camera))
    palette = PALETTE[:len(groups)]
    legend = {
        "study": "canada singular-curve whole-merge tests (pair_probe3)",
        "mesh": analysis["mesh"],
        "sources": {"curve_analysis": args.analysis,
                    "curve_merge_test": args.merge_test,
                    "block_baseline_context": args.blocks_vtk},
        "outputs": {
            "curve_cells_vtp": OUT_VTP,
            "multiblock_vtm": OUT_VTM,
            "block_curves_vtp": OUT_BLOCKS["curves"],
            "block_context_vtp": OUT_BLOCKS["context"] if P_ctx is not None else "(skipped)",
            "block_tested_pairs_vtp": OUT_BLOCKS["tested_pairs"],
            "render_png": OUT_PNG},
        "counts": {
            "singular_edges": int(len(sing)),
            "curves": int(n_curves),
            "curve_line_cells": int(n_curves),
            "parallel_curve_pairs_in_analysis": int(len(parallel)),
            "whole_curve_merge_attempts": len(tested),
            "whole_curve_merges_passed": int(sum(1 for a in tested if a["ok"]))},
        "colour_rule": "one group_id per group; two curves share group/colour "
                       "iff they would be merged into one (applied as a "
                       "partition because a curve can hold only one colour: "
                       "tested pairs claim their curves in attempt order 1-5, "
                       "0-2, then 0-5 only if still free; remaining curves "
                       "pair greedily along the analysis parallel list; "
                       "leftovers are singleton groups)",
        "status_map": STATUS_NAME,
        "groups": [
            {
                "group_id": gi,
                "colour": palette[gi],
                "curves": g["curves"],
                "status": STATUS_NAME[g["status"]],
                "status_code": g["status"],
                "partner_of": ({str(a): b for a, b in zip(g["curves"], reversed(g["curves"]))}
                               if len(g["curves"]) == 2 else {}),
                "pair_distance": g["pair_distance"],
                "n_edges": {str(c): int(n_edges[c]) for c in g["curves"]},
                "tag": g["tag"],
                "what_failed_and_why": g["why"],
            } for gi, g in enumerate(groups)],
        "whole_curve_merge_attempts": [
            {"pair": a["pair"], "curves": [int(a["curve_i"]), int(a["curve_j"])],
             "own_colour_group": bool(partner_of.get(int(a["curve_i"]), -2)
                                      == int(a["curve_j"])),
             "ok": bool(a["ok"]),
             "reason": a.get("reason"),
             "hausdorff_max": (a.get("boundary") or {}).get("hausdorff_max"),
             "boundary_tolerance": (a.get("boundary") or {}).get("tol"),
             "boundary_vertex_count_change": (a.get("boundary") or {}).get(
                 "boundary_vertex_count_change"),
             "new_inverted_sj": (a.get("inversion") or {}).get("new_inverted_sj"),
             "tangent_deviation_deg": a.get("tangent_deviation_deg"),
             "min_dist": a.get("min_dist")}
            for a in tested],
        "verdict": merge_test.get("verdict"),
        "wall": merge_test.get("wall"),
        "camera_of_png": {**camera, "box_aspect": BOX_ASPECT, "domain": DOMAIN,
                          "note": "matplotlib Axes3D: ax.set_box_aspect("
                                  "camera['box_aspect']); ax.view_init(elev="
                                  "camera['elev'], azim=camera['azim'], roll="
                                  "camera['roll']); ax.set_proj_type(camera"
                                  "['proj_type']); ax.set_xlim/ylim/zlim("
                                  "camera['domain']['x'/'y'/'z'])"},
    }
    (outdir / OUT_LEGEND).write_text(json.dumps(legend, indent=2))
    print(f"[merge_curves] wrote {outdir / OUT_LEGEND}", flush=True)

    # ---- PNG ----------------------------------------------------------------
    if P_ctx is not None and not args.no_render:
        render_png(outdir, polylines, groups, tested, camera, P_ctx, surf_faces)


if __name__ == "__main__":
    main()
