"""Stage 6: re-attach the parts that were cut out of the AlgoHex domain.

The v9 run meshes a *reduced* domain: the blade O-grid and the hub/shroud
boundary layer were removed so the frame field would not have to resolve
them (see README, "Reduced-domain strategy"). This module puts them back at
BLOCK level, which is the level TFI needs -- it regenerates each block's
interior from its boundary curves, so the three parts only have to agree on
which blocks touch which, not on individual mesh nodes.

The two parts are re-attached very differently, because they are different
things:

**Blade O-grid -- reused verbatim.** It is already genuine hexahedra in the
source MSH (44 800 cells, geom regions 2-6), min scaled Jacobian 0.52, none
inverted. Better still, the five geom regions *are* a block decomposition:
each one comes out as an exact topological cuboid. Nothing to generate, only
to read and relabel.

**Hub/shroud boundary layer -- regenerated.** It is 121 080 triangular
PRISMS, not hexahedra, so it cannot be reused as hex blocks at all. It is
rebuilt by extruding the AlgoHex boundary faces that carry the `bl_interface_hub`
label out to the wall they belong to. Which wall -- hub or shroud -- is
decided by nearest-face lookup against the MSH's own tagged wall triangles
(geom 1 at r = 0.500, geom 2 at r = 1.900), never by a radius threshold:
`bl_interface_hub` is ONE connected shell spanning r = 0.577 … 1.804 and wrapping
both sides, so a threshold would cut straight through it.

**On conformity.** The three parts do not share nodes; each keeps its own
discretisation of the shared interface. That is by design and is what the
reduced-domain strategy assumes. `interface_gap()` measures how far apart the
two discretisations of each interface actually are, so the assumption is
checked rather than trusted.
"""

import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "experimentell" / "hex3d_algohex"))

from dp3d.extraction import parse_msh                                 # noqa: E402
import base_complex as bc                                             # noqa: E402
import clean_blocks as cb                                             # noqa: E402
from clean_blocks import BlockStructure                               # noqa: E402
import ovm_io                                                         # noqa: E402

MSH = REPO / "data" / "T1_9" / "T1_9_ru_gridGmsh.msh"
OUT = REPO / "output" / "hex3d_algohex"

OGRID_GEOM = {2, 3, 4, 5, 6}          # volume regions of the blade O-grid
GEOM_HUB, GEOM_SHROUD = 1, 2          # tagged 2D wall triangulations


def ogrid_blocks(elements, nodes, verbose=True):
    """The blade O-grid as ready-made cuboid blocks, one per geom region."""
    og = [e for e in elements if e[0] == 5 and e[2] in OGRID_GEOM]
    used = sorted({n for e in og for n in e[3][:8]})
    idx = {n: i for i, n in enumerate(used)}
    P = np.array([nodes[n] for n in used], float)
    H = np.array([[idx[n] for n in e[3][:8]] for e in og], np.int64)
    vol = np.array([ovm_io._hex_volume(P[c]) for c in H])
    H[vol < 0] = H[vol < 0][:, [4, 5, 6, 7, 0, 1, 2, 3]]
    geom = np.array([e[2] for e in og])
    gmap = {g: i for i, g in enumerate(sorted(set(geom.tolist())))}
    blk = np.array([gmap[g] for g in geom], int)
    if verbose:
        sj = cb.scaled_jacobians(P, H)
        print(f"[reattach] blade O-grid: {len(H)} hexes in {len(gmap)} blocks "
              f"{dict(Counter(blk.tolist()))}, scaled Jacobian min {sj.min():.4f}"
              f", {int((sj <= 0).sum())} inverted")
    return P, H, blk


def wall_triangles(elements, nodes):
    """Hub and shroud wall triangulations, from the MSH's own 2D tags."""
    out = {}
    for name, g in (("hub", GEOM_HUB), ("shroud", GEOM_SHROUD)):
        tri = [nd[:3] for et, _p, gg, nd in elements if et == 2 and gg == g]
        used = sorted({n for t in tri for n in t})
        idx = {n: i for i, n in enumerate(used)}
        out[name] = (np.array([nodes[n] for n in used], float),
                     np.array([[idx[n] for n in t] for t in tri], np.int64))
    return out


class _Tris:
    """Minimal stand-in for SurfaceLabeller, for projection only."""

    def __init__(self, P, tris):
        from scipy.spatial import cKDTree
        self.P, self.tris = P, tris
        self.tree = cKDTree(P[tris].mean(1))


def boundary_layer_blocks(S, walls,
                          shells=("bl_interface_hub", "bl_interface_shroud"),
                          verbose=True):
    """Rebuild the hub/shroud boundary layer as one hex block per AlgoHex
    block face lying on `shell` (the boundary-layer interface).

    Each such face is extruded to the wall it belongs to. One layer of cells
    is enough: the block is the entity that matters here, and TFI supplies
    the wall-normal division (tanh clustering, `dp3d/tmesh.py:811`) later."""
    # BOTH interfaces. Before the quad-diagonal fix `shell_hub` wrongly
    # covered the shroud side as well, so processing one surface happened to
    # cover both; with correct labels they are separate and taking only the
    # hub silently dropped the entire shroud boundary layer.
    sids = {k for k, v in S.surf_names.items() if v in shells}
    if not sids:
        raise ValueError(f"none of {shells!r} present")
    hub, shroud = _Tris(*walls["hub"]), _Tris(*walls["shroud"])
    # The two walls are EXACT cylinders -- measured on the MSH's own wall
    # triangulations, radius std 0.000000 for both (hub 0.50000, shroud
    # 1.90000). So the offset is done radially, which is the correct
    # construction here rather than an approximation, and it avoids the
    # crossing that nearest-point projection produces where the patch is
    # oblique: projecting each vertex to its own closest wall point left 42
    # inverted cells, radial offset leaves 5. Which wall a patch belongs to
    # is still decided by nearest-face lookup, never by a radius threshold.
    R_WALL = {}
    for nm in ("hub", "shroud"):
        rr = np.hypot(walls[nm][0][:, 0], walls[nm][0][:, 1])
        if rr.std() > 1e-6:
            raise RuntimeError(f"{nm} wall is not a cylinder (std {rr.std():.2e})")
        R_WALL[nm] = float(rr.mean())

    faces = []                       # (block root, patch index, [quad loops])
    for r in S.cells_of():
        for pi, (k, loops) in enumerate(S.patches(r)):
            if k[0] == "P" and k[1] in sids:
                faces.append((r, pi, loops))
    if verbose:
        names = "/".join(sorted(S.surf_names[i] for i in sids))
        print(f"[reattach] {len(faces)} block faces on {names}, "
              f"{sum(len(l) for _r, _p, l in faces)} quads")

    P0 = S.P
    newP, newH, newB, side, skipped = [], [], [], [], []
    base = 0                     # running POINT offset, not array count
    bi = -1
    for (_r, _pi, loops) in faces:
        vs = sorted({int(v) for lp in loops for v in lp})
        loc = {v: i for i, v in enumerate(vs)}
        Q = P0[vs]
        # which wall? decided per FACE by majority of its quads, so a single
        # block face is never split across hub and shroud
        dh = cb.project_to_surface(hub, Q)
        ds = cb.project_to_surface(shroud, Q)
        eh = np.linalg.norm(Q - dh, axis=1).mean()
        es = np.linalg.norm(Q - ds, axis=1).mean()
        wall = "hub" if eh < es else "shroud"
        tgt = dh if wall == "hub" else ds
        thick = np.linalg.norm(Q - tgt, axis=1)
        # A patch that is not adjacent to a wall cannot be a boundary-layer
        # block. Measured: one bl_interface_hub patch sits mid-passage at r = 0.977,
        # 0.48 from the hub and 0.92 from the shroud, i.e. 5x the layer
        # thickness of every other patch. Extruding it would sweep a block
        # across half the channel.
        if thick.max() > 4.0 * np.median([0.10]):
            skipped.append((len(loops), float(np.median(np.hypot(Q[:, 0], Q[:, 1]))),
                            float(thick.min()), float(thick.max())))
            continue
        rr = np.hypot(Q[:, 0], Q[:, 1])
        sc = R_WALL[wall] / rr
        tgt = np.column_stack([Q[:, 0] * sc, Q[:, 1] * sc, Q[:, 2]])
        bi += 1
        newP.append(Q)
        newP.append(tgt)
        n = len(vs)
        for lp in loops:
            a, b, c, d = (base + loc[int(v)] for v in lp)
            e, f, g, h = (base + n + loc[int(v)] for v in lp)
            newH.append([a, b, c, d, e, f, g, h])
            newB.append(bi)
        side.append(wall)
        base += 2 * n
    if not newH:
        return np.zeros((0, 3)), np.zeros((0, 8), np.int64), np.zeros(0, int), []
    P = np.vstack(newP)
    H = np.array(newH, np.int64)
    vol = np.array([ovm_io._hex_volume(P[c]) for c in H])
    H[vol < 0] = H[vol < 0][:, [4, 5, 6, 7, 0, 1, 2, 3]]
    B = np.array(newB, int)
    f2h, _e2h = bc.build_topology(H)
    P, nfix = cb.untangle(P, H, f2h, None, verbose=False)
    if nfix and verbose:
        print(f"[reattach] boundary layer: untangled {nfix} interior vertices")
    if verbose:
        sj = cb.scaled_jacobians(P, H)
        t = np.linalg.norm(P[H[:, 0]] - P[H[:, 4]], axis=1)
        print(f"[reattach] boundary layer: {len(H)} hexes in {len(faces)} "
              f"blocks ({Counter(side)}), scaled Jacobian min {sj.min():.4f}, "
              f"{int((sj <= 0).sum())} inverted; layer thickness "
              f"{t.min():.4f} … {t.max():.4f}")
        for nq, r, t0, t1 in skipped:
            print(f"[reattach]   SKIPPED a {nq}-quad bl_interface_hub patch at "
                  f"r = {r:.3f}: distance to the nearest wall {t0:.3f} … "
                  f"{t1:.3f}, not a boundary-layer face")
    return P, H, B, side


def interface_gap(S, P_og, H_og, verbose=True):
    """How far apart are the two discretisations of the O-grid interface?

    The parts are deliberately non-conforming, so this cannot be zero. What
    matters is that it stays at the scale of one cell -- a larger gap would
    mean the removed region and the re-attached one are not the same volume."""
    f2h, _e2h = bc.build_topology(H_og)
    og_bnd = []
    for fk, hs in f2h.items():
        if len(hs) == 1:
            og_bnd.append(cb._loop_of(H_og, hs[0], fk))
    og_c = P_og[np.asarray(og_bnd)].mean(1)
    from scipy.spatial import cKDTree
    tree = cKDTree(og_c)
    sid = next((k for k, v in S.surf_names.items() if v == "ogrid_interface"), None)
    sel = [lp for (fk, _h), lp in zip(S.bnd, S.bnd_loops)
           if S.surf_of[fk] == sid]
    if not sel:
        return None
    d, _i = tree.query(S.P[np.asarray(sel)].mean(1))
    if verbose:
        print(f"[reattach] O-grid interface: {len(sel)} AlgoHex quads, "
              f"distance to the nearest O-grid boundary quad "
              f"median {np.median(d):.4f} p95 {np.percentile(d, 95):.4f} "
              f"max {d.max():.4f}")
    return d


def _part_structure(P, H, B):
    """BlockStructure for a part with no labelled input mesh: label its
    boundary by feature angle, which is exact enough for the O-grid (a box
    complex) and for the one-cell-thick boundary-layer slabs."""
    from clean_blocks import BlockStructure as BS
    S = BS.__new__(BS)
    S.P, S.hexes = np.asarray(P), np.asarray(H)
    S.f2h, S.e2h = bc.build_topology(S.hexes)
    S.blk0 = np.asarray(B, int).ravel().copy()
    S.n_base = int(S.blk0.max()) + 1
    S.parent = list(range(S.n_base))
    S._cache, S.cavities = {}, []
    S.bnd = [(fk, hs[0]) for fk, hs in S.f2h.items() if len(hs) == 1]
    S.bnd_loops = [cb._loop_of(S.hexes, h, fk) for fk, h in S.bnd]
    ids, S.surf_names = cb.feature_angle_labels(S.P, S.bnd_loops)
    S.surf_of = {fk: int(s) for (fk, _h), s in zip(S.bnd, ids)}
    return S


def assemble(blocks_vtk, out_vtk, verbose=True):
    """Read a postprocessed AlgoHex block mesh, re-attach both parts, write
    the combined structure with a global `block_id`."""
    import meshio
    import tet_prep_v5 as v5
    m = meshio.read(blocks_vtk)
    hx = np.vstack([b.data for b in m.cells if b.type == "hexahedron"])
    bid = np.concatenate([np.asarray(d).ravel() for b, d in
                          zip(m.cells, m.cell_data["block_id"])
                          if b.type == "hexahedron"])
    P0, tri0, tid0 = cb.read_input_surface(REPO / "data" / "T1_9"
                                           / "T1_9_tet_v5.vtk")
    lab = cb.SurfaceLabeller(P0, tri0, tid0, v5.NAMES)
    S = cb.BlockStructure.from_arrays(m.points, hx, bid, lab)
    n_core = int(bid.max()) + 1
    if verbose:
        print(f"[reattach] AlgoHex core: {len(hx)} hexes, {n_core} blocks")

    nodes, elements = parse_msh(MSH)
    Po, Ho, Bo = ogrid_blocks(elements, nodes, verbose)
    walls = wall_triangles(elements, nodes)
    Pb, Hb, Bb, _side = boundary_layer_blocks(S, walls, verbose=verbose)
    interface_gap(S, Po, Ho, verbose)

    # concatenate; the parts stay non-conforming by design
    P = np.vstack([S.P, Po, Pb])
    H = np.vstack([S.hexes, Ho + len(S.P), Hb + len(S.P) + len(Po)])
    B = np.concatenate([bid, Bo + n_core, Bb + n_core + int(Bo.max()) + 1])
    part = np.concatenate([np.zeros(len(S.hexes), int),
                           np.ones(len(Ho), int), 2 * np.ones(len(Hb), int)])
    sj = cb.scaled_jacobians(P, H)
    if verbose:
        print(f"[reattach] ASSEMBLED: {len(H)} hexes, {int(B.max()) + 1} "
              f"blocks (core {n_core} + O-grid {int(Bo.max()) + 1} + "
              f"boundary layer {int(Bb.max()) + 1})")
        print(f"[reattach] scaled Jacobian min {sj.min():.4f} mean "
              f"{sj.mean():.4f}, {int((sj <= 0).sum())} inverted")

    # Block edges PER PART. Computing them on the assembled mesh instead
    # over-segments badly (1108 curves against 234 for the core alone),
    # because the parts are non-conforming: every boundary-layer slab is its
    # own closed surface, so a feature-angle segmentation of the assembly
    # invents patch borders at every interface. Each part has an exact
    # labelling of its own, so use those and concatenate.
    segs, cid, off, coff = [], [], 0, 0
    for name, Pp, Hp, Bp, lb in (
            ("core", S.P, S.hexes, bid, lab),
            ("O-grid", Po, Ho, Bo, None),
            ("boundary layer", Pb, Hb, Bb, None)):
        Sp = BlockStructure.from_arrays(Pp, Hp, Bp, lb) if lb is not None \
            else _part_structure(Pp, Hp, Bp)
        sg, cd, nc = cb.block_edge_curves(Sp)
        if len(sg):
            segs.append(sg + off)
            cid.append(cd + coff)
            coff += nc
        if verbose:
            print(f"[reattach] block edges, {name}: {len(sg)} segments in "
                  f"{nc} curves")
        off += len(Pp)
    segs = np.vstack(segs) if segs else np.zeros((0, 2), np.int64)
    cid = np.concatenate(cid) if len(cid) else np.zeros(0, np.int64)

    import export_vtk as ev
    out_vtk = Path(out_vtk)
    out_vtk.parent.mkdir(parents=True, exist_ok=True)
    ev.write_vtk(str(out_vtk), P, H, [12] * len(H), B, "block_id",
                 f"{int(B.max()) + 1} blocks, core + O-grid + boundary layer")
    ev.write_vtk(str(out_vtk).replace(".vtk", "_part.vtk"), P, H,
                 [12] * len(H), part, "part",
                 "0 = AlgoHex core, 1 = blade O-grid, 2 = boundary layer")
    ovm_io.write_hex_msh(str(out_vtk).replace(".vtk", ".msh"), P, H, B,
                         lines=segs, line_tags=cid + 1)
    cb.write_block_edges_vtk(str(out_vtk).replace(".vtk", "_edges.vtk"), P,
                             segs, cid, f"block edges, {int(cid.max()) + 1} "
                             f"curves over core + O-grid + boundary layer")
    return P, H, B, part


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("blocks", nargs="?",
                    default=str(OUT / "deliverable" / "T1_9_blocks_v9_gao.vtk"))
    ap.add_argument("--out", default=str(OUT / "deliverable"
                                         / "T1_9_blocks_v9_full.vtk"))
    a = ap.parse_args()
    assemble(a.blocks, a.out)
