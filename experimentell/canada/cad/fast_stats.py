"""Vectorised drop-in for clean_blocks._structure_stats (same result, numpy + scipy.sparse.csgraph).

The original runs pure-Python loops over ~3 M face/edge entries: 35 s on a 106k-cell mesh, of which
label_sheets 15 s, build_topology 8 s, _face_edges_map 4 s, singular_edges 3 s, block patches 5 s.
_structure_stats only needs the SET of separatrix faces (blocks_from_cut gets set(sheet_of)), so the
sheet propagation is a reachability question:
  * faces / edges: unique sorted vertex tuples per hex (np.unique)
  * singular edges: interior edges (not on a boundary face) with hex valence != 4
  * "straight across" from face f over its edge g (base_complex.sheet_faces / block_faces.label_sheets):
    in the first hex h0 of f (ascending index) that works, the other face f1 of h0 containing g, the
    other hex h1 of f1, the other face f2 of h1 containing g (f2 != f). Per hex this is a fixed table.
  * separatrix faces = faces reachable from the seed faces (faces with a singular edge) along
    straight-across steps over non-singular edges -> BFS in a directed sparse graph
  * blocks = connected components of the hex graph over non-separatrix interior faces
The per-block patch part (cuboid status, excess faces) reuses clean_blocks' own functions.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import breadth_first_order, connected_components

HEX = Path("/opt/stack/domain_partition_3D/experimentell/hex3d_algohex")
sys.path.insert(0, str(HEX.parent.parent)); sys.path.insert(0, str(HEX))
import clean_blocks as cb                                     # noqa: E402
import base_complex as bc                                     # noqa: E402

HF = np.array(bc.HEX_FACES)          # (6, 4) local vertices of the faces
HE = np.array(bc.HEX_EDGES)          # (12, 2)
# local edge index of each face side, and for (face, side) the other face of the hex with that edge
_FACE_EDGE = np.zeros((6, 4), int)
for f in range(6):
    for k in range(4):
        a, b = HF[f, k], HF[f, (k + 1) % 4]
        _FACE_EDGE[f, k] = next(i for i, (u, v) in enumerate(HE) if {u, v} == {a, b})
_OTHER_FACE = np.zeros((6, 4), int)
for f in range(6):
    for k in range(4):
        e = _FACE_EDGE[f, k]
        _OTHER_FACE[f, k] = next(g for g in range(6) if g != f and e in _FACE_EDGE[g])
# side index of local edge e within face g
_SIDE = -np.ones((6, 12), int)
for f in range(6):
    for k in range(4):
        _SIDE[f, _FACE_EDGE[f, k]] = k


def topology(hexes):
    H = np.asarray(hexes, np.int64)
    n = len(H)
    N = np.int64(H.max() + 1)
    # integer keys instead of np.unique(axis=0) (2.6 s of 6.6 s on 106k cells)
    fv = np.sort(H[:, HF], axis=2).reshape(-1, 4)                     # (6n, 4)
    k1, k2 = fv[:, 0] * N + fv[:, 1], fv[:, 2] * N + fv[:, 3]
    o = np.lexsort((k2, k1))
    new = np.ones(len(o), bool)
    new[1:] = (k1[o][1:] != k1[o][:-1]) | (k2[o][1:] != k2[o][:-1])
    fid = np.empty(len(o), np.int64)
    fid[o] = np.cumsum(new) - 1
    fid = fid.reshape(n, 6)
    nf = int(new.sum())
    ev = np.sort(H[:, HE], axis=2).reshape(-1, 2)
    eu, eid = np.unique(ev[:, 0] * N + ev[:, 1], return_inverse=True)
    eid = eid.reshape(n, 12)
    # face -> hexes (ascending hex index, like the dict of lists in build_topology)
    hx = np.repeat(np.arange(n), 6)
    order = np.lexsort((hx, fid.ravel()))
    fcount = np.bincount(fid.ravel(), minlength=nf)
    start = np.concatenate([[0], np.cumsum(fcount)[:-1]])
    f_h0 = hx[order][start]
    f_h1 = np.where(fcount == 2, hx[order][np.minimum(start + 1, len(order) - 1)], -1)
    # local face index of a face inside each of its hexes
    loc = np.tile(np.arange(6), n)[order]
    f_l0 = loc[start]
    f_l1 = np.where(fcount == 2, loc[np.minimum(start + 1, len(order) - 1)], -1)
    return dict(H=H, fid=fid, nf=nf, eid=eid, ne=len(eu), eu=eu, fcount=fcount,
                f_h0=f_h0, f_h1=f_h1, f_l0=f_l0, f_l1=f_l1)


def singular(T):
    eid, fid, fcount = T["eid"], T["fid"], T["fcount"]
    val = np.bincount(eid.ravel(), minlength=T["ne"])
    bfaces = fcount == 1
    bnd = np.zeros(T["ne"], bool)
    hi, lf = np.nonzero(bfaces[fid])
    bnd[eid[hi[:, None], _FACE_EDGE[lf]].ravel()] = True
    return (~bnd) & (val != 4)


def separatrix_faces(T, sing_e):
    H, fid, eid, nf = T["H"], T["fid"], T["eid"], T["nf"]
    f_h = np.stack([T["f_h0"], T["f_h1"]], 1)
    f_l = np.stack([T["f_l0"], T["f_l1"]], 1)
    # face -> its 4 global edges (via its first hex)
    fe = eid[T["f_h0"][:, None], _FACE_EDGE[T["f_l0"]]]                 # (nf, 4)
    src, dst = [], []
    for k in range(4):
        g = fe[:, k]
        res = np.full(nf, -1)
        for slot in (0, 1):                     # first hex of the face that yields a result wins
            h0, l0 = f_h[:, slot], f_l[:, slot]
            ok = (h0 >= 0) & (res < 0)
            if not ok.any():
                continue
            idx = np.nonzero(ok)[0]
            h0, l0 = h0[idx], l0[idx]
            le = _FACE_EDGE[l0, _side_of(eid, h0, l0, g[idx])]          # local edge of g in h0
            lf1 = _OTHER_FACE[l0, _SIDE[l0, le]]
            f1 = fid[h0, lf1]
            # other hex of f1
            a, b = T["f_h0"][f1], T["f_h1"][f1]
            h1 = np.where(a == h0, b, a)
            la, lb = T["f_l0"][f1], T["f_l1"][f1]
            l1 = np.where(a == h0, lb, la)
            good = h1 >= 0
            if not good.any():
                continue
            idx, h1, l1, f1_ = idx[good], h1[good], l1[good], f1[good]
            le1 = _edge_local(eid, h1, g[idx])
            lf2 = _OTHER_FACE[l1, _SIDE[l1, le1]]
            f2 = fid[h1, lf2]
            keep = f2 != idx
            res[idx[keep]] = f2[keep]
        m = (res >= 0) & ~sing_e[g]
        src.append(np.nonzero(m)[0]); dst.append(res[m])
    src = np.concatenate(src); dst = np.concatenate(dst)
    seeds = np.nonzero(sing_e[fe].any(1))[0]
    # BFS from a virtual source connected to all seeds
    vs = nf
    rows = np.concatenate([src, np.full(len(seeds), vs)])
    cols = np.concatenate([dst, seeds])
    G = csr_matrix((np.ones(len(rows), np.int8), (rows, cols)), shape=(nf + 1, nf + 1))
    reach = breadth_first_order(G, vs, directed=True, return_predecessors=False)
    cut = np.zeros(nf, bool)
    cut[reach[reach < nf]] = True
    return cut


def _side_of(eid, h, lf, g):
    """side index (0..3) of global edge g in local face lf of hex h"""
    E = eid[h[:, None], _FACE_EDGE[lf]]                                # (m, 4)
    return np.argmax(E == g[:, None], axis=1)


def _edge_local(eid, h, g):
    return np.argmax(eid[h] == g[:, None], axis=1)


def blocks(T, cut):
    inner = (T["fcount"] == 2) & ~cut
    a, b = T["f_h0"][inner], T["f_h1"][inner]
    n = len(T["H"])
    G = csr_matrix((np.ones(len(a), np.int8), (a, b)), shape=(n, n))
    nb, lab = connected_components(G, directed=False)
    return nb, lab


def structure_stats(P, hexes, labeller):
    """same dict as clean_blocks._structure_stats"""
    T = topology(hexes)
    sing_e = singular(T)
    cut = separatrix_faces(T, sing_e)
    nb, blk = blocks(T, cut)
    H, fid = T["H"], T["fid"]
    # boundary faces: surface label of the face centroid
    bmask = T["fcount"][fid] == 1                                       # (n, 6)
    bh, bl = np.nonzero(bmask)
    loops = H[bh[:, None], HF[bl]]
    surf = labeller.label(P[loops].mean(1)) if len(loops) else np.zeros(0, int)
    # block-boundary faces with their patch keys, per block
    other = np.where(T["f_h0"][fid] == np.arange(len(H))[:, None], T["f_h1"][fid], T["f_h0"][fid])
    rec = [[] for _ in range(nb)]
    for h, l, s in zip(bh, bl, surf):
        rec[blk[h]].append((("P", int(s)), [int(v) for v in H[h, HF[l]]]))
    ih, il = np.nonzero((~bmask) & (blk[np.maximum(other, 0)] != blk[:, None]))
    for h, l in zip(ih, il):
        rec[blk[h]].append((("B", int(blk[other[h, l]])), [int(v) for v in H[h, HF[l]]]))
    cub = exc = 0
    for r in rec:
        p = cb._patch_split(r)
        cub += cb.cuboid_status(p)[0] == "cuboid"
        exc += cb._excess(p)
    sizes = np.bincount(blk, minlength=nb)
    sj = cb.scaled_jacobians(P, H)
    return {"cells": len(H), "blocks": int(nb), "cuboids": int(cub), "excess": int(exc),
            "sing": int(sing_e.sum()), "tiny": int((sizes < 10).sum()),
            "inverted": int((sj <= 0).sum()), "min_sj": float(sj.min()),
            "hausdorff": float(labeller.max_dist)}
