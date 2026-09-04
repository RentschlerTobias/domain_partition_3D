"""Stage 4b: turn base-complex components into structured (i,j,k) blocks.

A base-complex component is only useful for TFI if it is a topological
cuboid, i.e. its cells can be given consistent (i,j,k) grid coordinates.
This module decides that by a grid walk and returns the block dimensions.

The frame has to be propagated correctly from cell to cell. Each hex has 3
local axes (3 pairs of opposite faces / 3 groups of 4 parallel edges). When
stepping through a shared face into a neighbour, the neighbour's local axes
are matched to the current ones VIA THE EDGES OF THE SHARED FACE: an edge of
that face belongs to a definite local axis in each of the two hexes, which
fixes the correspondence. Guessing the mapping instead (e.g. assuming the
remaining axes keep their order) silently produces garbage -- it reported
0 of 5 cuboids on the cylinder O-grid, which is impossible.
"""

import sys
from collections import Counter, deque
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "experimentell" / "hex3d_algohex"))

import base_complex as bc                                             # noqa: E402
import ovm_io                                                         # noqa: E402

HF = bc.HEX_FACES
# Local axis a is the direction NORMAL to the face pair AX_FACES[a], so the
# parallel edges of axis a must be the edges running along that normal.
# AX_FACES[0] = bottom/top -> normal is the 0-4 direction, NOT 0-1: the two
# tables have to be aligned by construction, otherwise a face's own edges get
# classified as its normal axis and the frame propagation collapses.
AX_FACES = [(0, 1), (2, 4), (3, 5)]
AX_EDGES = [[(0, 4), (1, 5), (2, 6), (3, 7)],   # normal of bottom/top
            [(0, 3), (1, 2), (4, 7), (5, 6)],   # normal of (0,1,5,4)/(2,3,7,6)
            [(0, 1), (3, 2), (4, 5), (7, 6)]]   # normal of (1,2,6,5)/(3,0,4,7)


def hex_maps(cell):
    """face-key -> local axis; edge-key -> (local axis, orientation).

    The orientation records whether the stored edge runs along or against the
    axis' canonical direction, which is what makes the sign of the frame
    transfer between two cells well defined."""
    c = [int(v) for v in cell]
    fk, ek = {}, {}
    for a, (p, q) in enumerate(AX_FACES):
        for fi in (p, q):
            fk[frozenset(c[i] for i in HF[fi])] = a
    for a, grp in enumerate(AX_EDGES):
        for i, j in grp:
            u, v = c[i], c[j]
            ek[(min(u, v), max(u, v))] = (a, 1 if u < v else -1)
    return fk, ek


def face_edges(cell, fi):
    c = [int(v) for v in cell]
    loop = [c[i] for i in HF[fi]]
    return [(min(loop[k], loop[(k + 1) % 4]), max(loop[k], loop[(k + 1) % 4]))
            for k in range(4)]


def block_grid(hexes, blk, f2h, blk_of, bi):
    """Assign (i,j,k) to every cell of the block. Returns dims or None."""
    FK, EK = {}, {}
    for h in blk:
        FK[h], EK[h] = hex_maps(hexes[h])

    seed = blk[0]
    coord = {seed: (0, 0, 0)}
    # perm[h][local_axis] = (global_axis, sign)
    perm = {seed: {0: (0, 1), 1: (1, 1), 2: (2, 1)}}
    dq = deque([seed])

    while dq:
        h = dq.popleft()
        ci = coord[h]
        for a, (fp, fq) in enumerate(AX_FACES):
            g, s = perm[h][a]
            for fi, step in ((fp, +1), (fq, -1)):
                k = frozenset(int(hexes[h][i]) for i in HF[fi])
                nb = [x for x in f2h.get(k, ()) if x != h]
                if not nb or blk_of[nb[0]] != bi:
                    continue
                h2 = nb[0]
                d = [0, 0, 0]
                d[g] = s * step
                nc = (ci[0] + d[0], ci[1] + d[1], ci[2] + d[2])
                if h2 in coord:
                    if coord[h2] != nc:
                        return None
                    continue
                # match h2's local axes to the global ones
                a2 = FK[h2].get(k)
                if a2 is None:
                    return None
                p2 = {a2: (g, s)}
                for ed in face_edges(hexes[h], fi):
                    ra, rb = EK[h].get(ed), EK[h2].get(ed)
                    if ra is None or rb is None:
                        continue
                    la, oa = ra
                    lb, ob = rb
                    if la == a:
                        continue
                    gg, ss = perm[h][la]
                    # the same physical edge may run along the axis in h and
                    # against it in h2 -> transfer the sign accordingly
                    val = (gg, ss * oa * ob)
                    if lb in p2 and p2[lb] != val:
                        return None
                    p2[lb] = val
                if len(p2) != 3:
                    return None
                perm[h2] = p2
                coord[h2] = nc
                dq.append(h2)

    if len(coord) != len(blk):
        return None
    C = np.array([coord[h] for h in blk])
    C -= C.min(0)
    dims = tuple(int(x) for x in (C.max(0) + 1))
    if len(set(map(tuple, C))) != len(blk):
        return None
    if int(np.prod(dims)) != len(blk):
        return None
    return dims


def analyse(ovm_path):
    P, e, f, poly = ovm_io.read_ovm(ovm_path)
    hexes, _ = ovm_io.ovm_to_cells(P, e, f, poly)
    f2h, e2h = bc.build_topology(hexes)
    sing = bc.singular_edges(hexes, P, f2h, e2h)
    cut = bc.sheet_faces(hexes, f2h, sing)
    blocks = bc.blocks_from_cut(hexes, f2h, cut)
    blk_of = np.zeros(len(hexes), int)
    for i, b in enumerate(blocks):
        for h in b:
            blk_of[h] = i
    dims, bad = [], 0
    for bi, b in enumerate(blocks):
        d = block_grid(hexes, b, f2h, blk_of, bi)
        if d is None:
            bad += 1
        else:
            dims.append(d)
    print(f"[block_structure] {len(blocks)} blocks: "
          f"{len(dims)} structured cuboids, {bad} not cuboid "
          f"({100 * len(dims) / max(1, len(blocks)):.0f}% ok)")
    if dims:
        cells = [int(np.prod(d)) for d in dims]
        print(f"[block_structure]   cells covered by cuboids: {sum(cells)}"
              f" of {len(hexes)} ({100 * sum(cells) / len(hexes):.0f}%)")
        print(f"[block_structure]   largest dims: "
              f"{sorted(dims, key=lambda d: -np.prod(d))[:8]}")
    return P, hexes, blocks, dims


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("ovm", nargs="?",
                    default=str(REPO / "output" / "hex3d_algohex"
                               / "T1_9_hex_v5.ovm"))
    a = ap.parse_args()
    analyse(a.ovm)
