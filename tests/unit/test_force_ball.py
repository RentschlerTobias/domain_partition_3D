"""A hex mesh without interior singular edges around a hole (a square ring) is one non-ball base block; the
force-ball cut splits it into cuboid blocks."""
import numpy as np

from hexdp.core import base_complex as bc
from hexdp.core import block_faces as bfx
from hexdp.core import clean_blocks as cb


def ring_mesh(n=4, hole=(1, 3), nz=2):
    """n x n x nz hex grid with the cells hole[0]..hole[1]-1 in x and y removed (a square ring)."""
    idx = -np.ones((n + 1, n + 1, nz + 1), int)
    P = []
    for i in range(n + 1):
        for j in range(n + 1):
            for k in range(nz + 1):
                idx[i, j, k] = len(P)
                P.append((i, j, k))
    H = []
    for i in range(n):
        for j in range(n):
            if hole[0] <= i < hole[1] and hole[0] <= j < hole[1]:
                continue
            for k in range(nz):
                H.append([idx[i, j, k], idx[i + 1, j, k], idx[i + 1, j + 1, k], idx[i, j + 1, k],
                          idx[i, j, k + 1], idx[i + 1, j, k + 1], idx[i + 1, j + 1, k + 1], idx[i, j + 1, k + 1]])
    return np.array(P, float), np.array(H)


def test_ring_is_one_non_ball_block_and_gets_cut():
    P, H = ring_mesh()
    f2h, e2h = bc.build_topology(H)
    sing = bc.singular_edges(H, P, f2h, e2h)
    sheet_of, _ = bfx.label_sheets(H, f2h, sing)
    blocks = bc.blocks_from_cut(H, f2h, set(sheet_of))
    assert len(blocks) == 1 and cb.non_ball_blocks(H, blocks) == 1
    cut, added = cb.force_ball_cut(H, f2h, set(sheet_of), log=lambda *a: None)
    assert added > 0
    after = bc.blocks_from_cut(H, f2h, cut)
    assert len(after) >= 2 and cb.non_ball_blocks(H, after) == 0
    assert not any(cb._not_cuboid(H, b) for b in after)            # every block a structured lattice
    assert sum(len(b) for b in after) == len(H)
