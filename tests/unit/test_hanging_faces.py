"""1:n hanging block faces (B12): detection on the fine index faces and the index-plane repair.

A block face conforms when all of its fine quads border the same neighbour (or the boundary) and that
neighbour's opposite face is exactly the same quad set. Otherwise the two blocks touch only along a
line in the coarse model: the block boundary pinches (chi 3, t_junction warnings). ``hanging_faces``
finds that; ``resolve_hanging_faces`` splits at the index planes where the neighbour changes.
"""
import os
from pathlib import Path

import numpy as np
import pytest

from hexdp.core import clean_blocks as cb
from hexdp.core.blockstructure import grid_cells
from hexdp.io.blocks_vtk import from_grids, structured_grids


class ZeroLabeller:
    """Every boundary face on surface 0; the hanging test does not look at labels."""
    names = {}
    max_dist = 0.0

    def label(self, Q):
        return np.zeros(len(Q), int)


def grid_mesh(nx, ny, nz):
    """Box of nx*ny*nz fine hexes on an integer lattice, VTK hex ordering (like ``grid_cells``)."""
    def node(i, j, k):
        return (i * (ny + 1) + j) * (nz + 1) + k

    P = np.array([(i, j, k) for i in range(nx + 1) for j in range(ny + 1) for k in range(nz + 1)],
                 float)
    H = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                H.append([node(i, j, k), node(i + 1, j, k), node(i + 1, j + 1, k), node(i, j + 1, k),
                          node(i, j, k + 1), node(i + 1, j, k + 1), node(i + 1, j + 1, k + 1),
                          node(i, j + 1, k + 1)])
    return P, np.array(H, np.int64)


def structure(dims, block_of):
    """BlockStructure over the box ``dims``; ``block_of(i, j, k)`` gives the block of a fine cell."""
    nx, ny, nz = dims
    P, H = grid_mesh(nx, ny, nz)
    bid = np.array([block_of(i, j, k)
                    for i in range(nx) for j in range(ny) for k in range(nz)], np.int64)
    return cb.BlockStructure.from_arrays(P, H, bid, ZeroLabeller())


def coarse(S):
    """The coarse model of the current partition (the read_blocks_vtk route)."""
    cells = S.cells_of()
    bid = np.zeros(len(S.hexes), np.int64)
    for r, c in cells.items():
        bid[np.asarray(c, np.int64)] = r
    grids, _ = structured_grids(S.P, S.hexes, bid)
    return grids, from_grids(S.P, grids)


def summarise(bs):
    kinds = sorted(i.kind for i in bs.check())
    return bs.boundary_topology()["chi"], bs.boundary_topology()["components"], kinds


def test_conforming_two_block_face_is_not_flagged():
    # x = 2 is a full 2x1 face on both sides: same quad set, same neighbour
    S = structure((4, 2, 1), lambda i, j, k: 0 if i < 2 else 1)
    assert cb.hanging_faces(S) == []


def test_one_two_hanging_face_is_flagged_and_split():
    # block 0 spans two fine quads on x = 2 and meets two blocks there (1 quad each): 1:2
    S = structure((4, 2, 1), lambda i, j, k: 0 if i < 2 else (1 if j < 1 else 2))
    bad = cb.hanging_faces(S)
    mixed = [h for h in bad if h["block"] == 0]
    assert len(mixed) == 1
    assert sorted(set(int(v) for v in mixed[0]["neighbours"].ravel())) == [1, 2]
    # the two small blocks see a uniform but too-large partner face; they are flagged too, and the
    # mixed face is the one the cut is read from
    assert {h["block"] for h in bad} == {0, 1, 2}
    assert cb._index_cuts(mixed[0]) == {(1, 1)}                # y = 1: 1 node index

    stats = cb.resolve_hanging_faces(S, verbose=False)
    assert stats == {"blocks_before": 3, "blocks_after": 4, "splits": 1, "rounds": 2,
                     "unresolved": 0, "stopped": "", "rolled_back": False}
    assert cb.hanging_faces(S) == []
    assert summarise(coarse(S)[1]) == (2, 1, [])


def test_one_three_hanging_face_is_flagged_and_split():
    # block 0 spans three fine quads on x = 1 and meets three blocks (1 quad each): 1:3
    S = structure((2, 3, 2), lambda i, j, k: 0 if i < 1 else 1 + j)
    bad = cb.hanging_faces(S)
    mixed = [h for h in bad if h["block"] == 0]
    assert len(mixed) == 1
    assert sorted(set(int(v) for v in mixed[0]["neighbours"].ravel())) == [1, 2, 3]

    stats = cb.resolve_hanging_faces(S, verbose=False)
    assert stats["blocks_before"] == 4 and stats["blocks_after"] == 6
    assert stats["unresolved"] == 0
    assert cb.hanging_faces(S) == []
    grids, B = coarse(S)
    assert sorted(int(np.prod(np.array(g.shape) - 1)) for g in grids) == [2] * 6
    assert summarise(B) == (2, 1, [])


def test_guard_rolls_back_instead_of_looping():
    # the same 1:2 defect with a floor no cut can satisfy: one round, then the original partition
    S = structure((4, 2, 1), lambda i, j, k: 0 if i < 2 else (1 if j < 1 else 2))
    before = np.array([S.root(S.blk0[h]) for h in range(len(S.hexes))])
    stats = cb.resolve_hanging_faces(S, min_cells=100, verbose=False)
    assert stats["rolled_back"] is True
    assert stats["blocks_after"] == 3 and stats["splits"] == 0
    assert stats["unresolved"] == len(cb.hanging_faces(S)) == 3
    assert np.array_equal(before, np.array([S.root(S.blk0[h]) for h in range(len(S.hexes))]))


def test_max_blocks_guard_also_rolls_back():
    S = structure((4, 2, 1), lambda i, j, k: 0 if i < 2 else (1 if j < 1 else 2))
    stats = cb.resolve_hanging_faces(S, max_blocks=3, verbose=False)
    assert stats["rolled_back"] is True and stats["blocks_after"] == 3
    assert stats["stopped"].startswith("guard")


# --- integration: the real s10b cut, from the stored pre-steps -------------------------------

STEPS = Path(os.environ.get("HEXDP_B12_STEPS", "/work/scratch/topo/s10b/steps"))


@pytest.mark.skipif(not (STEPS / "cut_non_cuboids_in.npz").exists(),
                    reason="s10b pre-steps not on this machine (HEXDP_B12_STEPS)")
def test_s10b_hanging_face_from_stored_steps():
    d = np.load(STEPS / "cut_non_cuboids_in.npz")
    S = cb.BlockStructure.from_arrays(d["P"], d["H"], d["bid"], ZeroLabeller())
    cb.cut_non_cuboids(S, verbose=False)
    assert len(S.cells_of()) == 4                       # the cut itself ends on the singular edge
    assert cb.hanging_faces(S)

    stats = cb.resolve_hanging_faces(S, verbose=False)
    assert stats["blocks_before"] == 4
    assert stats["blocks_after"] == 7                   # 4 -> 7, the cascade after the first cut
    assert stats["unresolved"] == 0 and not stats["rolled_back"]
    assert cb.hanging_faces(S) == []

    grids, B = coarse(S)
    assert summarise(B) == (2, 1, [])                   # chi 3 -> 2, no t_junction
    # the split is an index cut on the existing fine mesh: the fine cells did not move
    fine = [tuple(sorted(c)) for g in grids for c in grid_cells(g).tolist()]
    assert sorted(fine) == sorted(tuple(sorted(c)) for c in np.asarray(d["H"]).tolist())
