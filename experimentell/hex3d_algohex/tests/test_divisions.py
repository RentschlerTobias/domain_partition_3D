"""Acceptance checks for the conforming-division MILP (tfi.solve_block_divisions).

Run: PY -m pytest experimentell/hex3d_algohex/tests/test_divisions.py -q
or simply: PY experimentell/hex3d_algohex/tests/test_divisions.py
"""

import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent.parent
sys.path.insert(0, str(REPO / "experimentell" / "hex3d_algohex"))

import tfi                                                            # noqa: E402

BLOCKS = REPO / "output" / "hex3d_algohex" / "deliverable" / "T1_9_blocks_v11.vtk"


def _load():
    P, H, B, f2h = tfi.load_blocks(str(BLOCKS))
    lat, missing = tfi.lattices(P, H, B, f2h, verbose=False)
    classes = tfi.direction_classes(lat, f2h, H, B, verbose=False)
    return P, H, B, f2h, lat, missing, classes


def test_counts_meet_their_lower_bound_and_conform():
    P, H, B, f2h, lat, missing, classes = _load()
    t0 = time.time()
    counts = tfi.solve_block_divisions(lat, classes, P, 0.05, verbose=False)
    dt = time.time() - t0
    assert dt < 60, f"solver took {dt:.1f}s"

    for ci, c in enumerate(classes):
        lb = max(max(1, int(round(tfi.axis_length(P, lat[r][1], ax) / 0.05)))
                 for r, ax in c)
        assert counts[ci] >= lb, f"class {ci}: {counts[ci]} < lower bound {lb}"

    # the counts, fed back in as if they were a mesh, must be conforming:
    # every axis of a class carries the same number
    new = tfi.block_counts(lat, classes, counts)
    synthetic = {r: (np.array(new[r]), lat[r][1]) for r in lat}
    assert tfi.check_conformity(synthetic, classes)


def test_coarser_than_the_domain_stays_feasible():
    """--target-h 10 is coarser than any axis: every count falls to 1 and the
    solver must say so rather than come back infeasible."""
    P, H, B, f2h, lat, missing, classes = _load()
    counts = tfi.solve_block_divisions(lat, classes, P, 10.0, verbose=False)
    assert set(counts.values()) == {1}, counts


def test_frozen_counts_are_honoured():
    """A pinned class keeps the neighbour's existing count even when the
    target h asks for something else."""
    P, H, B, f2h, lat, missing, classes = _load()
    r, ax = classes[0][0]
    pin = int(lat[r][0][ax]) + 3
    counts = tfi.solve_block_divisions(lat, classes, P, 0.05,
                                       frozen_counts={(r, ax): pin},
                                       verbose=False)
    assert counts[0] == pin, (counts[0], pin)


if __name__ == "__main__":
    for fn in (test_counts_meet_their_lower_bound_and_conform,
               test_coarser_than_the_domain_stays_feasible,
               test_frozen_counts_are_honoured):
        fn()
        print(f"PASS {fn.__name__}")
