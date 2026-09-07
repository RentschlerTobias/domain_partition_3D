"""Acceptance checks for the block refill (tfi.refill_block / refill_complex).

The fold gate of TFI_RESEARCH.md, at complex level: refilling at the counts a
block already has must give the boundary back to machine precision and must
leave the mesh watertight. If that fails, no prescribed-h refill is worth
looking at.

Run: PY -m pytest experimentell/hex3d_algohex/tests/test_refill.py -q
"""

import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent.parent
sys.path.insert(0, str(REPO / "experimentell" / "hex3d_algohex"))

import base_complex as bc                                              # noqa: E402
import clean_blocks as cb                                              # noqa: E402
import tfi                                                             # noqa: E402

BLOCKS = REPO / "output" / "hex3d_algohex" / "deliverable" / "T1_9_blocks_v11.vtk"


def _identity():
    P, H, B, f2h = tfi.load_blocks(str(BLOCKS))
    lat, _missing = tfi.lattices(P, H, B, f2h, verbose=False)
    classes = tfi.direction_classes(lat, f2h, H, B, verbose=False)
    counts = {ci: int(lat[classes[ci][0][0]][0][classes[ci][0][1]])
              for ci in range(len(classes))}
    out = tfi.refill_complex(P, H, B, f2h, lat, classes, counts, verbose=False)
    return P, H, B, f2h, out


def test_identity_refill_keeps_boundary_and_watertightness():
    P, H, B, f2h, (Pn, Hn, Bn, rep) = _identity()
    assert len(Hn) == len(H)

    ok, _bnd = tfi.check_watertight(Pn, Hn, verbose=False)
    assert ok, "refilled complex is not watertight"

    bnd_old = sorted({v for fk, hs in f2h.items() if len(hs) == 1 for v in fk})
    f2hn, _e = bc.build_topology(Hn)
    bnd_new = sorted({v for fk, hs in f2hn.items() if len(hs) == 1 for v in fk})
    assert len(bnd_new) == len(bnd_old)
    d, _ = cKDTree(Pn[bnd_new]).query(P[bnd_old])
    assert d.max() < 1e-12, f"boundary moved by {d.max():.2e}"

    sj = cb.scaled_jacobians(Pn, Hn)
    assert int((sj <= 0).sum()) == 0, "identity refill inverted cells"
    assert sj.min() >= 0.11, f"min scaled Jacobian {sj.min():.4f}"


def test_resample_is_orientation_invariant():
    """The weld's precondition: two blocks see a shared face transposed and
    mirrored relative to each other, and must still land on the same points."""
    rng = np.random.default_rng(0)
    Q = rng.normal(size=(7, 5, 3))
    fu, fv = tfi.uniform_fractions(9), tfi.uniform_fractions(4)
    A = tfi.resample_face_grid(Q, fu, fv)
    B1 = tfi.resample_face_grid(Q.transpose(1, 0, 2), fv, fu)
    assert np.allclose(A, B1.transpose(1, 0, 2), atol=1e-15)
    B2 = tfi.resample_face_grid(Q[::-1], fu[::-1], fv)
    assert np.allclose(A, B2, atol=1e-15)


def test_refill_block_hits_the_requested_dimensions():
    P, H, B, f2h = tfi.load_blocks(str(BLOCKS))
    lat, _m = tfi.lattices(P, H, B, f2h, verbose=False)
    r = sorted(lat)[0]
    X = tfi.refill_block(P, lat[r][1], (4, 5, 6))
    assert X.shape == (5, 6, 7, 3)


if __name__ == "__main__":
    for fn in (test_identity_refill_keeps_boundary_and_watertightness,
               test_resample_is_orientation_invariant,
               test_refill_block_hits_the_requested_dimensions):
        fn()
        print(f"PASS {fn.__name__}")
