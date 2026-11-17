"""2D quad layouts: rectilinear abstraction + separatrices, plane and periodic, with the fine Winslow mesh."""
import numpy as np
import pytest

from hexdp.core.layout2d import Domain2D, LayoutError, Loop, best_shear, layout, mesh, sheared


def sq(a, b, ccw=True):
    return np.array([(a, a), (b, a), (b, b), (a, b)], float)[:: 1 if ccw else -1]


def band(P, y0, y1, n=8):
    bot = Loop(np.c_[np.linspace(0, P, n, endpoint=False), np.full(n, y0)], shift=P)
    top = Loop(np.c_[np.linspace(P, 0, n, endpoint=False), np.full(n, y1)], shift=-P)
    return bot, top


def test_plane_domains():
    t = np.linspace(0, 2 * np.pi, 64, endpoint=False)
    disc = Domain2D([Loop(np.c_[np.cos(t), np.sin(t)])])
    assert layout(disc).n_blocks == 1                       # four virtual corners, no seed vertex
    assert mesh(layout(disc), disc, 0.1).min_sj > 0
    ring = Domain2D([Loop(sq(0, 3)), Loop(sq(1, 2, ccw=False))])
    assert layout(ring).n_blocks == 8
    assert mesh(layout(ring), ring, 0.25).min_sj == pytest.approx(1.0)
    L = np.array([(0, 0), (2, 0), (2, 1), (1, 1), (1, 2), (0, 2)], float)
    assert layout(Domain2D([Loop(L)])).n_blocks == 3


def test_pins_become_vertices():
    d = Domain2D([Loop(np.array([(0, 0), (1, 0), (2, 0), (2, 1), (0, 1)], float), pins=[1])])
    lay = layout(d)
    assert lay.n_blocks == 2                                # the pin on the bottom edge sends a ray up


def test_periodic_band_with_a_hole():
    P = 4.0
    d = Domain2D([*band(P, 0, 3), Loop(sq(1, 2, ccw=False))], period=P)
    lay = layout(d)
    assert lay.n_blocks == 5                                # 3 rows x 2 columns around the cylinder - the hole
    assert max(b[1] for b in lay.blocks) > len(lay.X) - 1    # one column crosses the seam
    m = mesh(lay, d, 0.25)
    assert m.min_sj == pytest.approx(1.0)


def test_periodic_band_without_vertices_is_refused():
    with pytest.raises(LayoutError):
        layout(Domain2D([*band(4.0, 0, 1)], period=4.0))


def test_best_shear_aligns_a_slanted_parallelogram():
    s = 0.6
    para = np.array([(0, 0), (1, 0), (1 + s, 1), (s, 1)], float)
    assert best_shear(Domain2D([Loop(para)])) == pytest.approx(s, abs=0.03)
    assert layout(sheared(Domain2D([Loop(para)]), s)).n_blocks == 1
