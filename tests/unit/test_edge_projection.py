"""B14: boundary block edges are projected onto the surface, because the collapse leaves 2-5 nodes per edge."""
from types import SimpleNamespace

import numpy as np

from hexdp.core import clean_blocks as cb
from hexdp.export.dataset import (_boundary_corner_pairs, _boundary_triangles, _project_edge, _seam_curve,
                                  block_edges)
from hexdp.io.blocks_vtk import from_grids


def _cube(tmp_path=None):
    """One block, 2x2x2 fine cells; every fine node is a surface point."""
    P, index = [], {}

    def node(x, y, z):
        key = (x, y, z)
        if key not in index:
            index[key] = len(P)
            P.append(key)
        return index[key]

    G = np.empty((3, 3, 3), np.int64)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                G[i, j, k] = node(i, j, k)
    return from_grids(np.array(P, float), [G])


def test_the_boundary_edges_of_a_cube_are_its_twelve_edges():
    assert len(_boundary_corner_pairs(_cube())) == 12


def test_an_edge_that_would_move_too_far_is_left_alone():
    # a projection that moves a point more than half the chord probably crossed a gap
    q = np.array([[0.5, 0.0, 0.0]])
    far = lambda lab, Q: Q + np.array([[0.0, 0.0, 2.0]])
    a, b = np.zeros(3), np.array([1.0, 0, 0])
    assert np.allclose(_project_edge(far, None, q, a, b), q)
    near = lambda lab, Q: Q + np.array([[0.0, 0.0, 0.2]])
    one_patch = SimpleNamespace(label=lambda P: np.zeros(len(P), int))
    assert np.allclose(_project_edge(near, one_patch, q, a, b), [[0.5, 0.0, 0.2]])


def test_projection_onto_the_own_surface_does_not_move_a_straight_edge():
    bs = _cube()
    V, T, tl = _boundary_triangles(bs)
    lab = cb.SurfaceLabeller(V, T, tl)
    e0, c0 = block_edges(bs, 9, None)
    e1, c1 = block_edges(bs, 9, lab)
    assert (e0 == e1).all() and np.allclose(c0, c1, atol=1e-9)


def test_each_edge_keeps_its_two_endpoints():
    bs = _cube()
    V, T, tl = _boundary_triangles(bs)
    ends, curves = block_edges(bs, 7, cb.SurfaceLabeller(V, T, tl))
    assert curves.shape == (len(ends), 7, 3)
    for (a, b), c in zip(ends, curves):
        assert np.allclose(c[0], bs.points[a]) and np.allclose(c[-1], bs.points[b])


# --- B14: feature edges follow the seam curve between their two patches -------------------------------

def _seam():
    poly = np.c_[np.linspace(0, 1, 11), np.zeros(11), np.zeros(11)]
    return [(poly, False, (1, 2))]


def test_a_crease_edge_follows_its_seam():
    q = np.array([[0.2, 0.4, 0.0], [0.5, 0.0, 0.0], [0.8, -0.4, 0.0]])
    out = _seam_curve(_seam(), {1, 2}, q, 5)
    assert out is not None and out.shape == (5, 3)
    assert np.allclose(out[:, 2], 0.0) and np.allclose(out[1:-1, 1], 0.0)
    assert np.allclose(out[0], q[0]) and np.allclose(out[-1], q[-1])


def test_a_seam_of_another_patch_pair_is_not_used():
    q = np.array([[0.0, 0.1, 0.0], [1.0, -0.1, 0.0]])
    assert _seam_curve(_seam(), {1, 3}, q, 5) is None


def test_a_seam_far_from_the_corners_is_refused():
    q = np.array([[0.0, 2.0, 0.0], [1.0, 2.0, 0.0]])
    assert _seam_curve(_seam(), {1, 2}, q, 5) is None
