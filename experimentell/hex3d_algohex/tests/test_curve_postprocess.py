"""Acceptance checks for curve_postprocess.py on synthetic samples (no data files needed).

Run: PY -m pytest experimentell/hex3d_algohex/tests/test_curve_postprocess.py -q
or simply: PY experimentell/hex3d_algohex/tests/test_curve_postprocess.py
"""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

import curve_postprocess as cp                                                # noqa: E402


def _box_surface(lo, hi, n=8):
    """Triangulated box surface, one label per face (1..6)."""
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    pts, tris, labs = [], [], []
    lab = 0
    for ax in range(3):
        for side in (0, 1):
            lab += 1
            u, v = [d for d in range(3) if d != ax]
            g = np.linspace(0, 1, n + 1)
            base = len(pts)
            for i in g:
                for j in g:
                    p = np.empty(3)
                    p[ax] = hi[ax] if side else lo[ax]
                    p[u] = lo[u] + i * (hi[u] - lo[u]); p[v] = lo[v] + j * (hi[v] - lo[v])
                    pts.append(p)
            for i in range(n):
                for j in range(n):
                    a = base + i * (n + 1) + j
                    tris += [[a, a + n + 1, a + 1], [a + 1, a + n + 1, a + n + 2]]
                    labs += [lab, lab]
    pts = np.asarray(pts)
    # weld the shared box edges so the label boundaries are connected seams
    key = np.round(pts, 9)
    _, first, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    return pts[first], inv.ravel()[np.asarray(tris)], np.asarray(labs)


def _two_block_sample():
    """Two unit hexes along x, filling the box [0,2]x[0,1]x[0,1]; straight edges."""
    V = np.array([[x, y, z] for z in (0, 1) for y in (0, 1) for x in (0, 1, 2)], float)
    idx = {tuple(p): i for i, p in enumerate(map(tuple, V))}

    def hexa(x0):
        c = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
        return [idx[(x0 + a, b, cc)] for a, b, cc in c]
    blocks = np.array([hexa(0), hexa(1)])
    E = set()
    for b in blocks:
        for f in cp.FACES:
            for i in range(4):
                E.add(tuple(sorted((b[f[i]], b[f[(i + 1) % 4]]))))
    edges, ctrl = [], []
    for a, b in sorted(E):
        for p, q in ((a, b), (b, a)):
            edges.append([p, q])
            ctrl.append([(2 * V[p] + V[q]) / 3, (V[p] + 2 * V[q]) / 3])
    pts, tris, labs = _box_surface((0, 0, 0), (2, 1, 1))
    return {"vertices": V, "blocks": blocks, "edges": np.array(edges), "edge_ctrl": np.array(ctrl),
            "surface_points": pts, "surface_tris": tris, "surface_tri_label": labs}


def test_box_corners_feature_lines_and_pairs():
    s = _two_block_sample()
    out = cp.postprocess(s, verbose=False)
    import json
    q = json.loads(out["post_quality"])
    # 12 box edges, the 4 along x split by the block interface -> 16 feature edges, all on a seam
    assert q["feature_edges"] == 16 and q["feature_arcs"] == 16
    # the 8 box corners touch 3 patches; the 4 interface vertices on the x-edges are feature-line vertices
    role = out["post_vertex_role"]
    assert (role == cp.ROLE["corner"]).sum() == 8
    assert (role == cp.ROLE["feature"]).sum() == 4
    # one G1 pair per feature-line vertex, none at corners, no interior grid lines in a single layer
    assert q["pairs"] == {"feature": 4}
    assert all(role[v] == cp.ROLE["feature"] for v in out["post_pairs"][:, 0])
    # straight edges stay straight, nothing moves
    assert q["vertex_moved_max"] < 1e-9
    assert q["kink_after_max"]["feature"] < 1e-6
    assert len(out["post_edge_ctrl"]) == len(s["edges"]) == len(out["post_arc_offset"]) - 1


def test_g1_cubic_on_quarter_circle():
    t = np.linspace(0, np.pi / 2, 200)
    P = np.stack([np.cos(t), np.sin(t), 0 * t], 1)
    B1, B2 = cp.g1_cubic(P, np.array([0.0, 1, 0]), np.array([-1.0, 0, 0]))
    # tangents exact by construction; least squares in chord-length parameters gives arms ~0.545
    # (the max-error optimum is 0.5523) and a radial error below 0.5 %
    assert abs(B1[0] - 1) < 1e-12 and abs(B2[1] - 1) < 1e-12
    assert abs(B1[1] - 0.5523) < 0.02 and abs(B2[0] - 0.5523) < 0.02
    Q = cp.bezier(P[0], B1, B2, P[-1], 400)
    assert np.abs(np.linalg.norm(Q[:, :2], axis=1) - 1).max() < 5e-3


def test_seam_spline_keeps_closed_loop_smooth():
    # closed polygon around an ellipse whose tip radius (0.25) is about two node spacings: the spline
    # turn grows ~linearly with the window (smooth, no corner). A tip radius far below the node spacing
    # is unresolved by the triangulation and reads as a corner -- by design.
    t = np.linspace(0, 2 * np.pi, 41)[:-1]
    loop = np.stack([np.cos(t), 0.5 * np.sin(t), 0 * t], 1)
    seams = cp.Seams([(np.vstack([loop, loop[:1]]), True)])
    tip = np.array([1.0, 0.0, 0.0])
    h = seams.spacing
    a1, _ = seams.turn(tip, h)
    a2, _ = seams.turn(tip, 2 * h)
    assert a2 > 1.5 * a1 > 0


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn(); print("ok", name)
