"""snap_to_surface: an inflated hex cube is projected back onto the input cube, edges and corners included."""
from types import SimpleNamespace

import numpy as np

from hexdp.core import base_complex as bc
from hexdp.core import clean_blocks as cb


def hex_cube(n=4):
    g = np.linspace(0, 1, n + 1)
    idx = np.arange((n + 1) ** 3).reshape(n + 1, n + 1, n + 1)
    P = np.array([(x, y, z) for x in g for y in g for z in g])
    H = [[idx[i, j, k], idx[i + 1, j, k], idx[i + 1, j + 1, k], idx[i, j + 1, k],
          idx[i, j, k + 1], idx[i + 1, j, k + 1], idx[i + 1, j + 1, k + 1], idx[i, j + 1, k + 1]]
         for i in range(n) for j in range(n) for k in range(n)]
    return P, np.array(H)


def input_cube():
    """The unit cube as 12 triangles, one label per side (1..6)."""
    V = np.array([(x, y, z) for x in (0, 1) for y in (0, 1) for z in (0, 1)], float)
    quads = {1: (0, 1, 3, 2), 2: (4, 6, 7, 5), 3: (0, 4, 5, 1), 4: (2, 3, 7, 6), 5: (0, 2, 6, 4), 6: (1, 5, 7, 3)}
    T, ids = [], []
    for lab, (a, b, c, d) in quads.items():
        T += [(a, b, c), (a, c, d)]
        ids += [lab, lab]
    return V, np.array(T), np.array(ids)


def test_snap_inflated_cube():
    P, H = hex_cube()
    f2h, _e2h = bc.build_topology(H)
    bnd = [(fk, hs[0]) for fk, hs in f2h.items() if len(hs) == 1]
    loops = []
    for fk, h in bnd:
        for fc in cb.HF:
            lp = [int(H[h][i]) for i in fc]
            if frozenset(lp) == fk:
                loops.append(lp)
                break
    V, T, ids = input_cube()
    lab = cb.SurfaceLabeller(V, T, ids)
    surf = lab.label(P[np.asarray(loops)].mean(1))
    on_b = np.unique(np.asarray(loops))
    Q = P.copy()
    Q[on_b] = 0.5 + 1.05 * (P[on_b] - 0.5)                  # inflate the boundary by 5 %
    S = SimpleNamespace(P=Q, hexes=H, f2h=f2h, bnd_loops=loops,
                        surf_of={fk: int(s) for (fk, _h), s in zip(bnd, surf)})
    r = cb.snap_to_surface(S, lab, verbose=False)
    assert r["kept"] and r["before"] > 0.02 and r["after"] < 1e-9
    # every boundary vertex on its side(s): as many coordinates at 0 or 1 as the original vertex had
    on_side = lambda X: (np.isclose(X, 0) | np.isclose(X, 1)).sum(1)
    assert (on_side(S.P[on_b]) == on_side(P[on_b])).all()
    corners = on_b[on_side(P[on_b]) == 3]
    assert np.allclose(S.P[corners], P[corners])            # the 8 corners exactly at the cube's corners
