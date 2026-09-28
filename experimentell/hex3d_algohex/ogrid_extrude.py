"""Stage 6b: a CONFORMING blade O-grid, extruded from the AlgoHex core.

`reattach.py` reuses the dtOO O-grid verbatim. Its five blocks do not match
the core's block faces on the cut surface (measured on the tistos template:
4 core corners on the interface, 1.0 ... 1.3 from the nearest O-grid seam at
a block size of ~1.5), so the assembly is non-conforming at block level --
T-junctions everywhere on the interface, and no shared nodes for the CFD.

Here the core decides the partition instead. Every core block face on
`ogrid_interface` is extruded inward to the blade wall, one block per face,
exactly like the hub/shroud boundary layer in `reattach.boundary_layer_blocks`.
The dtOO O-grid is kept only as a MAP: its structured wall-normal grid lines
run from the cut surface to the blade wall, so a point on the interface is
carried to the wall along the same line the O-grid would use. That also
carries the O-grid's own wall clustering, and the sharp trailing edge is
handled by the O-grid's parametrisation rather than by AlgoHex.

Conforming by construction: extruded columns are computed per interface
VERTEX, so two faces sharing an edge share every extruded node.
"""

import sys
from collections import Counter, defaultdict, deque
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parent.parent)]

import base_complex as bc                                             # noqa: E402
import clean_blocks as cb                                             # noqa: E402
import ovm_io                                                         # noqa: E402
import reattach as ra                                                 # noqa: E402
from dp3d.extraction import parse_msh                                 # noqa: E402

HEX_EDGES = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
             (0, 4), (1, 5), (2, 6), (3, 7)]


class OgridMap:
    """Wall-normal structure of the dtOO O-grid: layer index per node
    (0 = blade wall), a parent pointer one layer down, and the outermost
    layer's quads for locating interface points."""

    def __init__(self, P, H, verbose=True):
        self.P, self.H = P, H
        f2h, _ = bc.build_topology(H)
        bf = [tuple(int(v) for v in fk) for fk, hs in f2h.items() if len(hs) == 1]
        r = np.hypot(P[:, 0], P[:, 1])
        on_wall = lambda f, R: np.all(np.abs(r[list(f)] - R) < 1e-6)  # noqa: E731
        rh, rs = r.min(), r.max()
        tube = [f for f in bf if not (on_wall(f, rh) or on_wall(f, rs))]
        comps = _face_components(tube)
        if len(comps) != 2:
            raise RuntimeError(f"O-grid side boundary has {len(comps)} "
                               "components, expected inner + outer tube")
        # the blade is the tube with the smaller extent
        ext = [np.ptp(P[np.unique(np.ravel(c))], 0).sum() for c in comps]
        inner, outer = (comps[0], comps[1]) if ext[0] < ext[1] else (comps[1], comps[0])
        wall = np.unique(np.ravel(inner))

        adj = defaultdict(set)
        for h in H:
            for a, b in HEX_EDGES:
                adj[int(h[a])].add(int(h[b]))
                adj[int(h[b])].add(int(h[a]))
        lay = -np.ones(len(P), int)
        lay[wall] = 0
        q = deque(int(v) for v in wall)
        while q:
            v = q.popleft()
            for w in adj[v]:
                if lay[w] < 0:
                    lay[w] = lay[v] + 1
                    q.append(w)
        parent = -np.ones(len(P), int)
        for v in range(len(P)):
            if lay[v] > 0:
                down = [w for w in adj[v] if lay[w] == lay[v] - 1]
                if len(down) != 1:
                    raise RuntimeError(f"node {v} at layer {lay[v]} has "
                                       f"{len(down)} parents, not structured")
                parent[v] = down[0]
        self.lay, self.parent = lay, parent
        self.L = int(lay[np.unique(np.ravel(outer))].min())
        if lay[np.unique(np.ravel(outer))].max() != self.L:
            raise RuntimeError("outer tube is not a single O-grid layer")
        self.outer = np.array([_quad_loop(f, H) for f in outer], np.int64)
        self.tree = cKDTree(P[self.outer].mean(1))
        if verbose:
            print(f"[ogrid] {len(H)} hexes, {self.L} wall-normal layers, "
                  f"{len(wall)} blade-wall nodes, {len(self.outer)} outer quads")

    def ancestor(self, v, depth):
        for _ in range(depth):
            v = self.parent[v]
        return v

    def column(self, Q, k=12):
        """Wall-normal polyline [L+1, 3] through interface point Q, from the
        outer layer (index 0) to the blade wall (index L), plus the distance
        from Q to the outer layer."""
        _d, cand = self.tree.query(Q, k=k)
        best = None
        for qi in np.atleast_1d(cand):
            quad = self.outer[qi]
            s, t, res = _bilinear_inverse(self.P[quad], Q)
            out = max(0.0, -s, s - 1, -t, t - 1)
            key = (out > 1e-6, res)
            if best is None or key < best[0]:
                best = (key, quad, s, t, res)
        _key, quad, s, t, res = best
        s, t = np.clip(s, 0, 1), np.clip(t, 0, 1)
        C = np.empty((self.L + 1, 3))
        for depth in range(self.L + 1):
            a, b, c, d = (self.P[self.ancestor(int(v), depth)] for v in quad)
            C[depth] = ((1 - s) * (1 - t) * a + s * (1 - t) * b
                        + s * t * c + (1 - s) * t * d)
        return C, res


def _face_components(faces):
    par = list(range(len(faces)))

    def find(i):
        while par[i] != i:
            par[i] = par[par[i]]
            i = par[i]
        return i
    by_v = defaultdict(list)
    for i, f in enumerate(faces):
        for v in f:
            by_v[v].append(i)
    for fs in by_v.values():
        for j in fs[1:]:
            par[find(j)] = find(fs[0])
    out = defaultdict(list)
    for i, f in enumerate(faces):
        out[find(i)].append(f)
    return list(out.values())


def _quad_loop(face, H):
    """Order a face's 4 vertex ids as a loop, using any hex that owns it."""
    fs = set(face)
    for h in H:
        if fs <= set(int(x) for x in h):
            for loop in ((0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
                         (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)):
                if set(int(h[i]) for i in loop) == fs:
                    return [int(h[i]) for i in loop]
    raise RuntimeError("face not found")


def _bilinear_inverse(X, Q, iters=12):
    """(s, t) minimising |X(s,t) - Q| on a bilinear quad X[4,3]."""
    a, b, c, d = X
    s = t = 0.5
    for _ in range(iters):
        F = (1 - s) * (1 - t) * a + s * (1 - t) * b + s * t * c + (1 - s) * t * d
        Fs = (1 - t) * (b - a) + t * (c - d)
        Ft = (1 - s) * (d - a) + s * (c - b)
        J = np.column_stack([Fs, Ft])
        step, *_ = np.linalg.lstsq(J, Q - F, rcond=None)
        s, t = s + step[0], t + step[1]
    F = (1 - s) * (1 - t) * a + s * (1 - t) * b + s * t * c + (1 - s) * t * d
    return s, t, float(np.linalg.norm(F - Q))


def interface_faces(S):
    """[(block, quad loops)] for every block SIDE (lattice face) lying on
    `ogrid_interface`.

    Per lattice face, not per label patch: a block can sit on the interface
    with two of its sides (measured on the tistos template, a 55-quad patch
    with 5 corners), and extruding that patch as one piece gives a block that
    is not a lattice. Two sides extrude to two blocks sharing a radial face."""
    import tfi
    sid = {k for k, v in S.surf_names.items() if v == "ogrid_interface"}
    B = np.empty(len(S.hexes), int)
    for r, cells in S.cells_of().items():
        B[list(cells)] = r
    lat, _missing = tfi.lattices(S.P, S.hexes, B, S.f2h, verbose=False)
    out = []
    for r, (_dims, vert) in lat.items():
        for ax in (0, 1, 2):
            for side in (0, 1):
                G = tfi.face_grid(vert, ax, side)
                loops, on = [], 0
                for i in range(G.shape[0] - 1):
                    for j in range(G.shape[1] - 1):
                        q = [int(G[i, j]), int(G[i + 1, j]),
                             int(G[i + 1, j + 1]), int(G[i, j + 1])]
                        loops.append(q)
                        fk = frozenset(q)
                        if len(S.f2h.get(fk, ())) == 1 and S.surf_of.get(fk) in sid:
                            on += 1
                # Majority, not all: the nearest-face labels leave single
                # stray quads (tistos template: 1 of 6 on a neighbouring side
                # of block 14), which must neither add a side nor drop one.
                if 2 * on > len(loops):
                    out.append((r, loops))
                if on and on != len(loops):
                    print(f"[ogrid] block {r} side ({ax},{side}): {on}/"
                          f"{len(loops)} quads labelled ogrid_interface, "
                          f"{'taken' if 2 * on > len(loops) else 'ignored'}")
    return out


def extrude_ogrid(S, omap, n_layers=None, verbose=True):
    """One hex block per core block face on `ogrid_interface`, extruded to
    the blade wall along the O-grid's grid lines. `n_layers` cells in the
    wall-normal direction, sampled at the O-grid's own layer positions
    (default: all of them, i.e. the O-grid's wall clustering verbatim)."""
    faces = interface_faces(S)
    verts = sorted({int(v) for _r, loops in faces for lp in loops for v in lp})
    L = omap.L
    m = n_layers or L
    idx = np.round(np.linspace(0, L, m + 1)).astype(int)
    loc = {v: i for i, v in enumerate(verts)}
    cols, res = [], []
    for v in verts:
        C, rr = omap.column(S.P[v])
        # The core vertex sits up to ~0.06 off the O-grid's outer layer (two
        # discretisations of one surface). Fading that offset out towards the
        # wall was tried and measured worse (124 inverted against 74).
        C[0] = S.P[v]
        cols.append(C[idx])
        res.append(rr)
    cols = np.asarray(cols)          # [nv, m+1, 3]
    n = len(verts)
    # layer 0 are the core's own vertices; layers 1..m are new points
    newP = cols[:, 1:, :].transpose(1, 0, 2).reshape(-1, 3)
    base = len(S.P)

    def pid(v, layer):
        return int(v) if layer == 0 else base + (layer - 1) * n + loc[int(v)]
    H, B = [], []
    for bi, (_r, loops) in enumerate(faces):
        for lp in loops:
            for li in range(m):
                H.append([pid(v, li) for v in lp] + [pid(v, li + 1) for v in lp])
                B.append(bi)
    P = np.vstack([S.P, newP])
    H = np.asarray(H, np.int64)
    vol = np.array([ovm_io._hex_volume(P[c]) for c in H])
    H[vol < 0] = H[vol < 0][:, [4, 5, 6, 7, 0, 1, 2, 3]]
    if verbose:
        res = np.asarray(res)
        sj = cb.scaled_jacobians(P, H)
        print(f"[ogrid] {len(faces)} core faces on ogrid_interface, "
              f"{sum(len(l) for _r, l in faces)} quads, {n} vertices; "
              f"interface -> O-grid outer layer: median {np.median(res):.5f} "
              f"max {res.max():.5f}")
        print(f"[ogrid] extruded {len(H)} hexes in {len(faces)} blocks, "
              f"{m} layers; scaled Jacobian min {sj.min():.4f} mean "
              f"{sj.mean():.4f}, {int((sj <= 0).sum())} inverted")
    return P, H, np.asarray(B, int), faces


if __name__ == "__main__":
    import argparse
    import meshio
    import tet_prep_v5 as v5
    ap = argparse.ArgumentParser()
    ap.add_argument("blocks")
    ap.add_argument("--msh", required=True)
    ap.add_argument("--input-vtk", required=True)
    ap.add_argument("--layers", type=int, default=None)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    m = meshio.read(a.blocks)
    hx = np.vstack([c.data for c in m.cells if c.type == "hexahedron"])
    bid = np.concatenate([np.asarray(d).ravel() for c, d in
                          zip(m.cells, m.cell_data["block_id"])
                          if c.type == "hexahedron"])
    P0, tri0, tid0 = cb.read_input_surface(a.input_vtk)
    lab = cb.SurfaceLabeller(P0, tri0, tid0, v5.NAMES)
    S = cb.BlockStructure.from_arrays(m.points, hx, bid, lab)
    nodes, elements = parse_msh(a.msh)
    Po, Ho, _Bo = ra.ogrid_blocks(elements, nodes, verbose=False)
    omap = OgridMap(Po, Ho)
    P, Hx, Bx, _f = extrude_ogrid(S, omap, n_layers=a.layers)
    n_core = int(bid.max()) + 1
    Hall = np.vstack([S.hexes, Hx])
    Ball = np.concatenate([bid, Bx + n_core])
    part = np.concatenate([np.zeros(len(S.hexes), int), np.ones(len(Hx), int)])
    sj = cb.scaled_jacobians(P, Hall)
    print(f"[ogrid] core + O-grid: {len(Hall)} hexes, {int(Ball.max()) + 1} "
          f"blocks, {int((sj <= 0).sum())} inverted")
    import export_vtk as ev
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    ev.write_vtk(a.out, P, Hall, [12] * len(Hall), Ball, "block_id",
                 "core + extruded conforming O-grid")
    ev.write_vtk(a.out.replace(".vtk", "_part.vtk"), P, Hall, [12] * len(Hall),
                 part, "part", "0 = AlgoHex core, 1 = extruded O-grid")
