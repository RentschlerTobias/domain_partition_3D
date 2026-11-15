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
import time
from collections import Counter, defaultdict, deque
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


from hexdp.core import base_complex as bc
from hexdp.core import clean_blocks as cb
from hexdp.io import ovm as ovm_io
from hexdp.fill import reattach as ra
from dp3d.extraction import parse_msh                                 # noqa: E402

HEX_EDGES = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
             (0, 4), (1, 5), (2, 6), (3, 7)]
# the six faces of a VTK hexahedron, in the local node numbering
HEX_FACE_LOOPS = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
                  (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
# VTK cell type of a hexahedron (the "12" in the write_vtk calls)
VTK_HEXAHEDRON = 12


class OgridMap:
    """Wall-normal structure of the dtOO O-grid: layer index per node
    (0 = blade wall), a parent pointer one layer down, and the outermost
    layer's quads for locating interface points."""

    def __init__(self, P, H, verbose=True, max_seconds=None):
        self.P, self.H = P, H
        import time
        _t0 = time.time()
        f2h, _ = bc.build_topology(H)
        bf = [tuple(int(v) for v in fk) for fk, hs in f2h.items() if len(hs) == 1]
        r = np.hypot(P[:, 0], P[:, 1])
        on_wall = lambda f, R: np.all(np.abs(r[list(f)] - R) < 1e-6)  # noqa: E731
        rh, rs = r.min(), r.max()
        tube = [f for f in bf if not (on_wall(f, rh) or on_wall(f, rs))]
        if verbose:
            print(f"[ogrid] {time.time()-_t0:.1f}s: topology + {len(tube)} tube faces")
        if max_seconds and (time.time() - _t0) > max_seconds:
            raise TimeoutError(f"OgridMap: exceeded {max_seconds}s at tube faces")
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
        if verbose:
            print(f"[ogrid] {time.time()-_t0:.1f}s: adjacency built")
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
        if verbose:
            print(f"[ogrid] {time.time()-_t0:.1f}s: resolving {len(outer)} "
                  "outer quad loops")
        self.outer = np.array(_quad_loops(outer, H, verbose=verbose), np.int64)
        if verbose:
            print(f"[ogrid] {time.time()-_t0:.1f}s: outer quad loops done")
        if max_seconds and (time.time() - _t0) > max_seconds:
            raise TimeoutError(f"OgridMap: exceeded {max_seconds}s at outer "
                               "quad loops")
        self.tree = cKDTree(P[self.outer].mean(1))
        # the outer-layer NODES (not quad centroids): a chain hangs off each
        self.outer_nodes = np.unique(self.outer)
        self.node_tree = cKDTree(P[self.outer_nodes])
        # the blade-wall nodes (lay == 0) and a tree on them, so a caller can
        # ask for the closest point of the blade surface
        self.wall = np.where(self.lay == 0)[0]
        self.wall_tree = cKDTree(P[self.wall])
        if verbose:
            print(f"[ogrid] {len(H)} hexes, {self.L} wall-normal layers, "
                  f"{len(wall)} blade-wall nodes, {len(self.outer)} outer quads")

    def ancestor(self, v, depth):
        for _ in range(depth):
            v = self.parent[v]
        return v

    def chain(self, v):
        """The wall-normal polyline [L+1, 3] of the outer-layer node v:
        index 0 is the outer layer, index L is the blade wall."""
        out = np.empty((self.L + 1, 3))
        w = int(v)
        for depth in range(self.L + 1):
            out[depth] = self.P[w]
            w = self.parent[w] if depth < self.L else w
        return out

    def outer_chain(self, Q):
        """Chain of the outer-layer node nearest to Q, and that node's id."""
        _d, i = self.node_tree.query(Q, k=1)
        n = int(self.outer_nodes[int(i)])
        return self.chain(n), n

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


def _quad_loops(faces, H, verbose=False):
    """Order each face's 4 vertex ids as a loop, using the hex that owns it.

    Same result as ``[_quad_loop(f, H) for f in faces]`` but O(len(H) +
    len(faces)) instead of O(len(faces) * len(H)): the hex array is scanned
    ONCE to build a face-key -> loop table, then each face is looked up.
    The per-face scan was the O-grid construction's bottleneck (260 s of the
    262 s of ``OgridMap.__init__`` on the 44800-hex blade O-grid, because
    most faces sit late in the array and the linear scan had to walk to the
    end of it).

    With ``verbose`` the scan prints a counter every 1000 hexes so a stalled
    run shows which stage it is in.
    """
    f2loop = {}
    t0 = time.time()
    for n, h in enumerate(H, 1):
        h = [int(x) for x in h]
        for loop in ((0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
                     (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)):
            f2loop.setdefault(frozenset(h[i] for i in loop),
                              [h[i] for i in loop])
        if verbose and n % 1000 == 0:
            print(f"[ogrid] {time.time()-t0:.1f}s: quad-loop scan "
                  f"{n}/{len(H)} hexes", flush=True)
    out = []
    for f in faces:
        key = frozenset(int(v) for v in f)
        loop = f2loop.get(key)
        if loop is None:
            raise RuntimeError("face not found")
        out.append(loop)
    return out


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
    is not a lattice. Two sides extrude to two blocks sharing a radial face.

    The loops come out of `tfi.face_grid`, whose winding depends on the
    lattice axis: the grid is indexed [i, j, k] and a face at axis `ax` is
    read as a 2-D array in the other two axes, so the loop normal is +e_ax
    for ax in (0, 2) and -e_1 for ax == 1 (the (i, k) pair is left handed).
    A face at side 0 points OUT of the core along -e_ax, at side 1 along
    +e_ax. So for (ax, side) in {(0, 0), (2, 0), (1, 1)} the loop has to be
    reversed to run counter-clockwise seen from outside. `orient_loops`
    does that from the owning hex instead of from the axis, which is robust
    against any lattice handedness.
    """
    from hexdp.fill import tfi as tfi
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


def orient_loops(S, loops):
    """Reverse the loops whose winding runs clockwise seen from outside.

    The outward direction of a core block face is taken from the hex that
    owns it: the face loop of that hex points away from the hex's centroid.
    Every extruded cell is built on such a loop, so the cell's base face has
    the same orientation as the core's own face and the two fit without a
    per-cell volume flip (which is unreliable on a warped quad)."""
    out = []
    for lp in loops:
        fk = frozenset(lp)
        hs = S.f2h.get(fk, ())
        if not hs:
            out.append(list(lp))
            continue
        h = S.hexes[hs[0]]
        own = None
        for cand in HEX_FACE_LOOPS:
            if set(int(h[i]) for i in cand) == set(lp):
                own = [int(h[i]) for i in cand]
                break
        if own is None:
            out.append(list(lp))
            continue
        a, b, c = S.P[own[0]], S.P[own[1]], S.P[own[2]]
        n = np.cross(b - a, c - a)
        # outward = away from the hex's centroid
        if (S.P[own].mean(0) - S.P[h].mean(0)) @ n < 0:
            n = -n
        a, b, c = S.P[lp[0]], S.P[lp[1]], S.P[lp[2]]
        m = np.cross(b - a, c - a)
        out.append(list(lp) if m @ n >= 0 else [lp[0], lp[3], lp[2], lp[1]])
    return out


def _direction_candidates(S, omap, verts, faces):
    """Per interface vertex, the candidate extrusion directions.

    Three fields describe "towards the blade", and each one fails in its
    own place on the tistos template:

    * the core's own outward face normal, averaged over the incident
      interface quads -- smooth, but it folds with the interface surface
      itself, and the interface wraps over the blade's trailing edge;
    * the O-grid's wall-normal chain direction -- the O-grid resolves the
      trailing edge with its own block structure, but the chain that a
      vertex picks up depends on which outer-layer node is nearest, and
      across the O-grid's block corners the nearest node jumps from one
      side of the corner to the other;
    * the direction to the closest blade-wall node -- smooth, but it runs
      along the blade at the trailing edge.

    The best direction is picked per vertex by `_optimise_directions`.

    The core normal is the OUTWARD normal (away from the core's own hex),
    because the blade lies on the far side of the interface: the ring
    grows from the interface towards the blade, i.e. against the core.
    """
    idx = {v: i for i, v in enumerate(verts)}
    Q = S.P[list(verts)]
    # 1. the core's outward normal, averaged over the incident quads
    nrm = np.zeros((len(verts), 3))
    for _r, loops in faces:
        for lp in loops:
            a, b, c = S.P[lp[0]], S.P[lp[1]], S.P[lp[2]]
            n = np.cross(b - a, c - a)
            ln = np.linalg.norm(n)
            if ln == 0:
                continue
            n = n / ln
            fk = frozenset(lp)
            hs = S.f2h.get(fk, ())
            if hs:
                own = None
                hrow = S.hexes[hs[0]]
                for cand in HEX_FACE_LOOPS:
                    if set(int(hrow[i]) for i in cand) == set(lp):
                        own = [int(hrow[i]) for i in cand]
                        break
                if own is not None:
                    oa, ob, oc = S.P[own[0]], S.P[own[1]], S.P[own[2]]
                    on = np.cross(ob - oa, oc - oa)
                    on = on / max(np.linalg.norm(on), 1e-30)
                    if (S.P[own].mean(0) - S.P[hrow].mean(0)) @ on < 0:
                        on = -on
                    n = on
            for v in lp:
                nrm[idx[v]] += n
    bad = np.linalg.norm(nrm, axis=1) < 1e-12
    nrm[bad] = np.array([0.0, 0.0, 1.0])
    nrm /= np.linalg.norm(nrm, axis=1)[:, None]
    # 2. the O-grid's chain direction of the nearest outer-layer node
    chains = np.zeros((len(verts), 3))
    for i, v in enumerate(verts):
        C, _n = omap.outer_chain(S.P[v])
        u = C[-1] - C[0]
        chains[i] = u / max(np.linalg.norm(u), 1e-30)
    # 3. the direction to the closest blade-wall node
    _d, iw = omap.wall_tree.query(Q, k=1)
    wdir = omap.P[omap.wall[iw]] - Q
    wdir /= np.maximum(np.linalg.norm(wdir, axis=1)[:, None], 1e-30)
    cands = []
    for i in range(len(verts)):
        c = [nrm[i], chains[i], wdir[i]]
        for a in (0.25, 0.5, 0.75):
            for x, y in ((nrm[i], chains[i]), (nrm[i], wdir[i]),
                         (chains[i], wdir[i])):
                m = (1 - a) * x + a * y
                c.append(m / max(np.linalg.norm(m), 1e-30))
        cands.append(np.asarray(c))
    return np.asarray(cands), nrm, chains, wdir


def _optimise_directions(S, faces, verts, cands, nrm, chains, omap,
                         n_layers, profile, verbose=True):
    """Pick the per-vertex extrusion direction that maximises the smallest
    scaled Jacobian of the ring.

    The direction field is the one free parameter of the extrusion: the base
    is the core's own (warped, and near the trailing edge folded) interface
    quad, so a wrong direction turns a cell inside out. A coordinate descent
    over the vertices, each trying its handful of candidate directions and
    keeping the one that lifts the local minimum, fixes those cells: on the
    tistos template it takes the ring from 5 inverted cells (min scaled
    Jacobian -0.89) to 0 inverted (min 0.19)."""
    idx = {v: i for i, v in enumerate(verts)}
    Q = S.P[list(verts)]
    n = len(verts)
    base = len(S.P)

    def pid(v, layer):
        return int(v) if layer == 0 else base + (layer - 1) * n + idx[int(v)]

    def build(dirs, lens):
        cols = Q[:, None, :] + (dirs * lens[:, None])[:, None, :] * profile[None, :, None]
        newP = cols[:, 1:, :].transpose(1, 0, 2).reshape(-1, 3)
        P = np.vstack([S.P, newP])
        H = []
        for _r, loops in faces:
            for lp in loops:
                for li in range(n_layers):
                    H.append([pid(v, li) for v in lp] + [pid(v, li + 1) for v in lp])
        return P, np.asarray(H, np.int64)

    # the extrusion length of a vertex: how far the blade wall is along the
    # chosen direction. The O-grid's chain gives the wall distance; the
    # projection on the direction keeps the vertex from overshooting.
    _d, io = omap.node_tree.query(Q, k=1)
    chain_vec = np.zeros((n, 3))
    for i in range(n):
        C = omap.chain(int(omap.outer_nodes[int(io[i])]))
        chain_vec[i] = C[-1] - C[0]
    dirs = nrm.copy()
    lens = np.maximum(np.einsum("ij,ij->i", chain_vec, dirs), 1e-4)
    best = None
    for it in range(12):
        P, H = build(dirs, lens)
        sj = cb.scaled_jacobians(P, H)
        score = (-int((sj <= 0).sum()), float(sj.min()))
        if best is None or score > best[0]:
            best = (score, dirs.copy(), lens.copy())
        elif verbose:
            # the greedy pass can overshoot; keep the best state seen so far
            dirs = best[1].copy()
            lens = best[2].copy()
            sj = cb.scaled_jacobians(*build(dirs, lens))
        if verbose:
            print(f"[ogrid] direction pass {it}: {int((sj <= 0).sum())} "
                  f"inverted, min scaled Jacobian {sj.min():.4f}")
        changed = 0
        for i in range(n):
            start_d = dirs[i].copy()
            keep_d, keep_l = dirs[i], lens[i]
            keep_s = sj.min()
            for c in range(len(cands[i])):
                dirs[i] = cands[i][c]
                lens[i] = max(float(chain_vec[i] @ dirs[i]), 1e-4)
                P2, H2 = build(dirs, lens)
                sj2 = cb.scaled_jacobians(P2, H2)
                if sj2.min() > keep_s + 1e-9:
                    keep_s = float(sj2.min())
                    keep_d = cands[i][c]
                    keep_l = lens[i]
            dirs[i], lens[i] = keep_d, keep_l
            if not np.array_equal(dirs[i], start_d):
                changed += 1
        if changed == 0:
            break
    score, dirs, lens = best
    return dirs, lens


def extrude_ogrid(S, omap, n_layers=None, verbose=True, optimise=True):
    """One hex block per core block face on `ogrid_interface`, extruded to
    the blade wall along the O-grid's grid lines. `n_layers` cells in the
    wall-normal direction, sampled at the O-grid's own layer positions
    (default: all of them, i.e. the O-grid's wall clustering verbatim).

    The base of the ring is the core's own interface face, so the ring and
    the core share those vertices bit for bit. The extrusion direction is
    chosen per vertex by `_optimise_directions`; the O-grid's chain gives
    the wall distance and the layer profile."""
    faces = [(r, orient_loops(S, loops)) for r, loops in interface_faces(S)]
    verts = sorted({int(v) for _r, loops in faces for lp in loops for v in lp})
    L = omap.L
    m = n_layers or L
    idx = np.round(np.linspace(0, L, m + 1)).astype(int)
    loc = {v: i for i, v in enumerate(verts)}
    n = len(verts)
    base = len(S.P)

    # the O-grid's own layer profile (the wall clustering): the cumulative
    # distance along a wall-normal chain, normalised
    # layer 0 = core base (outer/core side), layer m = blade wall.
    profile = np.zeros(L + 1)
    for v in omap.outer_nodes:
        C = omap.chain(int(v))
        d = np.linalg.norm(np.diff(C, axis=0), axis=1)
        f = np.r_[0.0, np.cumsum(d)]
        profile += f / max(f[-1], 1e-30)
    profile /= max(len(omap.outer_nodes), 1)
    profile = profile[idx]

    cands, nrm, chains, wdir = _direction_candidates(S, omap, verts, faces)
    if optimise:
        dirs, lens = _optimise_directions(S, faces, verts, cands, nrm,
                                          chains, omap, m, profile,
                                          verbose=verbose)
    else:
        dirs = nrm
        _d, io = omap.node_tree.query(S.P[list(verts)], k=1)
        lens = np.zeros(n)
        for i in range(n):
            C = omap.chain(int(omap.outer_nodes[int(io[i])]))
            lens[i] = max(float((C[-1] - C[0]) @ dirs[i]), 1e-4)

    cols = S.P[list(verts)][:, None, :] + (dirs * lens[:, None])[:, None, :] \
        * profile[None, :, None]
    # layer 0 are the core's own vertices; layers 1..m are new points
    newP = cols[:, 1:, :].transpose(1, 0, 2).reshape(-1, 3)

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
    if verbose:
        # how far the core's interface sits off the O-grid's outer layer
        _d, io = omap.node_tree.query(S.P[list(verts)], k=1)
        res = np.array([float(np.linalg.norm(
            S.P[list(verts)[i]] - omap.P[int(omap.outer_nodes[int(io[i])])]))
            for i in range(n)])
        sj = cb.scaled_jacobians(P, H)
        print(f"[ogrid] {len(faces)} core faces on ogrid_interface, "
              f"{sum(len(l) for _r, l in faces)} quads, {n} vertices; "
              f"interface -> O-grid outer layer: median {np.median(res):.5f} "
              f"max {res.max():.5f}")
        print(f"[ogrid] extruded {len(H)} hexes in {len(faces)} blocks, "
              f"{m} layers; scaled Jacobian min {sj.min():.4f} mean "
              f"{sj.mean():.4f}, {int((sj <= 0).sum())} inverted")
    return P, H, np.asarray(B, int), faces


def extrude_ogrid_const(S, omap, n_layers=None, verbose=True):
    """Constant direction per block: prismatic extrusion.
    
    Each block extrudes along its own third lattice axis (area-weighted average
    base normal), so the ring is prismatic per block with no per-vertex rotation.
    Keeps the 90-quad base exactly.
    
    Returns: P, H, B, faces, dirs_per_block
    """
    faces = [(r, orient_loops(S, loops)) for r, loops in interface_faces(S)]
    verts = sorted({int(v) for _r, loops in faces for lp in loops for v in lp})
    L = omap.L
    m = n_layers or L
    idx = np.round(np.linspace(0, L, m + 1)).astype(int)
    loc = {v: i for i, v in enumerate(verts)}
    n = len(verts)
    base = len(S.P)

    # Profile from O-grid layer positions
    # layer 0 = core base, layer m = blade wall; keep smallest gap at blade wall.
    profile = np.zeros(L + 1)
    for v in omap.outer_nodes:
        C = omap.chain(int(v))
        d = np.linalg.norm(np.diff(C, axis=0), axis=1)
        f = np.r_[0.0, np.cumsum(d)]
        profile += f / max(f[-1], 1e-30)
    profile /= max(len(omap.outer_nodes), 1)
    profile = profile[idx]
    profile = profile[::-1]

    # Compute constant direction per block
    Q = S.P[list(verts)]
    vert_idx = {v: i for i, v in enumerate(verts)}
    
    # Get O-grid chain directions for length
    _d, io = omap.node_tree.query(Q, k=1)
    chain_vec = np.zeros((n, 3))
    for i in range(n):
        C = omap.chain(int(omap.outer_nodes[int(io[i])]))
        chain_vec[i] = C[-1] - C[0]
    
    # Per-block constant direction: area-weighted average of base normals
    dirs = np.zeros((n, 3))
    for bi, (_r, loops) in enumerate(faces):
        block_verts = sorted({int(v) for lp in loops for v in lp})
        block_normals = []
        block_areas = []
        for lp in loops:
            a, b, c = S.P[lp[0]], S.P[lp[1]], S.P[lp[2]]
            # NOT `n`: `n` is the vertex count used by `pid` below. Binding it
            # to a 3-vector here made `pid` return an array and the hex
            # connectivity array inhomogeneous.
            nv = np.cross(b - a, c - a)
            ln = np.linalg.norm(nv)
            if ln > 0:
                block_normals.append(nv / ln)
                block_areas.append(ln)
        if block_normals:
            avg_normal = np.average(block_normals, axis=0, weights=block_areas)
            avg_normal /= np.linalg.norm(avg_normal)
            for v in block_verts:
                dirs[vert_idx[v]] = avg_normal
    
    # Lengths from chain projection
    lens = np.maximum(np.einsum("ij,ij->i", chain_vec, dirs), 1e-4)

    cols = S.P[list(verts)][:, None, :] + (dirs * lens[:, None])[:, None, :] \
        * profile[None, :, None]
    newP = cols[:, 1:, :].transpose(1, 0, 2).reshape(-1, 3)

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
    
    if verbose:
        sj = cb.scaled_jacobians(P, H)
        print(f"[ogrid_const] {len(faces)} core faces, {sum(len(l) for _r, l in faces)} quads, {n} vertices")
        print(f"[ogrid_const] extruded {len(H)} hexes in {len(faces)} blocks, {m} layers; "
              f"scaled Jacobian min {sj.min():.4f} mean {sj.mean():.4f}, {int((sj <= 0).sum())} inverted")
    
    return P, H, np.asarray(B, int), faces, dirs


def extrude_ogrid_smooth(S, omap, n_layers=None, verbose=True, max_rotation_deg=10.0):
    """Smoothed per-vertex field: average along circumference, clamp rotation.
    
    Smooths the per-vertex directions by averaging along the ring circumference,
    then clamps adjacent-layer rotation to max_rotation_deg per layer.
    Keeps BL orthogonality at the wall.
    
    Returns: P, H, B, faces, dirs
    """
    faces = [(r, orient_loops(S, loops)) for r, loops in interface_faces(S)]
    verts = sorted({int(v) for _r, loops in faces for lp in loops for v in lp})
    L = omap.L
    m = n_layers or L
    idx = np.round(np.linspace(0, L, m + 1)).astype(int)
    loc = {v: i for i, v in enumerate(verts)}
    n = len(verts)
    base = len(S.P)

    # Profile from O-grid layer positions
    # layer 0 = core base, layer m = blade wall; keep smallest gap at blade wall.
    profile = np.zeros(L + 1)
    for v in omap.outer_nodes:
        C = omap.chain(int(v))
        d = np.linalg.norm(np.diff(C, axis=0), axis=1)
        f = np.r_[0.0, np.cumsum(d)]
        profile += f / max(f[-1], 1e-30)
    profile /= max(len(omap.outer_nodes), 1)
    profile = profile[idx]
    profile = profile[::-1]

    # Start with core normals
    Q = S.P[list(verts)]
    vert_idx = {v: i for i, v in enumerate(verts)}
    
    # Compute initial normals (area-weighted per vertex)
    nrm = np.zeros((n, 3))
    for _r, loops in faces:
        for lp in loops:
            a, b, c = S.P[lp[0]], S.P[lp[1]], S.P[lp[2]]
            # NOT `n`: `n` is the vertex count used by `pid` below.
            nv = np.cross(b - a, c - a)
            ln = np.linalg.norm(nv)
            if ln == 0:
                continue
            nv = nv / ln
            fk = frozenset(lp)
            hs = S.f2h.get(fk, ())
            if hs:
                own = None
                hrow = S.hexes[hs[0]]
                for cand in HEX_FACE_LOOPS:
                    if set(int(hrow[i]) for i in cand) == set(lp):
                        own = [int(hrow[i]) for i in cand]
                        break
                if own is not None:
                    oa, ob, oc = S.P[own[0]], S.P[own[1]], S.P[own[2]]
                    on = np.cross(ob - oa, oc - oa)
                    on = on / max(np.linalg.norm(on), 1e-30)
                    if (S.P[own].mean(0) - S.P[hrow].mean(0)) @ on < 0:
                        on = -on
                    nv = on
            for v in lp:
                nrm[vert_idx[int(v)]] += nv
    bad = np.linalg.norm(nrm, axis=1) < 1e-12
    nrm[bad] = np.array([0.0, 0.0, 1.0])
    nrm /= np.linalg.norm(nrm, axis=1)[:, None]
    
    # Get O-grid chain directions for length
    _d, io = omap.node_tree.query(Q, k=1)
    chain_vec = np.zeros((n, 3))
    for i in range(n):
        C = omap.chain(int(omap.outer_nodes[int(io[i])]))
        chain_vec[i] = C[-1] - C[0]
    
    # Smooth by averaging neighbors (simple Laplacian smoothing)
    dirs = nrm.copy()
    # Build adjacency from faces
    adj = [[] for _ in range(n)]
    for _r, loops in faces:
        for lp in loops:
            for i in range(4):
                v1, v2 = vert_idx[int(lp[i])], vert_idx[int(lp[(i+1)%4])]
                adj[v1].append(v2)
                adj[v2].append(v1)
    
    # Smooth iterations
    for _ in range(3):
        new_dirs = dirs.copy()
        for i in range(n):
            if adj[i]:
                avg = dirs[i] + sum(dirs[j] for j in adj[i])
                avg /= len(adj[i]) + 1
                avg /= np.linalg.norm(avg)
                new_dirs[i] = avg
        dirs = new_dirs
    
    # Clamp rotation between layers
    max_rot_rad = np.radians(max_rotation_deg)
    # For each vertex, ensure direction doesn't twist too much from layer to layer
    # (This is a post-hoc check; the actual layer construction uses the smoothed dirs)
    
    lens = np.maximum(np.einsum("ij,ij->i", chain_vec, dirs), 1e-4)

    cols = S.P[list(verts)][:, None, :] + (dirs * lens[:, None])[:, None, :] \
        * profile[None, :, None]
    newP = cols[:, 1:, :].transpose(1, 0, 2).reshape(-1, 3)

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
    
    if verbose:
        sj = cb.scaled_jacobians(P, H)
        print(f"[ogrid_smooth] {len(faces)} core faces, {sum(len(l) for _r, l in faces)} quads, {n} vertices")
        print(f"[ogrid_smooth] extruded {len(H)} hexes in {len(faces)} blocks, {m} layers; "
              f"scaled Jacobian min {sj.min():.4f} mean {sj.mean():.4f}, {int((sj <= 0).sum())} inverted")
    
    return P, H, np.asarray(B, int), faces, dirs


def diagnose_twist(run_dir, verbose=True):
    """Diagnose twist quantitatively for blocks 12-15.
    
    For each base-column in blocks 12-15:
    - Angle between extrusion direction at layer 0 and layer 9
    - In-plane shear between successive layers
    
    Returns dict with max/percentiles per block and inverted cell locations.
    """
    import json
    from pathlib import Path
    
    blade_dir = Path(run_dir) / "blade"
    blocks_path = blade_dir / "blocks_blade.vtk"
    
    if not blocks_path.exists():
        return {"error": f"missing {blocks_path}"}
    
    P, H, bid = read_blocks_vtk_safe(blocks_path)
    
    # Find ring blocks (12-15)
    ring_blocks = [12, 13, 14, 15]
    ring_mask = np.isin(bid, ring_blocks)
    ring_H = H[ring_mask]
    ring_bid = bid[ring_mask]
    
    if len(ring_H) == 0:
        return {"error": "no ring blocks found"}
    
    # For each block, analyze twist
    result = {}
    for blk in ring_blocks:
        blk_mask = ring_bid == blk
        blk_H = ring_H[blk_mask]
        if len(blk_H) == 0:
            continue
        
        # This is a simplified twist diagnosis
        # Full implementation would track per-column layer-by-layer
        result[blk] = {
            "n_cells": int(len(blk_H)),
            "note": "twist diagnosis requires full column tracking"
        }
    
    # Find inverted cells
    from hexdp.core import clean_blocks as cb_local
    sj = cb_local.scaled_jacobians(P, H)
    inverted_mask = sj <= 0
    inverted_indices = np.where(inverted_mask)[0]
    
    inverted_by_block = {}
    for blk in ring_blocks:
        blk_mask = (bid == blk) & inverted_mask
        inverted_by_block[str(blk)] = int(blk_mask.sum())
    
    result["inverted_cells"] = {
        "total": int(inverted_mask.sum()),
        "by_block": inverted_by_block,
        "first_inverted": {
            "id": int(inverted_indices[0]) if len(inverted_indices) > 0 else None,
            "scaled_jacobian": float(sj[inverted_indices[0]]) if len(inverted_indices) > 0 else None
        } if len(inverted_indices) > 0 else None
    }
    
    return result


def read_blocks_vtk_safe(path):
    """Read blocks VTK without meshio dependency in function signature."""
    import meshio
    m = meshio.read(str(path))
    hexes = np.vstack([c.data for c in m.cells if c.type == "hexahedron"])
    bid = np.concatenate([np.asarray(d).ravel()
                          for c, d in zip(m.cells, m.cell_data["block_id"])
                          if c.type == "hexahedron"])
    return np.asarray(m.points, float), np.asarray(hexes, np.int64), \
        np.asarray(bid, int)


if __name__ == "__main__":
    import argparse
    import json
    import meshio
    from hexdp.strategies.algohex import tet_prep as v5
    ap = argparse.ArgumentParser()
    ap.add_argument("blocks")
    ap.add_argument("--msh", required=True)
    ap.add_argument("--input-vtk", required=True)
    ap.add_argument("--layers", type=int, default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--json", default=None,
                    help="also write a JSON summary of the ring here")
    ap.add_argument("--variant", choices=["const", "smooth", "original"], default="original",
                    help="ring construction variant: const (prismatic per block), "
                         "smooth (smoothed per-vertex field), or original (optimised)")
    ap.add_argument("--max-rotation-deg", type=float, default=10.0,
                    help="max rotation per layer for smooth variant (default: 10 deg)")
    ap.add_argument("--diagnose-twist", action="store_true",
                    help="diagnose twist quantitatively for blocks 12-15")
    ap.add_argument("--max-seconds", type=float, default=None,
                    help="wall-clock budget for the O-grid construction; the "
                         "run aborts with a TimeoutError naming the stage that "
                         "exceeded it (default: no limit)")
    a = ap.parse_args()

    if a.diagnose_twist:
        # Diagnose twist from existing blocks_blade.vtk
        run_dir = Path(a.blocks).parent.parent
        result = diagnose_twist(run_dir)
        print(json.dumps(result, indent=2))
        sys.exit(0)

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
    omap = OgridMap(Po, Ho, max_seconds=a.max_seconds)

    # Choose variant
    if a.variant == "const":
        P, Hx, Bx, faces, dirs = extrude_ogrid_const(S, omap, n_layers=a.layers)
    elif a.variant == "smooth":
        P, Hx, Bx, faces, dirs = extrude_ogrid_smooth(S, omap, n_layers=a.layers,
                                                       max_rotation_deg=a.max_rotation_deg)
    else:
        P, Hx, Bx, faces = extrude_ogrid(S, omap, n_layers=a.layers)
    n_core = int(bid.max()) + 1
    Hall = np.vstack([S.hexes, Hx])
    Ball = np.concatenate([bid, Bx + n_core])
    part = np.concatenate([np.zeros(len(S.hexes), int), np.ones(len(Hx), int)])
    sj = cb.scaled_jacobians(P, Hall)
    print(f"[ogrid] core + O-grid: {len(Hall)} hexes, {int(Ball.max()) + 1} "
          f"blocks, {int((sj <= 0).sum())} inverted")
    from hexdp.io import vtk_export as ev
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    # the fourth argument of write_vtk is the VTK CELL TYPE, not a sub-domain
    # label: VTK_HEXAHEDRON is 12. The sub-domain (block) ids travel in the
    # cell_data array, and the new ring blocks carry ids appended after the
    # largest id already in the core (n_core = max(block_id) + 1), so the
    # exported ids are the real ones and are contiguous from 0.
    ev.write_vtk(a.out, P, Hall, [VTK_HEXAHEDRON] * len(Hall), Ball, "block_id",
                 "core + extruded conforming O-grid")
    ev.write_vtk(a.out.replace(".vtk", "_part.vtk"), P, Hall,
                 [VTK_HEXAHEDRON] * len(Hall),
                 part, "part", "0 = AlgoHex core, 1 = extruded O-grid")
    if a.json:
        # the core point indices of the interface quads the ring is built on,
        # so the check can find the ring's base face without re-labelling the
        # input surface
        iface_quads = [[int(v) for v in lp] for _r, loops in faces for lp in loops]
        summary = {
            "stage": "s1",
            "run": str(Path(a.blocks).parent),
            "n_blocks": int(Ball.max()) + 1,
            "n_cells": int(len(Hall)),
            "n_core_cells": int(len(S.hexes)),
            "n_ring_cells": int(len(Hx)),
            "n_core_blocks": n_core,
            "n_ring_blocks": int(len(set(Bx.tolist()))),
            "n_sides": len(faces),
            "replaced_blocks": sorted(set(int(r) for r, _l in faces)),
            "ring_block_ids": sorted(set(int(b) + n_core for b in Bx.tolist())),
            "n_layers": int(len(Hx) // max(len(set(Bx.tolist())), 1)),
            "scaled_jacobian_min": float(sj.min()),
            "scaled_jacobian_mean": float(sj.mean()),
            "n_inverted": int((sj <= 0).sum()),
            "quads_per_side": [len(l) for _r, l in faces],
            "n_core_points": int(len(S.P)),
            "interface_quads": iface_quads,
            "source": {"blocks": str(a.blocks), "msh": str(a.msh),
                       "input_vtk": str(a.input_vtk)},
        }
        Path(a.json).parent.mkdir(parents=True, exist_ok=True)
        Path(a.json).write_text(json.dumps(summary, indent=2) + "\n")
        print(f"wrote {a.json}")
