"""B3: the blade O-ring as block structure on the snapped core (blocks_b2.vtk).

    ring.py <run_dir> [--layers 2] [--out blade/blocks_b3.vtk] [--json blade/b3_ring.json]

Topology (decided 2026-11-11): 12 core blocks + one closed O-layer between the
cutout (label 7) and the blade wall. The cutout is already split into one side
per core block by the core's own block edges; every such side gets exactly one
ring block (ids 12, 13, ...), whose outer face IS that core side (same vertex
ids -> conformal, no T-nodes) and whose inner face lies on the blade.

Wall nodes: every cutout node goes to its foot point on the exact dtOO blade
S(u, v) (blade/machine_scratch/machine_blade.npz), the offset pairing of v3g.
Near the trailing-edge cap the foot points crowd onto the TE line, so per
segment between two split nodes the chordwise parameter u is blended towards
the arc-length-proportional map, with the smallest weight that keeps every
wall spacing >= --min-spacing of the segment mean. Rows on hub / shroud keep
v = 0 / 1 exactly (the blade ends on the hub/shroud cylinders).

Inside, the ring gets --layers straight coarse layers (this is block structure;
the fine ring mesh with BL is B6). Cells are classified with quality.py.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))

import tfi                                                            # noqa: E402
import quality as qual                                                # noqa: E402
import clean_blocks as cb                                             # noqa: E402
from geometry import closest_on_tris                                 # noqa: E402

CUT_LABEL = 7          # tet_prep_v5.S_OGRID_IFACE


class BladeSurface:
    """The dtOO blade sample S[i, j] = S(U[i], V[j]) as triangles with (u, v), u
    periodic. Also any other closed structured surface given as arrays (the cutout)."""

    def __init__(self, npz=None, S=None, U=None, V=None):
        if npz is not None:
            d = np.load(npz, allow_pickle=True)
            S, U, V = d["S"], d["U"], d["V"]
        self.S, self.U, self.V = S, U, V
        nu, nv = self.S.shape[:2]
        I, J = np.meshgrid(np.arange(nu - 1), np.arange(nv - 1), indexing="ij")
        a = (I * nv + J).ravel()
        b, c, e = a + nv, a + nv + 1, a + 1
        self.T = np.r_[np.c_[a, b, c], np.c_[a, c, e]]
        self.X = self.S.reshape(-1, 3)
        uu, vv = np.meshgrid(self.U, self.V, indexing="ij")
        self.UV = np.c_[uu.ravel(), vv.ravel()]
        self.kd = cKDTree(self.X[self.T].mean(1))
        self.period = float(self.U[-1] - self.U[0])

    def foot(self, x, k=32):
        """Closest point on the blade and its (u, v)."""
        _, nn = self.kd.query(x, k=k)
        A, B, C = (self.X[self.T[nn, i]] for i in range(3))
        Q = closest_on_tris(x[None, :], A, B, C)
        i = int(np.argmin(((Q - x) ** 2).sum(1)))
        t = self.T[nn[i]]
        # barycentric of Q[i] in the triangle -> (u, v); u never straddles the seam
        # within one triangle because the grid has the seam as a row
        M = np.c_[self.X[t[1]] - self.X[t[0]], self.X[t[2]] - self.X[t[0]]]
        lam = np.linalg.lstsq(M, Q[i] - self.X[t[0]], rcond=None)[0]
        uv = self.UV[t[0]] + lam[0] * (self.UV[t[1]] - self.UV[t[0]]) + lam[1] * (self.UV[t[2]] - self.UV[t[0]])
        return Q[i], uv

    def edges_le_te(self, window=1.0):
        """u of the trailing and leading edge: the two windows of largest total
        turning of the u-lines (averaged over v); TE = the one further downstream
        (larger z, flow from inlet z_min to outlet z_max)."""
        t = np.diff(self.S, axis=0)
        t /= np.linalg.norm(t, axis=2, keepdims=True)
        tt = np.concatenate([t, t[:1]], axis=0)          # closed: turn across the seam too
        ang = np.arccos(np.clip((tt[1:] * tt[:-1]).sum(2), -1, 1)).mean(1)   # at U[1:]
        uu = self.U[1:]
        per = self.period
        dist = lambda a, b: np.abs((a - b + per / 2) % per - per / 2)
        score = np.array([ang[dist(uu, c) < window].sum() for c in uu])
        c1 = uu[np.argmax(score)]
        score2 = np.where(dist(uu, c1) > 4 * window, score, -1)
        c2 = uu[np.argmax(score2)]
        z = lambda c: self.S[np.argmin(dist(self.U, c)), :, 2].mean()
        te, le = (c1, c2) if z(c1) > z(c2) else (c2, c1)
        return float(te), float(le)

    def at(self, u, v):
        """Bilinear S(u, v); u periodic."""
        u = self.U[0] + np.mod(u - self.U[0], self.period)
        i = np.clip(np.searchsorted(self.U, u) - 1, 0, len(self.U) - 2)
        j = np.clip(np.searchsorted(self.V, v) - 1, 0, len(self.V) - 2)
        s = (u - self.U[i]) / (self.U[i + 1] - self.U[i])
        t = (v - self.V[j]) / (self.V[j + 1] - self.V[j])
        S = self.S
        return ((1 - s) * (1 - t))[..., None] * S[i, j] + (s * (1 - t))[..., None] * S[i + 1, j] + \
            ((1 - s) * t)[..., None] * S[i, j + 1] + (s * t)[..., None] * S[i + 1, j + 1]


def cut_sides(P, lat, patches):
    """Cutout block sides as (block, grid[span, chord]) with span ascending in r."""
    out = []
    for p in patches:
        if p["label"] != CUT_LABEL:
            continue
        g = tfi.face_grid(lat[p["block"]][1], *p["side"])
        r = np.hypot(P[g][..., 0], P[g][..., 1])
        if np.ptp(r, axis=0).mean() < np.ptp(r, axis=1).mean():
            g = g.T                                   # make axis 0 the span direction
        if np.hypot(*P[g[0, 0], :2]) > np.hypot(*P[g[-1, 0], :2]):
            g = g[::-1]
        out.append((p["block"], g))
    return out


def chain_sides(sides):
    """Order the sides around the blade: consecutive sides share a span line.
    Returns [(block, grid)] with each grid's chord axis running in loop direction."""
    rest = list(sides)
    loop = [rest.pop(0)]
    while rest:
        end = set(loop[-1][1][:, -1].tolist())
        for i, (b, g) in enumerate(rest):
            if set(g[:, 0].tolist()) == end:
                loop.append(rest.pop(i))
                break
            if set(g[:, -1].tolist()) == end:
                loop.append((b, g[:, ::-1]))
                rest.pop(i)
                break
        else:
            raise SystemExit("[b3] cutout sides do not form a closed loop")
    if set(loop[-1][1][:, -1].tolist()) != set(loop[0][1][:, 0].tolist()):
        raise SystemExit("[b3] cutout loop not closed")
    return loop


def _bb(X):
    """Blade-to-blade coordinates (z, r*theta) of 3D points."""
    r = np.hypot(X[..., 0], X[..., 1])
    return np.stack([X[..., 2], r.mean() * np.arctan2(X[..., 1], X[..., 0])], -1)


def quads_ok(cut, wall):
    """All quads (cut_j, cut_j+1, wall_j+1, wall_j) of one row segment have the
    same orientation in the blade-to-blade plane (no crossing offset lines)."""
    c, w = _bb(cut), _bb(np.asarray(wall))
    Q = np.stack([c[:-1], c[1:], w[1:], w[:-1]], 1)
    x, y = Q[..., 0], Q[..., 1]
    A = 0.5 * (x * np.roll(y, -1, 1) - np.roll(x, -1, 1) * y).sum(1)
    # every corner turns the same way as the quad (convex, so also not twisted)
    e1 = np.roll(Q, -1, 1) - Q
    e0 = Q - np.roll(Q, 1, 1)
    cr = e0[..., 0] * e1[..., 1] - e0[..., 1] * e1[..., 0]
    return bool((np.sign(cr) == np.sign(A)[:, None]).all() and len(set(np.sign(A))) == 1)


def blend_segment(u_foot, u_prop, min_frac, valid=None):
    """Smallest w in [0, 1] with monotone spacing >= min_frac * mean and, when given,
    valid(u) true (no crossing offset lines at the TE/LE cap)."""
    for w in np.linspace(0, 1, 21):
        u = (1 - w) * u_foot + w * u_prop
        d = np.diff(u)
        if (d > 0).all() and d.min() >= min_frac * d.mean() and (valid is None or valid(u)):
            return w, u
    return 1.0, u_prop


def optimize_ring(P, RC, bs, lid, Uf, Vf, hub_shroud, G=None, eps=0.3, verbose=True):
    """Block-level smoothing of the coarse ring: wall nodes slide ON the blade in
    (u, v) (u only on hub / shroud rows), inner layer nodes move freely (in the
    tangent plane of hub / shroud on those rows, re-projected at the end), cutout
    nodes (= core) stay. Objective: one-sided barrier sum max(0, eps - sJ)^2 plus
    the same barrier on the 24 centroid tets (so folded cells are pushed out too),
    over the ring cells. Kept only if (bad cells, -min sJ) improves."""
    from scipy.optimize import minimize
    L = lid.shape[0] - 1
    ns, nl = Uf.shape
    wall = lid[0].ravel()
    inner = lid[1:L].ravel()
    fixv = np.zeros((ns, nl), bool)
    fixv[0] |= bool(hub_shroud[0])
    fixv[-1] |= bool(hub_shroud[1])
    # inner nodes on the hub / shroud rows: 2 tangential dof on that surface
    onsurf = {}
    for row, name, flag in ((0, "hub", hub_shroud[0]), (ns - 1, "shroud", hub_shroud[1])):
        if flag:
            for l in range(1, L):
                for v in lid[l, row]:
                    onsurf[int(v)] = name
    Tb = {}
    for v, name in onsurf.items():
        S = G.surfaces[name]
        _, i = S.kd.query(P[v])
        n = np.cross(S.B[i] - S.A[i], S.C[i] - S.A[i])
        n /= np.linalg.norm(n)
        Tb[v] = n
    u0, v0 = Uf.ravel().copy(), Vf.ravel().copy()
    freev = ~fixv.ravel()
    sc = np.median(np.abs(np.diff(Uf, axis=1)))       # u step ~ one cell
    dv = np.median(np.abs(np.diff(Vf, axis=0)))
    Pin0 = P[inner].copy()
    h = np.median(np.linalg.norm(P[RC[:, 1]] - P[RC[:, 0]], axis=1))
    nu, nvf, ni = len(u0), int(freev.sum()), len(inner)

    def place(x):
        X = P.copy()
        u = u0 + sc * x[:nu]
        v = v0.copy()
        v[freev] += dv * x[nu:nu + nvf]
        X[wall] = bs.at(u, np.clip(v, bs.V[0], bs.V[-1]))
        D = h * x[nu + nvf:].reshape(-1, 3)
        for k, v in enumerate(inner):
            if int(v) in Tb:                       # drop the normal component
                n = Tb[int(v)]
                D[k] -= np.dot(D[k], n) * n
        X[inner] = Pin0 + D
        return X

    def reproject(X):
        for v, name in onsurf.items():
            X[v] = G.project(name, X[v])[0]
        return X

    def score(X):
        c = qual.classify(X, RC)
        return (int((~c["ok"]).sum()), -float(c["sj"].min()))

    def obj(x):
        X = place(x)
        sj = cb.scaled_jacobians(X, RC)
        t = qual.centroid_tets(X, RC)
        tv = t / np.abs(t).mean()
        return float((np.maximum(0, eps - sj) ** 2).sum() + (np.maximum(0, 0.05 - tv) ** 2).sum())

    x0 = np.zeros(nu + nvf + 3 * ni)
    g0 = score(place(x0))
    bnd = [(-0.45, 0.45)] * nu + [(-0.45, 0.45)] * nvf + [(-0.45, 0.45)] * (3 * ni)
    res = minimize(obj, x0, method="L-BFGS-B", bounds=bnd, options={"maxiter": 400})
    X = reproject(place(res.x))
    P = reproject(P.copy())
    g0 = score(P)
    g1 = score(X)
    if verbose:
        print(f"[b3] ring smoothing: bad cells {g0[0]} -> {g1[0]}, min sJ {-g0[1]:.4f} -> {-g1[1]:.4f}")
    return X if g1 < g0 else P


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run")
    ap.add_argument("--blocks", default=None, help="default blade/blocks_b2.vtk")
    ap.add_argument("--layers", type=int, default=2)
    ap.add_argument("--min-spacing", type=float, default=0.3)
    ap.add_argument("--eps", type=float, default=0.3, help="barrier target of the ring smoothing")
    ap.add_argument("--out", default=None)
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    run = Path(a.run)
    blocks = Path(a.blocks) if a.blocks else run / "blade" / "blocks_b2.vtk"
    out = Path(a.out) if a.out else run / "blade" / "blocks_b3.vtk"
    js = Path(a.json) if a.json else run / "blade" / "b3_ring.json"
    patches = json.loads((run / "blade" / "b2_snap.json").read_text())["patches"]

    P, H, B, f2h = tfi.load_blocks(str(blocks))
    lat, missing = tfi.lattices(P, H, B, f2h, verbose=False)
    loop = chain_sides(cut_sides(P, lat, patches))
    nspan = loop[0][1].shape[0]
    if any(g.shape[0] != nspan for _, g in loop):
        raise SystemExit("[b3] cutout sides disagree on the span count")
    # closed loop per span row: node ids, split positions (side starts)
    ids = np.concatenate([g[:, :-1] for _, g in loop], axis=1)          # (nspan, nloop)
    split = np.cumsum([0] + [g.shape[1] - 1 for _, g in loop])[:-1]
    nloop = ids.shape[1]

    bs = BladeSurface(run / "blade" / "machine_scratch" / "machine_blade.npz")
    U = np.zeros(ids.shape)
    V = np.zeros(ids.shape)
    gap = np.zeros(ids.shape)
    for i in range(nspan):
        for j in range(nloop):
            x = P[ids[i, j]]
            q, uv = bs.foot(x)
            U[i, j], V[i, j] = uv
            gap[i, j] = np.linalg.norm(q - x)
    # unwrap u along each row in loop direction, check the winding
    per = bs.period
    for i in range(nspan):
        U[i] = U[i, 0] + np.unwrap((U[i] - U[i, 0]) * 2 * np.pi / per) * per / (2 * np.pi)
    sgn = np.sign(U[:, -1] - U[:, 0]).mean()
    if sgn < 0:                                        # run u forward along the loop
        U = -U
    sgn_u = (lambda uu: -uu) if sgn < 0 else (lambda uu: uu)   # back to the blade's u
    weights = np.zeros((nspan, len(split)))
    Uf = U.copy()
    for i in range(nspan):
        row = np.r_[U[i], U[i, 0] + per]                # closed
        cut = P[np.r_[ids[i], ids[i, :1]]]
        vclosed = np.r_[V[i], V[i, :1]]
        seg_len = np.linalg.norm(np.diff(cut, axis=0), axis=1)
        bounds = list(split) + [nloop]

        def seg(k, rw):
            """blend segment k of row rw; (ok, w, u)."""
            j0, j1 = bounds[k], bounds[k + 1]
            uf = rw[j0:j1 + 1]
            fr = np.r_[0, np.cumsum(seg_len[j0:j1])] / seg_len[j0:j1].sum()
            up = uf[0] + fr * (uf[-1] - uf[0])
            vrow, cpts = vclosed[j0:j1 + 1], cut[j0:j1 + 1]
            ok = lambda uu: quads_ok(cpts, bs.at(sgn_u(uu), vrow))
            w, u = blend_segment(uf, up, a.min_spacing, valid=ok)
            return ok(u) and np.all(np.diff(u) > 0), w, u

        # The wall end of a ring edge (split node) is free: only its cutout end is
        # fixed by the core. Where the foot point of a split node lands on the wrong
        # side of the TE (the cap problem) no blend inside the segments can repair
        # it, so the split node's u is moved towards its neighbours, smallest move first.
        nsplit = len(split)
        for kk in range(nsplit):
            res = [seg(k, row) for k in range(nsplit)]
            bad = [k for k in range(nsplit) if not res[k][0]]
            if not bad:
                break
            for k in bad:
                for jn in (bounds[k], bounds[k + 1]):          # the two ring edges of segment k
                    jl, jr = (jn - 1) % nloop, (jn + 1) % nloop
                    u0 = row[jn]
                    lo = row[jl] - (per if jl > jn else 0)
                    hi = row[jr] + (per if jr < jn else 0)
                    best = None
                    for t in sorted(np.linspace(-0.9, 0.9, 37), key=abs):
                        cand = row.copy()
                        un = u0 + t * ((hi - u0) if t > 0 else (u0 - lo))
                        cand[jn] = un
                        if jn == 0:
                            cand[-1] = un + per
                        if jn == nloop:
                            cand[0] = un - per
                        kl, kr = (k - 1) % nsplit if jn == bounds[k] else k, k if jn == bounds[k] else (k + 1) % nsplit
                        if seg(kl, cand)[0] and seg(kr, cand)[0]:
                            best = cand
                            break
                    if best is not None:
                        row = best
                        break
        for k in range(nsplit):
            okk, w, u = seg(k, row)
            j0, j1 = bounds[k], bounds[k + 1]
            weights[i, k] = w
            row[j0:j1 + 1] = u
        Uf[i] = row[:-1]
    if sgn < 0:
        Uf = -Uf
    # hub / shroud rows end exactly on the blade's v = 0 / 1 edge
    # (detected by projection onto the reference hub / shroud: the facetted msh
    # surface sits ~2e-4 inside the analytic cylinder, so no radius test)
    from geometry import Geometry
    G = Geometry(run)
    on = lambda name, X: bool(np.all(np.linalg.norm(G.project(name, X) - X, axis=1) < 1e-8))
    hub_row, shroud_row = on("hub", P[ids[0]]), on("shroud", P[ids[-1]])
    Vf = V.copy()
    if hub_row:
        Vf[0] = bs.V[0]
    if shroud_row:
        Vf[-1] = bs.V[-1]
    W = bs.at(Uf, Vf)                                  # (nspan, nloop, 3) wall nodes

    # ring points: layer 0 = wall, layer L = cutout (existing core ids)
    L = a.layers
    newP = [P]
    lid = np.zeros((L + 1, nspan, nloop), int)
    lid[L] = ids
    base = len(P)
    for l in range(L):
        t = l / L
        X = (1 - t) * W + t * P[ids]
        lid[l] = base + np.arange(nspan * nloop).reshape(nspan, nloop)
        base += nspan * nloop
        newP.append(X.reshape(-1, 3))
    P2 = np.vstack(newP)
    ring_cells, ring_bid = [], []
    bounds = list(split) + [nloop]
    nb = int(B.max()) + 1
    for k in range(len(split)):
        cells = []
        for i in range(nspan - 1):
            for jj in range(bounds[k], bounds[k + 1]):
                j, j2 = jj % nloop, (jj + 1) % nloop
                for l in range(L):
                    cells.append([lid[l, i, j], lid[l, i, j2], lid[l, i + 1, j2], lid[l, i + 1, j],
                                  lid[l + 1, i, j], lid[l + 1, i, j2], lid[l + 1, i + 1, j2], lid[l + 1, i + 1, j]])
        cells = np.array(cells)
        vol = qual.centroid_tets(P2, cells).sum(1)
        if np.median(vol) < 0:                          # one orientation per block
            cells = cells[:, [4, 5, 6, 7, 0, 1, 2, 3]]
        ring_cells.append(cells)
        ring_bid += [nb + k] * len(cells)
    RC = np.vstack(ring_cells)
    P2 = optimize_ring(P2, RC, bs, lid, Uf, Vf, hub_shroud=(hub_row, shroud_row), G=G, eps=a.eps)
    W = P2[lid[0]]
    H2 = np.vstack([H, RC])
    B2 = np.r_[B, ring_bid]

    # checks
    core_cut = set(ids.ravel().tolist())
    ring_outer = set(lid[L].ravel().tolist())
    wall_res = []
    for x in W.reshape(-1, 3)[:: max(1, W.size // 3 // 300)]:
        q, _ = bs.foot(x)
        wall_res.append(float(np.linalg.norm(q - x)))
    ring_q = qual.summary(P2, RC)
    hs_dev = {}
    for row, name, flag in ((0, "hub", hub_row), (nspan - 1, "shroud", shroud_row)):
        if flag:
            X = P2[lid[1:, row].ravel()]
            hs_dev[name + "_inner_nodes"] = float(np.linalg.norm(G.project(name, X) - X, axis=1).max())
            X = P2[lid[0, row]]
            hs_dev[name + "_wall_nodes"] = float(np.linalg.norm(G.project(name, X) - X, axis=1).max())
    all_q = qual.summary(P2, H2)
    cls = qual.classify(P2, RC)
    worst = np.argsort(cls["sj"])[:5]
    # surfaces per vertex, by construction: core from B2 (the cutout is interior now),
    # ring wall = blade (+ hub / shroud on those rows), ring inner rows on hub / shroud
    vs = {}
    for v, names in json.loads((run / "blade" / "b2_snap.json").read_text())["vertex_surfaces"].items():
        n = [x for x in names if x != "cutout"]
        if n:
            vs[v] = n
    for l in range(L):
        for i in range(nspan):
            extra = ["hub"] if (i == 0 and hub_row) else ["shroud"] if (i == nspan - 1 and shroud_row) else []
            for v in lid[l, i]:
                names = (["blade"] if l == 0 else []) + extra
                if names:
                    vs[str(int(v))] = sorted(names)
    rep = {"blocks_in": str(blocks), "blocks_out": str(out),
           "vertex_surfaces": vs,
           "core_blocks": int(nb), "ring_blocks": list(range(nb, nb + len(split))),
           "ring_block_of_core_side": {str(nb + k): int(b) for k, (b, _) in enumerate(loop)},
           "span_rows": int(nspan), "loop_nodes": int(nloop), "split_positions": split.tolist(),
           "layers": L,
           "conformal_outer_face": core_cut == ring_outer,
           "t_nodes": len(core_cut ^ ring_outer),
           "wall_residual_max": max(wall_res),
           "hub_shroud_rows": {"hub": hub_row, "shroud": shroud_row, "deviation_max": hs_dev},
           "gap_cutout_blade": {"min": float(gap.min()), "median": float(np.median(gap)), "max": float(gap.max())},
           "u_blend_weight": {"per_segment_max": weights.max(0).tolist(), "per_row": weights.tolist()},
           "quality_ring": ring_q, "quality_all": all_q,
           "worst_ring_cells": [{"cell": int(c), "block": int(ring_bid[c]), "sj": float(cls["sj"][c]),
                                 "inverted": bool(cls["inverted"][c]), "concave": bool(cls["concave"][c])}
                                for c in worst],
           "gate": {"conformal": core_cut == ring_outer, "inverted_0": ring_q["inverted"] == 0,
                    "wall_on_blade_lt_1e-4": max(wall_res) < 1e-4}}   # bilinear S(u,v) vs its triangulation
    import meshio
    cl = qual.classify(P2, H2)
    part = np.r_[np.zeros(len(H)), np.ones(len(RC))]
    meshio.write(str(out), meshio.Mesh(P2, [("hexahedron", H2)],
                                       cell_data={"block_id": [B2], "sj": [cl["sj"]], "part": [part],
                                                  "inverted": [cl["inverted"].astype(float)],
                                                  "concave": [cl["concave"].astype(float)]}))
    js.write_text(json.dumps(rep, indent=1))
    print(json.dumps({k: rep[k] for k in ("ring_block_of_core_side", "span_rows", "loop_nodes", "t_nodes", "hub_shroud_rows",
                                          "wall_residual_max", "gap_cutout_blade", "quality_ring", "gate")},
                     indent=1))
    print("u-blend max per segment:", np.round(weights.max(0), 2).tolist())
    return 0


if __name__ == "__main__":
    sys.exit(main())
