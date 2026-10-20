"""Can the periodic boundary of a tistos block sample be made conforming in post-processing?

Per sample.npz (block structure after beam collapse, as in the dataset):
 1. combinatorics: block-face layouts on periodic_A and periodic_B isomorphic (same quads, same
    arrangement, same neighbouring patches along the border), orientation consistent with the
    rotation by the pitch
 2. divisions: corresponding block edges in direction classes with equal counts
 3. geometry: corner displacement needed (B := R.A, symmetric average projected back onto the
    surfaces) and the corner-based scaled Jacobian of every block before/after

usage: periodic_test.py <out.json> <sample.npz>...
"""
import json, sys
from collections import Counter, defaultdict
from multiprocessing import Pool

import numpy as np
import networkx as nx
from networkx.algorithms import isomorphism as iso
from scipy.spatial import cKDTree

sys.path.insert(0, "/opt/stack/domain_partition_3D/experimentell/hex3d_algohex")
import curve_postprocess as cp          # noqa: E402

PER_A, PER_B, NB = 3, 4, 4
CORNERS = [(0, 1, 3, 4), (1, 2, 0, 5), (2, 3, 1, 6), (3, 0, 2, 7),
           (4, 7, 5, 0), (5, 4, 6, 1), (6, 5, 7, 2), (7, 6, 4, 3)]


def rot(a):
    return np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])


def closest_on(X, pts, tris, k=16):
    """closest points of X on the triangle set (exact over the k nearest centroids)"""
    A, B, C = (pts[tris[:, i]] for i in range(3))
    _, cand = cKDTree((A + B + C) / 3).query(X, k=min(k, len(tris)))
    best, bq = np.full(len(X), np.inf), X.copy()
    for j in range(cand.shape[1]):
        a, b, c = A[cand[:, j]], B[cand[:, j]], C[cand[:, j]]
        q = closest_tri(X, a, b, c)
        d = np.linalg.norm(X - q, axis=1)
        m = d < best
        best[m], bq[m] = d[m], q[m]
    return bq


def closest_tri(p, a, b, c):
    ab, ac, ap = b - a, c - a, p - a
    d1, d2 = (ab * ap).sum(1), (ac * ap).sum(1)
    bp = p - b; d3, d4 = (ab * bp).sum(1), (ac * bp).sum(1)
    cpp = p - c; d5, d6 = (ab * cpp).sum(1), (ac * cpp).sum(1)
    va, vb, vc = d3 * d6 - d5 * d4, d5 * d2 - d1 * d6, d1 * d4 - d3 * d2
    den = va + vb + vc; den = np.where(np.abs(den) < 1e-300, 1e-300, den)
    q = a + ab * (vb / den)[:, None] + ac * (vc / den)[:, None]
    sd = lambda x, y: np.where(np.abs(y) < 1e-300, 0.0, x / np.where(np.abs(y) < 1e-300, 1, y))
    m = (vc <= 0) & (d1 >= 0) & (d3 <= 0); q[m] = (a + ab * sd(d1, d1 - d3)[:, None])[m]
    m = (vb <= 0) & (d2 >= 0) & (d6 <= 0); q[m] = (a + ac * sd(d2, d2 - d6)[:, None])[m]
    m = (va <= 0) & (d4 - d3 >= 0) & (d5 - d6 >= 0)
    q[m] = (b + (c - b) * sd(d4 - d3, (d4 - d3) + (d5 - d6))[:, None])[m]
    m = (d1 <= 0) & (d2 <= 0); q[m] = a[m]
    m = (d3 >= 0) & (d4 <= d3); q[m] = b[m]
    m = (d6 >= 0) & (d5 <= d6); q[m] = c[m]
    return q


def block_sj(V, blocks):
    out = []
    for b in blocks:
        P8 = V[b]
        s = []
        for a, i, j, k in CORNERS:
            u, v, w = P8[i] - P8[a], P8[j] - P8[a], P8[k] - P8[a]
            n = np.linalg.norm(u) * np.linalg.norm(v) * np.linalg.norm(w)
            s.append(np.dot(u, np.cross(v, w)) / n if n > 0 else 0.0)
        out.append(min(s))
    return np.array(out)


def quad_labels(s):
    """boundary quads and their patch (as in curve_postprocess.postprocess)"""
    V0 = np.asarray(s["vertices"], float)
    blocks = np.asarray(s["blocks"], np.int64)
    pts, tris, lab = (np.asarray(s[k]) for k in ("surface_points", "surface_tris", "surface_tri_label"))
    E, ctrl, _ = cp.undirected(s)
    faces, nown = cp.faces_of(blocks)
    bq = [q for q, n in zip(faces, nown) if n == 1]
    tol = cp.VERTEX_TOL * np.linalg.norm(np.ptp(pts, 0))
    Dv, labs = cp.patch_dists(V0, pts, tris, lab)
    vpatch = [{labs[j] for j in np.flatnonzero(Dv[i] < tol)} for i in range(len(V0))]
    eid = {tuple(e): k for k, e in enumerate(map(tuple, E))}
    qlab = []
    for q in bq:
        common = set.intersection(*(vpatch[v] for v in q))
        if len(common) == 1:
            qlab.append(common.pop()); continue
        mids = [cp.bezier(V0[E[k][0]], ctrl[k, 0], ctrl[k, 1], V0[E[k][1]], 3)[1]
                for k in (eid[tuple(sorted((q[i], q[(i + 1) % 4])))] for i in range(4))]
        D, _ = cp.patch_dists(np.mean(mids, 0)[None], pts, tris, lab)
        cand = [j for j, l in enumerate(labs) if l in common] or list(range(len(labs)))
        qlab.append(labs[min(cand, key=lambda j: D[0, j])])
    return V0, blocks, pts, tris, lab, bq, qlab, vpatch


def layout_graph(quads, attr):
    G = nx.Graph()
    for fi, q in enumerate(quads):
        G.add_node(("f", fi), t="f", a=None)
        for i, v in enumerate(q):
            G.add_node(("v", v), t="v", a=attr[v])
            G.add_edge(("f", fi), ("v", v))
            G.add_edge(("v", v), ("v", q[(i + 1) % 4]))
    return G


def run(path):
    s = dict(np.load(path, allow_pickle=True))
    name = path.split("/")[-2]
    V0, blocks, pts, tris, lab, bq, qlab, vpatch = quad_labels(s)
    qa = [q for q, l in zip(bq, qlab) if l == PER_A]
    qb = [q for q, l in zip(bq, qlab) if l == PER_B]
    inc = defaultdict(set)
    for q, l in zip(bq, qlab):
        for v in q:
            inc[v].add(l)
    attr = {v: frozenset(inc[v] - {PER_A, PER_B}) for v in inc}
    res = {"sample": name, "n_quads_A": len(qa), "n_quads_B": len(qb),
           "n_verts_A": len({v for q in qa for v in q}), "n_verts_B": len({v for q in qb for v in q})}
    # rotation that maps the A surface onto the B surface
    SA, SB = pts[np.unique(tris[lab == PER_A])], pts[np.unique(tris[lab == PER_B])]
    tb = cKDTree(SB)
    R = min((rot(sg * 2 * np.pi / NB) for sg in (1, -1)), key=lambda R: np.median(tb.query(SA @ R.T)[0]))
    GA, GB = layout_graph(qa, attr), layout_graph(qb, attr)
    gm = iso.GraphMatcher(GA, GB, node_match=lambda x, y: x["t"] == y["t"] and x["a"] == y["a"])
    best = None
    for k, m in enumerate(gm.isomorphisms_iter()):
        vm = {a[1]: b[1] for a, b in m.items() if a[0] == "v"}
        cost = sum(np.linalg.norm(R @ V0[a] - V0[b]) for a, b in vm.items())
        if best is None or cost < best[0]:
            best = (cost, vm, m)
        if k > 2000:
            break
    res["isomorphic"] = best is not None
    if best is None:
        res["attr_A"] = sorted(Counter(tuple(sorted(attr[v])) for q in qa for v in set(q)).items())
        res["attr_B"] = sorted(Counter(tuple(sorted(attr[v])) for q in qb for v in set(q)).items())
        return res
    _, vm, m = best
    # orientation: R maps the outward normal of A to the inward side of B
    sgn = []
    for fi, q in enumerate(qa):
        qq = qb[m[("f", fi)][1]]
        nA = np.cross(V0[q[2]] - V0[q[0]], V0[q[3]] - V0[q[1]])
        nB = np.cross(V0[qq[2]] - V0[qq[0]], V0[qq[3]] - V0[qq[1]])
        sgn.append(np.sign(np.dot(R @ nA, nB)))
    res["orientation_consistent"] = len(set(sgn)) == 1
    # divisions of corresponding edges
    cls = {}
    for (u, v), c in zip(np.asarray(s["edges"]), np.asarray(s["dir_class"])):
        cls[(int(u), int(v))] = int(c)
    cnt = np.asarray(s["dir_class_count"])
    sides = {(q[i], q[(i + 1) % 4]) for q in qa for i in range(4)}
    mism = [(cnt[cls[e]], cnt[cls[(vm[e[0]], vm[e[1]])]]) for e in sides]
    res["edges_A"] = len(sides)
    res["division_mismatch"] = int(sum(a != b for a, b in mism))
    # geometry: displacement, symmetric fix, block quality
    elen = np.median([np.linalg.norm(V0[a] - V0[b]) for a, b in cls])
    disp = np.array([np.linalg.norm(R @ V0[a] - V0[b]) for a, b in vm.items()])
    res["disp_rel_median"], res["disp_rel_max"] = float(np.median(disp) / elen), float(disp.max() / elen)
    V1 = V0.copy()
    trA = tris[lab == PER_A]
    for a, b in vm.items():
        p = 0.5 * (V0[a] + R.T @ V0[b])
        others = [l for l in attr[a]]
        for _ in range(20):                       # alternating projection onto A and its border patches
            p = closest_on(p[None], pts, trA)[0]
            for l in others:
                p = closest_on(p[None], pts, tris[lab == l])[0]
        p = closest_on(p[None], pts, trA)[0]
        V1[a], V1[b] = p, R @ p
    s0, s1 = block_sj(V0, blocks), block_sj(V1, blocks)
    res.update(sj_min_before=float(s0.min()), sj_min_after=float(s1.min()),
               inv_before=int((s0 <= 0).sum()), inv_after=int((s1 <= 0).sum()),
               periodic_residual_after=float(max(np.linalg.norm(R @ V1[a] - V1[b]) for a, b in vm.items())))
    # curved block edges: drag every edge with its end displacements, then make B = R.A on the
    # periodic surface (A and rotated-back B averaged)
    cur0 = edge_curves(s)
    dV = V1 - V0
    cur1 = {}
    for (u, v), pl in cur0.items():
        sl = np.r_[0, np.cumsum(np.linalg.norm(np.diff(pl, axis=0), axis=1))]
        t = (sl / sl[-1])[:, None] if sl[-1] > 0 else np.zeros((len(pl), 1))
        cur1[(u, v)] = pl + (1 - t) * dV[u] + t * dV[v]
    sideA = {(q[i], q[(i + 1) % 4]) for q in qa for i in range(4)}
    for (u, v) in list(sideA) + [(v, u) for (u, v) in sideA]:
        bu, bv = vm[u], vm[v]
        n = 64
        a_ = resample(cur1[(u, v)], n)
        b_ = resample(cur1[(bu, bv)], n) @ R
        avg = 0.5 * (a_ + b_)
        avg[0], avg[-1] = V1[u], V1[v]
        cur1[(u, v)], cur1[(bu, bv)] = avg, avg @ R.T
    a0, a1, all0, all1 = [], [], [], []
    for b in blocks:
        keys = [(int(b[x]), int(b[y])) for x, y in XE + YE + ZE]
        j0 = cells_sj(block_tfi(b, cur0, 8)); j1 = cells_sj(block_tfi(b, cur1, 8))
        all0.append(j0); all1.append(j1)
        if max(np.abs(resample(cur0[k], 8) - resample(cur1[k], 8)).max() for k in keys) > 1e-9:
            a0.append(j0); a1.append(j1)
    a0, a1 = np.concatenate(a0), np.concatenate(a1)
    all0, all1 = np.concatenate(all0), np.concatenate(all1)
    res.update(n_blocks_changed=len(a0) // 512,
               tfi_changed_minsj_before=float(a0.min()), tfi_changed_minsj_after=float(a1.min()),
               tfi_changed_mean_before=float(a0.mean()), tfi_changed_mean_after=float(a1.mean()),
               tfi_changed_inv_before=int((a0 <= 0).sum()), tfi_changed_inv_after=int((a1 <= 0).sum()),
               tfi_minsj_before=float(all0.min()), tfi_minsj_after=float(all1.min()),
               tfi_inv_before=int((all0 <= 0).sum()), tfi_inv_after=int((all1 <= 0).sum()))
    return res


# ---------------------------------------------------------------- TFI quality check
XE = [(0, 1), (3, 2), (4, 5), (7, 6)]
YE = [(0, 3), (1, 2), (4, 7), (5, 6)]
ZE = [(0, 4), (1, 5), (3, 7), (2, 6)]


def resample(pl, n):
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(pl, axis=0), axis=1))]
    t = np.linspace(0, s[-1], n + 1)
    return np.stack([np.interp(t, s, pl[:, k]) for k in range(3)], 1)


def edge_curves(s):
    E, off, poly = np.asarray(s["edges"]), np.asarray(s["edge_polyline_offset"]), np.asarray(s["edge_polyline"])
    return {(int(u), int(v)): poly[off[i]:off[i + 1]] for i, (u, v) in enumerate(E)}


def block_tfi(b, curves, n):
    """edge-based Gordon-Hall TFI of one block from its 12 curved edges -> [n+1]^3 points"""
    t = np.linspace(0, 1, n + 1)
    C = {k: resample(curves[(int(b[a]), int(b[c]))], n) for k, (a, c) in enumerate(XE + YE + ZE)}
    X = np.zeros((n + 1, n + 1, n + 1, 3))
    u, v, w = np.meshgrid(t, t, t, indexing="ij")
    U, Vv, W = u[..., None], v[..., None], w[..., None]
    # x-edges at (v,w) = (0,0),(1,0),(0,1),(1,1); y-edges at (u,w); z-edges at (u,v)
    ex = [C[0][:, None, None], C[1][:, None, None], C[2][:, None, None], C[3][:, None, None]]
    X += ((1 - Vv) * (1 - W) * C[0][:, None, None] + Vv * (1 - W) * C[1][:, None, None]
          + (1 - Vv) * W * C[2][:, None, None] + Vv * W * C[3][:, None, None])
    X += ((1 - U) * (1 - W) * C[4][None, :, None] + U * (1 - W) * C[5][None, :, None]
          + (1 - U) * W * C[6][None, :, None] + U * W * C[7][None, :, None])
    X += ((1 - U) * (1 - Vv) * C[8][None, None, :] + U * (1 - Vv) * C[9][None, None, :]
          + (1 - U) * Vv * C[10][None, None, :] + U * Vv * C[11][None, None, :])
    P = None
    corners = {}
    for i in (0, 1):
        for j in (0, 1):
            for k in (0, 1):
                idx = [0, 1, 2, 3, 4, 5, 6, 7][[(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                                                 (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)].index((i, j, k))]
                corners[(i, j, k)] = C[0][0] if False else None
    # subtract 2 x trilinear corner interpolation (Gordon-Hall edge form)
    cp_ = {}
    cp_[(0, 0, 0)], cp_[(1, 0, 0)] = C[0][0], C[0][-1]
    cp_[(0, 1, 0)], cp_[(1, 1, 0)] = C[1][0], C[1][-1]
    cp_[(0, 0, 1)], cp_[(1, 0, 1)] = C[2][0], C[2][-1]
    cp_[(0, 1, 1)], cp_[(1, 1, 1)] = C[3][0], C[3][-1]
    tri = sum(((U if i else 1 - U) * (Vv if j else 1 - Vv) * (W if k else 1 - W)) * cp_[(i, j, k)]
              for i in (0, 1) for j in (0, 1) for k in (0, 1))
    return X - 2 * tri


def cells_sj(X):
    n = X.shape[0] - 1
    c = [X[:-1, :-1, :-1], X[1:, :-1, :-1], X[1:, 1:, :-1], X[:-1, 1:, :-1],
         X[:-1, :-1, 1:], X[1:, :-1, 1:], X[1:, 1:, 1:], X[:-1, 1:, 1:]]
    out = np.full(X[:-1, :-1, :-1, 0].shape, np.inf)
    for a, i, j, k in CORNERS:
        u, v, w = c[i] - c[a], c[j] - c[a], c[k] - c[a]
        nn = np.linalg.norm(u, axis=-1) * np.linalg.norm(v, axis=-1) * np.linalg.norm(w, axis=-1)
        out = np.minimum(out, np.einsum("...i,...i", u, np.cross(v, w)) / np.maximum(nn, 1e-300))
    return out.ravel()


def tfi_quality(s, V, curves, n=8):
    sj = np.concatenate([cells_sj(block_tfi(b, curves, n)) for b in np.asarray(s["blocks"])])
    return float(sj.min()), int((sj <= 0).sum()), float(np.percentile(sj, 1))


if __name__ == "__main__":
    out, files = sys.argv[1], sys.argv[2:]
    with Pool(16) as pool:
        rows = []
        for r in pool.imap_unordered(run, files, chunksize=4):
            rows.append(r)
    json.dump(rows, open(out, "w"), indent=1, default=str)
    print("done", len(rows))
