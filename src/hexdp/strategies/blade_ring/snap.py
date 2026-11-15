"""B2: snap the boundary vertices of the core block mesh onto the exact dtOO geometry.

    snap.py <run_dir> [--blocks blocks.vtk] [--out blade/blocks_b2.vtk] [--json blade/b2_snap.json]

Every boundary vertex of blocks.vtk gets the set of tet.vtk surface labels of its
boundary quads (nearest input triangle, the clean_blocks labeller rule):
  1 surface   -> closest point on that surface,
  2 surfaces  -> closest point on the dtOO feature curves (1D .msh entities) that
                 lie on both surfaces; alternating projection if there is none,
  3+ surfaces -> the curve endpoint (dtOO corner) that lies on all of them.
The cutout (label 7) is snapped too: the ring is built on it later (B3), so its
nodes must sit on the exact outer surface of the removed O-grid.

Interior vertices follow by a harmonic extension of the boundary displacement (mesh
morphing), then clean_blocks.untangle (interior vertices only, boundary bit-identical).
"""
import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


from hexdp.core import clean_blocks as cb
from hexdp.strategies.blade_ring.geometry import Geometry, LABEL_SURF, VOL_FACES
from hexdp.quality import hexq as qual

# a msh node lies on a reference surface: exactly for the msh facets, to the
# sampling chord (<= 5e-6) for the analytic dtOO hub / shroud
ON_TOL = 1e-5


def curve_surfaces(G):
    """Which reference surfaces each dtOO feature curve lies on (all its nodes)."""
    out = {}
    for t, c in G.curves.items():
        out[t] = {n for n in G.surfaces
                  if np.linalg.norm(G.project(n, c.P) - c.P, axis=1).max() < ON_TOL}
    return out


def load_core(run, blocks, G=None):
    import meshio
    m = meshio.read(str(blocks))
    P = np.asarray(m.points, float)
    H = np.asarray(m.cells_dict["hexahedron"], int)
    B = np.asarray(m.cell_data["block_id"][0]).ravel().astype(int)
    f2h = defaultdict(list)
    quads = {}
    for i, h in enumerate(H):
        for f in VOL_FACES["hexahedron"]:
            q = tuple(int(x) for x in h[list(f)])
            k = tuple(sorted(q))
            f2h[k].append(i)
            quads[k] = q
    t = meshio.read(str(Path(run) / "tet.vtk"))
    tri, lab = [], []
    for c, d in zip(t.cells, t.cell_data["color"]):
        if c.type == "triangle":
            tri.append(c.data)
            lab.append(np.asarray(d).ravel())
    tri, lab = np.vstack(tri), np.concatenate(lab).astype(int)
    kd = cKDTree(t.points[tri].mean(1))
    bkeys = [k for k, hs in f2h.items() if len(hs) == 1]
    _, ii = kd.query(np.array([P[list(k)].mean(0) for k in bkeys]))
    qlab = {k: int(lab[i]) for k, i in zip(bkeys, ii)}
    # One label per boundary BLOCK SIDE (majority of its quads), not per quad: beam
    # leaves the boundary up to its Hausdorff tolerance (0.1) off the input surface,
    # and there the nearest-triangle rule mislabels single quads (0005: a hub quad at
    # z = 1.77 read as the O-grid cut, which then dragged its vertex 7 cm).
    from hexdp.fill import tfi as tfi
    lat, missing = tfi.lattices(P, H, B, {frozenset(k): v for k, v in f2h.items()}, verbose=False)
    if missing:
        raise SystemExit(f"[b2] blocks without lattice: {missing}")
    vlab = defaultdict(set)
    patches = []
    for blk, (dims, vert) in lat.items():
        for axis in (0, 1, 2):
            for side in (0, 1):
                g = tfi.face_grid(vert, axis, side)
                ks = [tuple(sorted(int(x) for x in (g[i, j], g[i + 1, j], g[i + 1, j + 1], g[i, j + 1])))
                      for i in range(g.shape[0] - 1) for j in range(g.shape[1] - 1)]
                labs = [qlab[k] for k in ks if k in qlab]
                if not labs:
                    continue
                if len(labs) != len(ks):
                    raise SystemExit(f"[b2] block {blk} side {axis}/{side} only partly on the boundary")
                vals, cnt = np.unique(labs, return_counts=True)
                major = int(vals[np.argmax(cnt)])
                rec = {"block": int(blk), "side": [axis, side], "label": major,
                       "agree": float(cnt.max() / len(labs)), "quads": len(labs)}
                if G is not None and len(vals) > 1:
                    # mixed side: take the candidate surface the whole side lies closest to
                    C = np.array([P[list(k)].mean(0) for k in ks])
                    dist = {int(l): float(np.mean(np.linalg.norm(G.project(LABEL_SURF[int(l)], C) - C, axis=1)))
                            for l in vals if int(l) in LABEL_SURF}
                    if dist:
                        major = min(dist, key=dist.get)
                        rec.update(label=major, mixed_dist={str(k): v for k, v in dist.items()})
                patches.append(rec)
                for v in np.unique(g):
                    vlab[int(v)].add(major)
    load_core.patches = patches
    return P, H, B, dict(f2h), vlab, m


def snap(G, P, vlab, csurf):
    Q = P.copy()
    kind = {}
    target = {}                                       # v -> ("surface", name) | ("curve", id)
    snap.target = target
    ends = {t: (c.P[0], c.P[-1]) for t, c in G.curves.items()}
    for v, labs in vlab.items():
        names = sorted({LABEL_SURF[l] for l in labs if l in LABEL_SURF})
        x = P[v]
        if not names:
            continue
        if len(names) == 1:
            Q[v] = G.project(names[0], x)[0]
            kind[v] = "surface"
            target[v] = ("surface", names[0])
            continue
        cs = [t for t, ss in csurf.items() if set(names) <= ss]
        if len(names) >= 3:
            # dtOO corner: a curve endpoint that lies on every named surface
            pts = np.array([e for t, ss in csurf.items() if len(set(names) & ss) >= 2 for e in ends[t]])
            if len(pts):
                on = np.all([np.linalg.norm(G.project(n, pts) - pts, axis=1) < ON_TOL for n in names], axis=0)
                pts = pts[on]
            if len(pts):
                Q[v] = pts[np.argmin(np.linalg.norm(pts - x, axis=1))]
                kind[v] = "corner"
                continue
        if cs:
            cand = np.array([G.project_curve(G.curves[t], names, x)[0] for t in cs])
            i = int(np.argmin(np.linalg.norm(cand - x, axis=1)))
            Q[v] = cand[i]
            kind[v] = "curve"
            target[v] = ("curve", cs[i])
            continue
        y = x.copy()                                  # no dtOO curve: alternating projection
        for _ in range(30):
            for n in names:
                y = G.project(n, y)[0]
        Q[v] = y
        kind[v] = "alt:" + "+".join(names)
    return Q, kind


def extend_displacement(P, Q, H, fixed):
    """Harmonic extension of the boundary displacement Q - P into the interior.

    Graph Laplacian over the hex edges with 1/length weights; boundary vertices are
    Dirichlet (their snapped displacement), interior vertices solve L d = 0. This is
    plain mesh morphing: it carries a large boundary correction (beam allows a
    Hausdorff of 0.1) into the volume instead of crushing the first cell layer."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import spsolve
    E = set()
    for h in H:
        for a, b in ((0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
                     (0, 4), (1, 5), (2, 6), (3, 7)):
            i, j = int(h[a]), int(h[b])
            E.add((min(i, j), max(i, j)))
    E = np.array(sorted(E))
    w = 1.0 / np.maximum(np.linalg.norm(P[E[:, 0]] - P[E[:, 1]], axis=1), 1e-12)
    n = len(P)
    W = coo_matrix((np.r_[w, w], (np.r_[E[:, 0], E[:, 1]], np.r_[E[:, 1], E[:, 0]])), shape=(n, n)).tocsr()
    L = (coo_matrix((np.asarray(W.sum(1)).ravel(), (np.arange(n), np.arange(n))), shape=(n, n)) - W).tocsr()
    fixed = np.zeros(n, bool) | fixed
    free = ~fixed
    D = Q - P
    out = Q.copy()
    if free.any():
        A = L[free][:, free]
        rhs = -L[free][:, fixed] @ D[fixed]
        out[free] = P[free] + np.column_stack([spsolve(A.tocsc(), rhs[:, k]) for k in range(3)])
    return out


def _frames(G, Q, target, free_int):
    """Per free vertex: its base point and the tangent basis it may move in
    (3 axes interior, 2 on a surface, 1 on a curve)."""
    out = []
    for v in free_int:
        out.append((v, np.eye(3)))
    for v, (k, n) in target.items():
        x = Q[v]
        if k == "surface":
            S = G.surfaces[n]
            _, i = S.kd.query(x)
            nrm = np.cross(S.B[i] - S.A[i], S.C[i] - S.A[i])
            nrm /= np.linalg.norm(nrm)
            t1 = np.cross(nrm, [1.0, 0, 0] if abs(nrm[0]) < 0.9 else [0, 1.0, 0])
            t1 /= np.linalg.norm(t1)
            out.append((v, np.array([t1, np.cross(nrm, t1)])))
        else:
            C = G.curves[n].P
            j = int(np.argmin(np.linalg.norm(0.5 * (C[1:] + C[:-1]) - x, axis=1)))
            t = C[j + 1] - C[j]
            out.append((v, (t / np.linalg.norm(t))[None, :]))
    return out


def slide_optimize(G, Q, H, target, interior, eps=0.35, rounds=4, verbose=True, csurf_names=None):
    """Raise the worst cells by moving interior vertices freely and boundary vertices
    TANGENTIALLY (surface: 2 dof, curve: 1 dof, corners fixed), re-projecting onto the
    exact surface/curve after every round. Objective = the one-sided barrier of
    clean_blocks.untangle, sum max(0, eps - sJ)^2, restricted to the 2-ring around the
    cells below eps. A round is kept only if (inverted, -min sJ) improves."""
    from scipy.optimize import minimize
    v2c = defaultdict(list)
    for i, h in enumerate(H):
        for v in h:
            v2c[int(v)].append(i)
    movable = set(interior) | set(target)
    g_best = None
    for rnd in range(rounds):
        sj = cb.scaled_jacobians(Q, H)
        g = (int((sj <= 0).sum()), -float(sj.min()))
        g_best = g if g_best is None else g_best
        low = np.where(sj < eps)[0]
        if not len(low):
            break
        seed = {int(v) for c in low for v in H[c]}
        seed = {int(v) for x in seed for c in v2c[x] for v in H[c]}
        free = sorted(v for v in seed if v in movable)
        fr = _frames(G, Q, {v: target[v] for v in free if v in target},
                     [v for v in free if v not in target])
        cells = sorted({c for v, _ in fr for c in v2c[v]})
        Hc = H[cells]
        sizes = [len(b) for _, b in fr]
        off = np.r_[0, np.cumsum(sizes)]
        h0 = np.array([min(np.linalg.norm(Q[H[c][1]] - Q[H[c][0]]) for c in v2c[v]) for v, _ in fr])
        bounds = [(-0.3 * h0[i], 0.3 * h0[i]) for i, (_, b) in enumerate(fr) for _ in range(len(b))]

        def place(x):
            X = Q.copy()
            for i, (v, b) in enumerate(fr):
                X[v] = Q[v] + x[off[i]:off[i + 1]] @ b
            return X

        def obj(x):
            s = cb.scaled_jacobians(place(x), Hc)
            return float((np.maximum(0.0, eps - s) ** 2).sum())

        res = minimize(obj, np.zeros(off[-1]), method="L-BFGS-B", bounds=bounds,
                       options={"maxiter": 300})
        X = place(res.x)
        for v, (k, n) in target.items():              # back onto the exact geometry
            if v in free:
                X[v] = G.project(n, X[v])[0] if k == "surface" else \
                    G.project_curve(G.curves[n], sorted(csurf_names[n]), X[v])[0]
        s2 = cb.scaled_jacobians(X, H)
        g2 = (int((s2 <= 0).sum()), -float(s2.min()))
        if verbose:
            print(f"[b2] slide round {rnd}: {len(free)} vertices ({sum(1 for v in free if v in target)} on the "
                  f"boundary), {len(cells)} cells: inverted {g[0]} -> {g2[0]}, min sJ {-g[1]:.4f} -> {-g2[1]:.4f}")
        if g2 < g:
            Q = X
            g_best = g2
        else:
            break
    return Q


def quality(P, H):
    """quality.summary (inverted = folded, concave reported apart) + corner sJ."""
    return qual.summary(P, H), cb.scaled_jacobians(P, H)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run")
    ap.add_argument("--blocks", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--json", default=None)
    ap.add_argument("--eps", type=float, default=0.35, help="barrier target of the slide optimisation")
    a = ap.parse_args()
    run = Path(a.run)
    blocks = Path(a.blocks) if a.blocks else run / "blocks.vtk"
    out = Path(a.out) if a.out else run / "blade" / "blocks_b2.vtk"
    js = Path(a.json) if a.json else run / "blade" / "b2_snap.json"
    G = Geometry(run)
    csurf = curve_surfaces(G)
    P, H, B, f2h, vlab, m = load_core(run, blocks, G)
    q0, _ = quality(P, H)
    Q, kind = snap(G, P, vlab, csurf)
    disp = np.linalg.norm(Q - P, axis=1)
    qs, _ = quality(Q, H)                              # boundary snapped, interior untouched
    # hierarchical morphing: corners/curves -> surfaces -> interior. A surface vertex
    # first follows the displacement its corners and curves imply (harmonic over the
    # mesh graph), then is projected onto its surface; only then the interior follows.
    # Snapping surface vertices straight to the nearest point shears the cells next
    # to a face that beam left 0.1 m inside the domain (0005: outlet at z = 2.38).
    fixed = np.zeros(len(P), bool)
    fixed[[v for v, k in kind.items() if k in ("corner", "curve")]] = True
    Qc = extend_displacement(P, np.where(fixed[:, None], Q, P), H, fixed)
    for v, k in kind.items():
        if k == "surface":
            n = LABEL_SURF[[l for l in vlab[v] if l in LABEL_SURF][0]]
            Q[v] = G.project(n, Qc[v])[0]
    fixed[list(vlab)] = True
    Q = extend_displacement(P, Q, H, fixed)
    q1, _ = quality(Q, H)
    Q2, moved = cb.untangle(Q, H, f2h, labeller=None, rings=2, verbose=True)
    qu, _ = quality(Q2, H)
    interior = [v for v in range(len(P)) if v not in vlab]
    Q2 = slide_optimize(G, Q2, H, snap.target, interior, eps=a.eps,
                        csurf_names=csurf)
    q2, sj = quality(Q2, H)
    resid = 0.0                                       # final vertex vs every surface it is labelled with
    for v, labs in vlab.items():
        for l in labs:
            if l in LABEL_SURF:
                resid = max(resid, float(np.linalg.norm(G.project(LABEL_SURF[l], Q2[v])[0] - Q2[v])))
    bdrift = max(float(np.linalg.norm(Q2[v] - Q[v])) for v in vlab)   # tangential slide length
    kinds = defaultdict(int)
    for k in kind.values():
        kinds[k] += 1
    rep = {"blocks_in": str(blocks), "blocks_out": str(out),
           "vertices": int(len(P)), "boundary_vertices": len(vlab),
           "snap_kinds": dict(kinds),
           "snap_displacement": {"max": float(disp.max()), "p95": float(np.percentile(disp[list(vlab)], 95))},
           "curves_on_surfaces": {int(t): sorted(s) for t, s in csurf.items()},
           "patches": load_core.patches,
           "vertex_surfaces": {str(v): sorted({LABEL_SURF[l] for l in labs if l in LABEL_SURF})
                               for v, labs in vlab.items()},
           "patches_label_agree_min": min(p["agree"] for p in load_core.patches),
           "quality_before": q0, "quality_snap_only": qs, "quality_after_snap": q1, "quality_after_untangle": qu,
           "quality_after_slide": q2,
           "residual_on_surfaces_max": resid,
           "untangle_moved_vertices": int(moved), "boundary_drift_after_untangle": bdrift,
           "gate": {"inverted_0": q2["inverted"] == 0, "concave_0": q2["concave"] == 0,
                    "min_sj_ge_0.2": q2["min_sj"] >= 0.2,
                    "residual_lt_1e-6": resid < 1e-6}}
    import meshio
    pd = {"snap_disp": disp, "snap_kind": np.array([{"surface": 1, "curve": 2, "corner": 3}.get(kind.get(v, ""), 4 if v in kind else 0) for v in range(len(P))], float)}
    meshio.write(str(out), meshio.Mesh(Q2, [("hexahedron", H)], point_data=pd,
                                       cell_data={"block_id": [B], "sj": [sj]}))
    js.write_text(json.dumps(rep, indent=1))
    print(json.dumps({k: rep[k] for k in ("snap_kinds", "snap_displacement", "residual_on_surfaces_max", "quality_before",
                                          "quality_after_snap", "quality_after_untangle", "quality_after_slide", "gate")}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
