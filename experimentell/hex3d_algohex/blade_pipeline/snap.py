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

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))

import clean_blocks as cb                                            # noqa: E402
from geometry import Geometry, LABEL_SURF, VOL_FACES                 # noqa: E402

ON_TOL = 1e-9


def curve_surfaces(G):
    """Which reference surfaces each dtOO feature curve lies on (all its nodes)."""
    kds = {n: cKDTree(s.P[np.unique(s.T)]) for n, s in G.surfaces.items()}
    out = {}
    for t, c in G.curves.items():
        out[t] = {n for n, kd in kds.items() if kd.query(c.P)[0].max() < ON_TOL}
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
    import tfi
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
    ends = {t: (c.P[0], c.P[-1]) for t, c in G.curves.items()}
    for v, labs in vlab.items():
        names = sorted({LABEL_SURF[l] for l in labs if l in LABEL_SURF})
        x = P[v]
        if not names:
            continue
        if len(names) == 1:
            Q[v] = G.project(names[0], x)[0]
            kind[v] = "surface"
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
            cand = np.array([G.curves[t].project(x)[0] for t in cs])
            Q[v] = cand[np.argmin(np.linalg.norm(cand - x, axis=1))]
            kind[v] = "curve"
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


def quality(P, H):
    sj = cb.scaled_jacobians(P, H)
    return {"inverted": int((sj <= 0).sum()), "min_sj": float(sj.min()),
            "p05_sj": float(np.percentile(sj, 5)), "mean_sj": float(sj.mean())}, sj


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run")
    ap.add_argument("--blocks", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--json", default=None)
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
    resid = 0.0                                       # snapped vertex vs every surface it is labelled with
    for v, labs in vlab.items():
        for l in labs:
            if l in LABEL_SURF:
                resid = max(resid, float(np.linalg.norm(G.project(LABEL_SURF[l], Q[v])[0] - Q[v])))
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
    q2, sj = quality(Q2, H)
    bdrift = max(float(np.linalg.norm(Q2[v] - Q[v])) for v in vlab)
    kinds = defaultdict(int)
    for k in kind.values():
        kinds[k] += 1
    rep = {"blocks_in": str(blocks), "blocks_out": str(out),
           "vertices": int(len(P)), "boundary_vertices": len(vlab),
           "snap_kinds": dict(kinds),
           "snap_displacement": {"max": float(disp.max()), "p95": float(np.percentile(disp[list(vlab)], 95))},
           "curves_on_surfaces": {int(t): sorted(s) for t, s in csurf.items()},
           "patches": load_core.patches,
           "patches_label_agree_min": min(p["agree"] for p in load_core.patches),
           "quality_before": q0, "quality_snap_only": qs, "quality_after_snap": q1, "quality_after_untangle": q2,
           "residual_on_surfaces_max": resid,
           "untangle_moved_vertices": int(moved), "boundary_drift_after_untangle": bdrift,
           "gate": {"inverted_0": q2["inverted"] == 0, "min_sj_ge_0.2": q2["min_sj"] >= 0.2,
                    "residual_lt_1e-6": resid < 1e-6}}
    import meshio
    pd = {"snap_disp": disp, "snap_kind": np.array([{"surface": 1, "curve": 2, "corner": 3}.get(kind.get(v, ""), 4 if v in kind else 0) for v in range(len(P))], float)}
    meshio.write(str(out), meshio.Mesh(Q2, [("hexahedron", H)], point_data=pd,
                                       cell_data={"block_id": [B], "sj": [sj]}))
    js.write_text(json.dumps(rep, indent=1))
    print(json.dumps({k: rep[k] for k in ("snap_kinds", "snap_displacement", "residual_on_surfaces_max", "quality_before",
                                          "quality_after_snap", "quality_after_untangle", "gate")}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
