"""B1: geometry oracle -- the exact dtOO surfaces of one run.

    geometry.py <run_dir> [--json out.json]

Reference surfaces come from the dtOO volume mesh case/mesh.msh, not from the
AlgoHex/beam output: the boundary faces of the volume mesh grouped by their
geometric entity. Physical tags are 0 on every element in the dtOO export, so
the entities are named the way tet_prep_v5 does it:
  - inlet / outlet / periodic_A / periodic_B: the tagged entities
    tet_prep_v5.GEOM_INLET / GEOM_OUTLET / GEOM_PER_A / GEOM_PER_B (and every
    other boundary entity coplanar with them, i.e. the O-grid footprint),
  - hub / shroud: boundary entities on the radius of GEOM_HUB / GEOM_SHROUD
    (radii read from those entities, never set),
  - blade: every remaining boundary entity (the wall of the dtOO O-grid),
  - cutout: the faces between a dtOO O-grid cell (tet_prep_v5.OGRID_GEOM) and
    any other cell = the exact outer surface of the removed O-grid (label 7).
The blade is additionally available as the dense dtOO sample S(u, v) from
blade/machine_scratch/machine_blade.npz when present.

Feature curves are the 1D entities of the .msh (no dihedral detection).

Projection: closest point on the facet set (exact point-triangle distance over
the k nearest facet centroids).  run as a script it projects the boundary
vertices of blocks.vtk onto the surface their tet.vtk label names and reports
the deviation, plus the self-residual of the reference nodes (gate < 1e-6).
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

import tet_prep_v5 as tp5                                             # noqa: E402

VOL_FACES = {
    "hexahedron": [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)],
    "wedge": [(0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)],
    "tetra": [(0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3)],
    "pyramid": [(0, 3, 2, 1), (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)],
}
# tet.vtk surface label (tet_prep_v5.S_*) -> reference surface
LABEL_SURF = {tp5.S_INLET: "inlet", tp5.S_OUTLET: "outlet", tp5.S_PER_A: "periodic_A",
              tp5.S_PER_B: "periodic_B", tp5.S_HUB: "hub", tp5.S_SHROUD: "shroud",
              tp5.S_OGRID_IFACE: "cutout"}


def closest_on_tris(x, A, B, C):
    """Closest points on triangles (A, B, C) to x, vectorised over triangles (Ericson)."""
    ab, ac, ap = B - A, C - A, x - A
    d1, d2 = (ab * ap).sum(1), (ac * ap).sum(1)
    bp = x - B
    d3, d4 = (ab * bp).sum(1), (ac * bp).sum(1)
    cp = x - C
    d5, d6 = (ab * cp).sum(1), (ac * cp).sum(1)
    va = d3 * d6 - d5 * d4
    vb = d5 * d2 - d1 * d6
    vc = d1 * d4 - d3 * d2
    den = va + vb + vc
    den = np.where(np.abs(den) < 1e-300, 1e-300, den)
    v, w = vb / den, vc / den
    out = A + ab * v[:, None] + ac * w[:, None]
    # regions (vertex / edge) -- evaluated last-wins in reverse priority
    def put(mask, val):
        out[mask] = val[mask]
    e_bc = (va <= 0) & ((d4 - d3) >= 0) & ((d5 - d6) >= 0)
    t = (d4 - d3) / np.where((d4 - d3) + (d5 - d6) == 0, 1, (d4 - d3) + (d5 - d6))
    put(e_bc, B + (C - B) * t[:, None])
    e_ac = (vb <= 0) & (d2 >= 0) & (d6 <= 0)
    t = d2 / np.where(d2 - d6 == 0, 1, d2 - d6)
    put(e_ac, A + ac * t[:, None])
    e_ab = (vc <= 0) & (d1 >= 0) & (d3 <= 0)
    t = d1 / np.where(d1 - d3 == 0, 1, d1 - d3)
    put(e_ab, A + ab * t[:, None])
    put((d6 >= 0) & (d5 <= d6), C)
    put((d3 >= 0) & (d4 <= d3), B)
    put((d1 <= 0) & (d2 <= 0), A)
    return out


class Surface:
    def __init__(self, name, points, tris):
        self.name, self.P, self.T = name, points, tris
        self.A, self.B, self.C = points[tris[:, 0]], points[tris[:, 1]], points[tris[:, 2]]
        self.kd = cKDTree((self.A + self.B + self.C) / 3.0)

    def project(self, X, k=24):
        X = np.atleast_2d(np.asarray(X, float))
        _, nn = self.kd.query(X, k=min(k, len(self.T)))
        nn = np.atleast_2d(nn)
        out = np.empty_like(X)
        for i, x in enumerate(X):
            c = nn[i]
            Q = closest_on_tris(x[None, :], self.A[c], self.B[c], self.C[c])
            out[i] = Q[np.argmin(((Q - x) ** 2).sum(1))]
        return out


class Curve:
    def __init__(self, ent, points):
        self.ent, self.P = ent, points

    def project(self, X):
        X = np.atleast_2d(np.asarray(X, float))
        A, B = self.P[:-1], self.P[1:]
        d = B - A
        L2 = np.maximum((d * d).sum(1), 1e-300)
        out = np.empty_like(X)
        for i, x in enumerate(X):
            t = np.clip(((x - A) * d).sum(1) / L2, 0, 1)
            Q = A + d * t[:, None]
            out[i] = Q[np.argmin(((Q - x) ** 2).sum(1))]
        return out


def _order_quad(P, ids):
    """Cyclic order of 4 (roughly planar) quad nodes."""
    ids = np.asarray(list(ids), int)
    X = P[ids]
    c = X.mean(0)
    n = np.cross(X[1] - X[0], X[2] - X[0])
    if np.linalg.norm(n) < 1e-14:
        n = np.cross(X[2] - X[0], X[3] - X[0])
    n /= np.linalg.norm(n)
    u = X[0] - c
    u /= np.linalg.norm(u)
    w = np.cross(n, u)
    ang = np.arctan2((X - c) @ w, (X - c) @ u)
    return ids[np.argsort(ang)]


def _chain(edges):
    """Order 2-node line elements of one entity into a polyline (node ids)."""
    adj = {}
    for a, b in edges:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    ends = [n for n, v in adj.items() if len(v) == 1]
    start = ends[0] if ends else edges[0][0]
    path, prev = [start], None
    while True:
        nxt = [n for n in adj[path[-1]] if n != prev]
        if not nxt or (nxt[0] == start):
            break
        prev = path[-1]
        path.append(nxt[0])
    return path


class Geometry:
    def __init__(self, run):
        import meshio
        run = Path(run)
        m = meshio.read(str(run / "case" / "mesh.msh"))
        P = np.asarray(m.points, float)
        self.points = P
        geo = m.cell_data["gmsh:geometrical"]
        count, side = {}, {}
        for c, g in zip(m.cells, geo):
            if c.type not in VOL_FACES:
                continue
            og = np.isin(np.asarray(g).ravel(), list(tp5.OGRID_GEOM))
            for f in VOL_FACES[c.type]:
                for e, o in zip(c.data[:, list(f)], og):
                    k = frozenset(e.tolist())
                    count[k] = count.get(k, 0) + 1
                    side.setdefault(k, set()).add(bool(o))
        # the O-grid cut = faces between a dtOO O-grid cell and a non-O-grid cell
        cut = [list(k) for k, v in side.items() if len(v) == 2]
        # boundary faces per 2D entity (tris + quads, quads split)
        ent_tris, self.curves = {}, {}
        for c, g in zip(m.cells, geo):
            g = np.asarray(g).ravel()
            if c.type in ("triangle", "quad"):
                for e, t in zip(c.data, g):
                    if count.get(frozenset(e.tolist()), 0) != 1:
                        continue                     # internal interface of the dtOO mesh
                    tl = ent_tris.setdefault(int(t), [])
                    if len(e) == 3:
                        tl.append(e)
                    else:
                        tl += [e[[0, 1, 2]], e[[0, 2, 3]]]
            elif c.type == "line":
                for t in np.unique(g):
                    path = _chain(c.data[g == t].tolist())
                    self.curves[int(t)] = Curve(int(t), P[path])
        r = np.hypot(P[:, 0], P[:, 1])
        tagged = {tp5.GEOM_INLET: "inlet", tp5.GEOM_OUTLET: "outlet",
                  tp5.GEOM_PER_A: "periodic_A", tp5.GEOM_PER_B: "periodic_B"}
        self.r_hub = float(np.median(r[np.unique(np.ravel(ent_tris[tp5.GEOM_HUB]))]))
        self.r_shroud = float(np.median(r[np.unique(np.ravel(ent_tris[tp5.GEOM_SHROUD]))]))
        rtol = 1e-6 * self.r_shroud

        # a plane / periodic surface of a tagged entity: fit plane (inlet/outlet) or
        # match by shared boundary nodes (periodic). Footprint entities of the O-grid
        # share a boundary curve with their tagged neighbour and lie on the same
        # surface -> attach an untagged entity to the tagged entity it shares most
        # nodes with, if all its nodes are within rtol of that entity's surface.
        groups = {n: [] for n in ("inlet", "outlet", "periodic_A", "periodic_B", "hub", "shroud", "blade")}
        for t, n in tagged.items():
            groups[n].append(t)
        groups["hub"].append(tp5.GEOM_HUB)
        groups["shroud"].append(tp5.GEOM_SHROUD)
        base = {n: Surface(n, P, np.asarray(np.vstack([ent_tris[t] for t in ts]), int))
                for n, ts in groups.items() if ts}
        self.entity_surface = {}
        for t, tl in ent_tris.items():
            ids = np.unique(np.ravel(tl))
            rr = r[ids]
            if t in tagged or t in (tp5.GEOM_HUB, tp5.GEOM_SHROUD):
                continue
            if np.all(np.abs(rr - self.r_hub) < rtol):
                groups["hub"].append(t)
                continue
            if np.all(np.abs(rr - self.r_shroud) < rtol):
                groups["shroud"].append(t)
                continue
            hit = None
            for n in ("inlet", "outlet", "periodic_A", "periodic_B"):
                d = np.linalg.norm(base[n].project(P[ids[:: max(1, len(ids) // 40)]]) - P[ids[:: max(1, len(ids) // 40)]], axis=1)
                if d.max() < 1e-6:
                    hit = n
                    break
            groups[hit or "blade"].append(t)
        self.groups = groups
        self.surfaces = {}
        for n, ts in groups.items():
            T = np.asarray(np.vstack([ent_tris[t] for t in ts]), int)
            self.surfaces[n] = Surface(n, P, T)
            for t in ts:
                self.entity_surface[t] = n
        cut_tris = []
        self.cut_quads = []
        for k in cut:
            if len(k) == 3:
                cut_tris.append(k)
            else:                                   # frozenset loses the quad order: rebuild it
                q = _order_quad(P, k)
                self.cut_quads.append(q)
                cut_tris += [q[[0, 1, 2]], q[[0, 2, 3]]]
        self.surfaces["cutout"] = Surface("cutout", P, np.asarray(cut_tris, int))
        npz = run / "blade" / "machine_scratch" / "machine_blade.npz"
        self.blade_S = np.load(npz, allow_pickle=True)["S"] if npz.exists() else None
        # hub / shroud from the analytic dtOO surfaces when exported (extract_blade.py):
        # the msh facets sit up to h^2 / 8r inside the cylinders (0.35 mm here), more
        # than a third of the first BL cell. Only the part of the full revolution near
        # the msh surface is kept.
        self.analytic = []
        wnpz = run / "blade" / "machine_scratch" / "machine_walls.npz"
        if wnpz.exists():
            W = np.load(wnpz)
            for key in ("hub", "shroud"):
                S = W[f"{key}_S"]
                nu, nv = S.shape[:2]
                I, J = np.meshgrid(np.arange(nu - 1), np.arange(nv - 1), indexing="ij")
                a0 = (I * nv + J).ravel()
                T = np.r_[np.c_[a0, a0 + nv, a0 + nv + 1], np.c_[a0, a0 + nv + 1, a0 + 1]]
                X = S.reshape(-1, 3)
                ref = self.surfaces[key]
                near = cKDTree(ref.P[np.unique(ref.T)]).query(X[T].mean(1))[0] < 0.1
                self.surfaces[key] = Surface(key, X, T[near])
                self.analytic.append(key)

    def project(self, name, X):
        return self.surfaces[name].project(X)

    def cutout_grid(self):
        """The cutout as the structured grid it is in the dtOO O-grid: S[i, k] with i
        around the blade (closed: S[-1] == S[0]) and k from hub to shroud. Rows are
        marched from the node loop on the hub across the quad faces: every node has
        two neighbours in its row, one below and one above. None if the cutout is not
        such a quad grid (then the caller has no (s, t) parameters)."""
        if not self.cut_quads:
            return None
        P = self.points
        adj = {}
        for q in self.cut_quads:
            for a_, b_ in zip(q, np.roll(q, -1)):
                adj.setdefault(int(a_), set()).add(int(b_))
                adj.setdefault(int(b_), set()).add(int(a_))
        r = np.hypot(P[:, 0], P[:, 1])
        rtol = 1e-6 * self.r_shroud
        hub = {v for v in adj if abs(r[v] - self.r_hub) < rtol}
        if not hub:
            return None
        row = _chain([(a_, b_) for a_ in hub for b_ in adj[a_] if b_ in hub and a_ < b_])
        rows, prev = [row], set()
        while not all(abs(r[v] - self.r_shroud) < rtol for v in rows[-1]):
            cur = set(rows[-1])
            nxt = []
            for v in rows[-1]:
                up = [w for w in adj[v] if w not in cur and w not in prev]
                if len(up) != 1:
                    return None
                nxt.append(up[0])
            prev = cur
            rows.append(nxt)
            if len(rows) > len(adj):
                return None
        G = P[np.array(rows)]                                  # (K, N, 3)
        return np.concatenate([G, G[:, :1]], axis=1).transpose(1, 0, 2)

    def theta_map(self, name):
        """theta(r, z) of a periodic surface (a graph over the meridian plane:
        rotational periodicity is circumferential), Clough-Tocher over its msh
        nodes; theta unwrapped about the surface's mean."""
        if not hasattr(self, "_theta"):
            self._theta = {}
        if name not in self._theta:
            from scipy.interpolate import CloughTocher2DInterpolator
            S = self.surfaces[name]
            ids = np.unique(S.T)
            X = S.P[ids]
            th = np.arctan2(X[:, 1], X[:, 0])
            th0 = float(np.angle(np.exp(1j * th).mean()))
            th = np.angle(np.exp(1j * (th - th0)))
            f = CloughTocher2DInterpolator(np.c_[np.hypot(X[:, 0], X[:, 1]), X[:, 2]], th)
            self._theta[name] = (f, th0)
        return self._theta[name]

    def project_circumferential(self, name, X):
        """Move X along the circle about the machine axis onto a periodic surface:
        r and z stay, only theta changes (normal projection onto a twisted periodic
        surface also moves r and z, which sheared the BL in B5). Outside the
        surface's (r, z) hull: normal projection."""
        X = np.atleast_2d(np.asarray(X, float))
        f, th0 = self.theta_map(name)
        r = np.hypot(X[:, 0], X[:, 1])
        th = f(np.c_[r, X[:, 2]]) + th0
        out = np.c_[r * np.cos(th), r * np.sin(th), X[:, 2]]
        bad = ~np.isfinite(th)
        if bad.any():
            out[bad] = self.project(name, X[bad])
        return out

    def project_any(self, name, X):
        """Circumferential for periodic surfaces, closest point otherwise."""
        if name.startswith("periodic"):
            return self.project_circumferential(name, X)
        return self.project(name, X)

    def project_curve(self, curve, names, X, iters=6):
        """Closest point on a dtOO curve that lies on the surfaces `names`: the msh
        polyline first, then alternating projection onto those surfaces (the
        polyline's chords leave a curved surface between its nodes)."""
        Y = curve.project(X)
        for _ in range(iters):
            for n in names:
                Y = self.project(n, Y)
        return Y


def _block_boundary_vertices(run):
    """blocks.vtk boundary vertices with the tet.vtk label of their faces."""
    import meshio
    m = meshio.read(str(Path(run) / "blocks.vtk"))
    P = np.asarray(m.points, float)
    H = m.cells_dict["hexahedron"]
    cnt = {}
    for h in H:
        for f in VOL_FACES["hexahedron"]:
            q = tuple(h[list(f)])
            k = frozenset(q)
            cnt.setdefault(k, [q, 0])[1] += 1
    bq = [q for q, n in cnt.values() if n == 1]
    t = meshio.read(str(Path(run) / "tet.vtk"))
    tri, lab = [], []
    for c, d in zip(t.cells, t.cell_data["color"]):
        if c.type == "triangle":
            tri.append(c.data)
            lab.append(np.asarray(d).ravel())
    tri, lab = np.vstack(tri), np.concatenate(lab).astype(int)
    kd = cKDTree(t.points[tri].mean(1))
    _, ii = kd.query(np.array([P[list(q)].mean(0) for q in bq]))
    vlab = {}
    for q, i in zip(bq, ii):
        for v in q:
            vlab.setdefault(int(v), set()).add(int(lab[i]))
    return P, vlab


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run")
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    G = Geometry(a.run)
    rep = {"r_hub": G.r_hub, "r_shroud": G.r_shroud, "analytic_surfaces": G.analytic,
           "groups": {n: sorted(ts) for n, ts in G.groups.items()},
           "facets": {n: int(len(s.T)) for n, s in G.surfaces.items()},
           "curves": len(G.curves)}
    # gate: reference nodes reproject onto their own surface
    res = {}
    for n, s in G.surfaces.items():
        ids = np.unique(s.T)[:: max(1, len(np.unique(s.T)) // 400)]
        res[n] = float(np.linalg.norm(s.project(G.points[ids]) - G.points[ids], axis=1).max())
    rep["self_residual_max"] = res
    rep["gate_self_residual"] = bool(max(res.values()) < 1e-6)
    # deviation of the AlgoHex/beam block boundary from the reference
    P, vlab = _block_boundary_vertices(a.run)
    dev = {}
    for v, labs in vlab.items():
        for l in labs:
            n = LABEL_SURF.get(l)
            if n is None:
                continue                    # 5/6 (BL cuts) do not exist with keep-BL
            dev.setdefault(n, []).append(float(np.linalg.norm(G.project(n, P[v])[0] - P[v])))
    rep["blocks_vtk_deviation"] = {n: {"n": len(d), "max": max(d), "p95": float(np.percentile(d, 95)),
                                       "median": float(np.median(d))} for n, d in dev.items()}
    if G.blade_S is not None:
        S = G.blade_S.reshape(-1, 3)[::7]
        dS = np.linalg.norm(G.project("blade", S) - S, axis=1)
        rep["blade_msh_vs_dtoo_S"] = {"max": float(dS.max()), "p95": float(np.percentile(dS, 95))}
    print(json.dumps(rep, indent=1))
    if a.json:
        Path(a.json).write_text(json.dumps(rep, indent=1))
    return 0 if rep["gate_self_residual"] else 1


if __name__ == "__main__":
    sys.exit(main())
