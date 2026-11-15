"""Detectors ``rotation_axis`` and ``sweep``.

rotation_axis
    Every normal line of a surface of revolution meets the axis. With the axis as a line (direction a,
    moment b = a x p), a surface point x with unit normal n lies on such a surface when
        a . (x x n) - b . n = 0,
    which is linear in (a, b). The axis is the smallest right singular vector of the stacked rows
    [x x n, -n], weighted by triangle area (Pottmann & Randrup's line-geometric fit). The fit is repeated
    on the surfaces that agree with it (residual below tol), so a blade or a cut plane does not pull it.
    Per surface the RMS distance of its normal lines to the axis classifies it:
        revolution    small residual (cylinders, cones, tori, planes normal to the axis)
        sector_plane  a plane that contains the axis (periodic / symmetry cut of a sector)
        other         anything else (blades, struts, non-rotational walls)

sweep
    r-z sections of the boundary triangulation in half planes phi = const around the axis. The domain
    is sweepable when every section has the same topology (number of loops, which surfaces they cross);
    how much the section shape changes over phi is reported as the Hausdorff distance to the middle
    section, relative to the domain size.

Both work on the triangulated boundary of any source (.msh, CAD, STL); nothing depends on surface
labels or numbers.
"""

from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

from hexdp.contracts import Feature, ParamSpec


# ------------------------------------------------------------------------------------------- helpers
def _triangles(geo):
    A, B, C, owner = [], [], [], []
    for t, s in geo.surfaces.items():
        if not len(s.tris):
            continue
        A.append(s.nodes[s.tris[:, 0]]); B.append(s.nodes[s.tris[:, 1]]); C.append(s.nodes[s.tris[:, 2]])
        owner.append(np.full(len(s.tris), t))
    return np.vstack(A), np.vstack(B), np.vstack(C), np.concatenate(owner)


def fit_axis(X, N, w):
    """Axis (point p, unit direction a) from surface points X with unit normals N and weights w, plus the
    distance of every normal line to the axis.

    The rows [x x n, -n] have a second null vector (a = 0, b = axis direction: the 'line at infinity')
    whenever no surface normal is parallel to the axis, e.g. a cylinder mantle alone. Of the two smallest
    singular vectors the combination with the largest direction part a is taken, which is the axis."""
    M = np.cross(X, N)
    R = np.hstack([M, -N]) * np.sqrt(w)[:, None]        # callers centre X first (conditioning)
    _, sv, Vt = np.linalg.svd(R, full_matrices=False)
    V = Vt[-1:] if sv[-2] > 10 * sv[-1] + 1e-12 * sv[0] else Vt[-2:]
    G = V[:, :3] @ V[:, :3].T                           # maximise |a| over the (near) null space
    c = np.linalg.eigh(G)[1][:, -1]
    v = c @ V
    a, b = v[:3], v[3:]
    na = np.linalg.norm(a)
    if na < 1e-9:
        return None
    a, b = a / na, b / na
    p = np.cross(b, a)                      # the axis point closest to the origin: b = a x p  =>  p = b x a
    return p, a, line_distances(X, N, p, a)


def line_distances(X, N, p, a):
    """Distance between each normal line (X, N) and the axis (p, a); parallel lines -> their distance."""
    c = np.cross(N, a)
    nc = np.linalg.norm(c, axis=1)
    d = X - p
    with np.errstate(invalid="ignore", divide="ignore"):
        skew = np.abs((d * c).sum(1)) / nc
    par = np.linalg.norm(np.cross(d, a), axis=1)             # normal parallel to the axis
    # a normal (nearly) parallel to the axis belongs to a plane normal to the axis: consistent with it
    return np.where(nc > 1e-3, skew, 0.0 * par)


def frame(p, a):
    """Orthonormal frame (e1, e2, a) with e1, e2 normal to the axis."""
    ref = np.array([1.0, 0, 0]) if abs(a[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = np.cross(a, ref); e1 /= np.linalg.norm(e1)
    return e1, np.cross(a, e1), a


def to_cyl(P, p, a):
    e1, e2, e3 = frame(p, a)
    d = P - p
    x, y, z = d @ e1, d @ e2, d @ e3
    return np.hypot(x, y), np.arctan2(y, x), z


# ------------------------------------------------------------------------------------------- detectors
class RotationAxisDetector:
    """Finds the rotation axis shared by the most curved boundary area and sorts the surfaces by it."""

    name = "rotation_axis"
    params = {
        "tol": ParamSpec(0.3, "max RMS normal-line distance of a revolution surface, relative to the facet size "
                              "(a flat facet of a curved surface tilts its normal by up to ~0.15)", lower=0),
        "min_curved_area": ParamSpec(0.005, "min fraction of the boundary area on curved revolution surfaces",
                                     lower=0, upper=1),
    }

    def detect(self, geo, params, ctx):
        A, B, C, owner = _triangles(geo)
        X = (A + B + C) / 3
        n = np.cross(B - A, C - A)
        area = 0.5 * np.linalg.norm(n, axis=1)
        N = n / np.maximum(2 * area, 1e-300)[:, None]
        h = np.maximum.reduce([np.linalg.norm(B - A, axis=1), np.linalg.norm(C - B, axis=1),
                               np.linalg.norm(A - C, axis=1)])
        lo, hi = X.min(0), X.max(0)
        scale = float(np.linalg.norm(hi - lo))
        shift = 0.5 * (lo + hi)                        # fit in centred coordinates (conditioning)
        Xc = X - shift
        tags = [int(t) for t in np.unique(owner)]
        idx = {t: np.where(owner == t)[0] for t in tags}
        tol = params["tol"]
        total = float(area.sum())

        def curved(t):
            m = idx[t]
            return len(m) >= 6 and np.linalg.norm(N[m] - N[m].mean(0), axis=1).max() > 1e-2

        def rms_all(p, a):
            d = line_distances(Xc, N, p, a) / np.maximum(h, 1e-300)
            return {t: float(np.sqrt(np.average(d[idx[t]] ** 2, weights=area[idx[t]]))) for t in tags}

        # candidates: the axis of every curved surface that is a surface of revolution on its own
        best = None
        for t in filter(curved, tags):
            f = fit_axis(Xc[idx[t]], N[idx[t]], area[idx[t]])
            if f is None:
                continue
            r = rms_all(f[0], f[1])
            if r[t] > tol:
                continue
            support = sum(area[idx[u]].sum() for u in tags if curved(u) and r[u] <= tol)
            if best is None or support > best[0]:
                best = (support, f[0], f[1])
        if best is None or best[0] / total < params["min_curved_area"]:
            return []
        # refit on all supporting surfaces (curved and planes normal to the axis), classify everything
        r = rms_all(best[1], best[2])
        sup = [u for u in tags if r[u] <= tol]
        sel = np.concatenate([idx[u] for u in sup])
        f = fit_axis(Xc[sel], N[sel], area[sel])
        p, a = (f[0], f[1]) if f is not None else (best[1], best[2])
        rms = rms_all(p, a)
        good = {u for u in tags if rms[u] <= tol}
        p = p + shift                                   # back to model coordinates
        p = p - np.dot(p, a) * a                        # the axis point closest to the origin
        a = a if a[np.argmax(np.abs(a))] > 0 else -a    # canonical sign
        frac = float(sum(area[idx[u]].sum() for u in good)) / total
        curved_frac = float(sum(area[idx[u]].sum() for u in good if curved(u))) / total
        # planes containing the axis: flat surfaces whose plane holds the axis
        sector, other = [], []
        for t in tags:
            if t in good:
                continue
            nn = N[idx[t]]
            flat = np.linalg.norm(nn - nn.mean(0), axis=1).max() < 1e-3
            nmean = nn.mean(0) / max(np.linalg.norm(nn.mean(0)), 1e-300)
            through = abs(np.dot(X[idx[t]].mean(0) - p, nmean)) < 1e-4 * scale and abs(np.dot(nmean, a)) < 1e-3
            (sector if flat and through else other).append(t)
        rr, th, z = to_cyl(np.vstack([geo.surfaces[t].nodes for t in tags]), p, a)
        span = _angular_span(th[rr > 1e-9 * scale])
        reasons = [f"axis direction {np.round(a, 4).tolist()} through {np.round(p, 4).tolist()}",
                   f"{len(good)} revolution surfaces ({100 * frac:.0f} % of the boundary area, "
                   f"{100 * curved_frac:.0f} % curved)",
                   f"{len(sector)} planes containing the axis", f"angular span {np.degrees(span):.1f} deg"]
        if other:
            reasons.append(f"{len(other)} non-rotational surfaces ({100 * (1 - frac) - 0:.0f} % of the area)")
        return [Feature("rotation_axis", {
            "point": p.tolist(), "direction": a.tolist(), "angular_span_deg": float(np.degrees(span)),
            "revolution_surfaces": sorted(good), "sector_planes": sector, "other_surfaces": other,
            "rms_rel_facet": rms, "area_fraction": frac, "curved_area_fraction": curved_frac},
            entities=tuple(sorted(good)), confidence=round(frac, 3), reasons=reasons)]


def _angular_span(theta):
    """Smallest angular interval holding all angles (radians, in [0, 2 pi])."""
    if len(theta) == 0:
        return 0.0
    t = np.sort(np.mod(theta, 2 * np.pi))
    gaps = np.diff(np.r_[t, t[0] + 2 * np.pi])
    return float(2 * np.pi - gaps.max())


def section(A, B, C, owner, p, a, phi, eps=0.0):
    """Intersection of the triangles with the half plane phi (around axis p, a): segments in (r, z) and
    the surface tag of each segment. Segments shorter than eps or lying on the axis (r < eps) are dropped:
    every half plane passes through the axis, so a surface that touches the axis (a sector plane) would
    otherwise show up as a point-sized cut in every section."""
    e1, e2, e3 = frame(p, a)
    er = np.cos(phi) * e1 + np.sin(phi) * e2
    nrm = np.cross(e3, er)                                  # plane normal (tangential direction)
    sa, sb, sc = ((V - p) @ nrm for V in (A, B, C))
    S = np.stack([sa, sb, sc], 1)
    cross = (S.min(1) < 0) & (S.max(1) > 0)
    segs, tags = [], []
    V = np.stack([A, B, C], 1)
    for i in np.where(cross)[0]:
        pts = []
        for j, k in ((0, 1), (1, 2), (2, 0)):
            if (S[i, j] < 0) != (S[i, k] < 0):
                t = S[i, j] / (S[i, j] - S[i, k])
                pts.append(V[i, j] + t * (V[i, k] - V[i, j]))
        if len(pts) == 2:
            q = np.array(pts) - p
            rz = np.c_[q @ er, q @ e3]
            if rz[:, 0].min() >= 0 and rz[:, 0].max() > eps and np.linalg.norm(rz[1] - rz[0]) > eps:
                segs.append(rz)
                tags.append(owner[i])
    return segs, tags


def loops(segs, tol):
    """Chain segments into polylines; returns the number of chains and of closed ones."""
    if not segs:
        return 0, 0
    E = np.array(segs)                                      # (m, 2, 2)
    P = E.reshape(-1, 2)
    tree = cKDTree(P)
    groups = tree.query_ball_point(P, tol)
    rep = np.arange(len(P))
    for i, g in enumerate(groups):
        rep[i] = min(g)
    ends = rep.reshape(-1, 2)
    parent = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    deg = {}
    for u, v in ends:
        deg[u] = deg.get(u, 0) + 1; deg[v] = deg.get(v, 0) + 1
        parent[find(u)] = find(v)
    comps = {}
    for u, v in ends:
        comps.setdefault(find(u), set()).update((u, v))
    closed = sum(all(deg[n] == 2 for n in c) for c in comps.values())
    return len(comps), closed


class SweepDetector:
    """Checks whether r-z sections around the rotation axis keep their topology over the angular span."""

    name = "sweep"
    params = {
        "sections": ParamSpec(9, "number of test sections over the span", lower=3),
        "margin_deg": ParamSpec(1.0, "keep the test sections this far from the ends of the span", unit="deg", lower=0),
        "hausdorff_scale": ParamSpec(0.2, "section shape change (Hausdorff / size) at which the confidence "
                                          "has dropped to one half; the decision is the topology", lower=1e-6),
    }

    def detect(self, geo, params, ctx):
        rd = RotationAxisDetector()
        axis = rd.detect(geo, {k: v.default for k, v in rd.params.items()}, ctx)
        # without an axis there is nothing to sweep around
        if not axis:
            return []
        ax = axis[0].data
        p, a = np.array(ax["point"]), np.array(ax["direction"])
        A, B, C, owner = _triangles(geo)
        Pall = np.vstack([A, B, C])
        r, th, _ = to_cyl(Pall, p, a)
        scale = float(np.linalg.norm(Pall.max(0) - Pall.min(0)))
        span = np.radians(ax["angular_span_deg"])
        t0 = _span_start(th[r > 1e-9 * scale], span)
        m = np.radians(params["margin_deg"])
        if span >= 2 * np.pi - 1e-6:
            phis = np.linspace(0, 2 * np.pi, params["sections"], endpoint=False)
        else:
            phis = np.linspace(t0 + m, t0 + span - m, params["sections"])
        phis = phis + 1.234567e-6                       # never exactly through mesh vertices
        sigs, shapes = [], []
        for ph in phis:
            segs, tags = section(A, B, C, owner, p, a, ph, eps=1e-6 * scale)
            n, closed = loops(segs, 1e-5 * scale)
            cut = set(int(t) for t in tags) - set(ax["sector_planes"])      # the sector planes are the ends
            sigs.append((n, closed, tuple(sorted(cut))))
            shapes.append(np.vstack([s.mean(0)[None] for s in segs]) if segs else np.zeros((0, 2)))
        mid = shapes[len(shapes) // 2]
        haus = []
        for S in shapes:
            if not len(S) or not len(mid):
                haus.append(np.inf); continue
            d1 = cKDTree(mid).query(S)[0].max(); d2 = cKDTree(S).query(mid)[0].max()
            haus.append(max(d1, d2) / scale)
        # topology = number of chains + surfaces cut; "closed" is fragile where a section grazes a mesh
        # edge (a node of degree 4), so it is reported but not compared
        same = all(s[0] == sigs[0][0] for s in sigs)
        same_faces = all(s[2] == sigs[0][2] for s in sigs)
        hmax = float(np.max(haus))
        reasons = [f"{len(phis)} sections over {ax['angular_span_deg']:.1f} deg",
                   f"loops per section: {sorted({s[0] for s in sigs})} (closed {sorted({s[1] for s in sigs})})",
                   "the same surfaces in every section" if same_faces else "the surfaces cut change with phi",
                   f"max section change (Hausdorff / size) {hmax:.4f}"]
        ok = same and same_faces
        data = {"axis_point": p.tolist(), "axis_direction": a.tolist(), "phi_deg": np.degrees(phis).round(3).tolist(),
                "angular_span_deg": ax["angular_span_deg"], "loops": [s[0] for s in sigs],
                "surfaces": [list(s[2]) for s in sigs], "hausdorff_rel": [float(h) for h in haus], "sweepable": ok}
        if not ok:
            return [Feature("sweep_rejected", data, confidence=0.0, reasons=reasons)]
        conf = round(1.0 / (1.0 + hmax / params["hausdorff_scale"]), 3)
        return [Feature("sweep", data, confidence=conf, reasons=reasons)]


def _span_start(theta, span):
    """Start angle of the smallest interval of length span that holds all angles."""
    t = np.sort(np.mod(theta, 2 * np.pi))
    gaps = np.diff(np.r_[t, t[0] + 2 * np.pi])
    i = int(np.argmax(gaps))
    return float(t[(i + 1) % len(t)])
