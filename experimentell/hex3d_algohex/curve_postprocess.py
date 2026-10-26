"""Curve post-processing of a block sample: continuity of the block edges and exact feature curves.

Works on the neutral sample (`export_sample.build` / `export_sample.load`) and needs numpy + scipy only
(decision C). It never re-runs AlgoHex; it is the cheap side of the seam. Everything is decided from
the block topology and the patch structure of the input surface, with no label numbers and no
geometry-specific rules, so it carries over to other geometry classes:

  feature edge   a boundary block edge whose two boundary quads lie on different surface patches
  corner         a vertex on >= 3 patches, at a seam junction, or where the seam curve has a real kink
                 (turn > KINK_DEG that does not shrink when the window is halved; a smooth hairpin turns
                 ~linearly with the window, a kink keeps its angle)
  pairs          two edges at a vertex that continue one curve (G1), mutually unique, by priority
                   feature  two feature edges on the same patch pair, no shared boundary quad, not at a corner
                   2d       two non-feature boundary edges in the same patch, no shared boundary quad
                   3d       two interior edges that share no block face
                 vertices with several candidates are singular: no forced pairing, listed as ambiguous
  feature arcs   each feature edge follows its seam curve of the input triangulation exactly. The seam
                 polylines are interpolated by C2 splines (the raw polyline turns by up to ~40 deg per
                 node at tight hairpins); vertices inside a feature line are moved onto their seam first
                 (two passes: assign, move, assign again)
  G1 cubic       one cubic per feature edge for consumers that need the Bezier form: end tangents = seam
                 tangents (so consecutive feature edges are G1 by construction), arm lengths by least
                 squares against the arc

The result is a dict of `post_*` arrays (see `postprocess`) that rides along in the sample (export
`--curve-postprocess`) or in a side file (`<sample>_curves.npz`, CLI below). The original `vertices` /
`edge_ctrl` stay untouched.

Run on existing samples:
    PY experimentell/hex3d_algohex/curve_postprocess.py runs/*/sample.npz [--out-dir DIR] [--inplace]
"""

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.interpolate import splev, splprep
from scipy.spatial import cKDTree

# VTK hex faces as corner slots (same as export_sample.FACES)
FACES = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
KINK_DEG = 30.0          # seam turn above this (and not shrinking with the window) = corner
VERTEX_TOL = 6e-3        # vertex on a patch: closer than this fraction of the surface bbox diagonal
EDGE_TOL = 0.2           # feature edge on a seam curve: both ends within this fraction of the chord
ARC_POINTS = 65          # samples per feature arc
SPLINE_REFINE = 16       # spline samples per seam polyline segment
ROLE = {"regular": 0, "feature": 1, "corner": 2, "singular": 3}
KIND = {"feature": 0, "2d": 1, "3d": 2}


# --------------------------------------------------------------------------
# surface: point-triangle distance, patch membership, seam curves
# --------------------------------------------------------------------------

def tri_dist(X, A, B, C, chunk=16):
    """Distance of points X [N,3] to the triangle set (A,B,C) [T,3] each -> [N]
    (closest point on triangle, Ericson, Real-Time Collision Detection 5.1.5)."""
    out = np.full(len(X), np.inf)
    for i0 in range(0, len(X), chunk):
        p = X[i0:i0 + chunk, None]
        ab, ac, ap = B - A, C - A, p - A
        d1, d2 = (ab * ap).sum(-1), (ac * ap).sum(-1)
        bp = p - B; d3, d4 = (ab * bp).sum(-1), (ac * bp).sum(-1)
        cp = p - C; d5, d6 = (ab * cp).sum(-1), (ac * cp).sum(-1)
        va, vb, vc = d3 * d6 - d5 * d4, d5 * d2 - d1 * d6, d1 * d4 - d3 * d2
        den = np.where(np.abs(va + vb + vc) < 1e-30, 1e-30, va + vb + vc)
        v, w = vb / den, vc / den
        q = A + ab * v[..., None] + ac * w[..., None]
        e_bc = (va <= 0) & (d4 - d3 >= 0) & (d5 - d6 >= 0)
        t = np.clip((d4 - d3) / np.where((d4 - d3) + (d5 - d6) == 0, 1e-30, (d4 - d3) + (d5 - d6)), 0, 1)
        q = np.where(e_bc[..., None], B + (C - B) * t[..., None], q)
        e_ac = (vb <= 0) & (d2 >= 0) & (d6 <= 0)
        t = np.clip(d2 / np.where(d2 - d6 == 0, 1e-30, d2 - d6), 0, 1)
        q = np.where(e_ac[..., None], A + ac * t[..., None], q)
        q = np.where(((d6 >= 0) & (d5 <= d6))[..., None], C, q)
        e_ab = (vc <= 0) & (d1 >= 0) & (d3 <= 0)
        t = np.clip(d1 / np.where(d1 - d3 == 0, 1e-30, d1 - d3), 0, 1)
        q = np.where(e_ab[..., None], A + ab * t[..., None], q)
        q = np.where(((d3 >= 0) & (d4 <= d3))[..., None], B, q)
        q = np.where(((d1 <= 0) & (d2 <= 0))[..., None], A, q)
        out[i0:i0 + chunk] = np.linalg.norm(p - q, axis=-1).min(1)
    return out


def patch_dists(X, pts, tris, labels):
    """[N, n_labels] distance of every point to every patch, and the label list."""
    labs = [int(l) for l in np.unique(labels)]
    D = np.stack([tri_dist(X, *(pts[tris[labels == l][:, k]] for k in range(3))) for l in labs], 1)
    return D, labs


def seam_polylines(pts, tris, labels):
    """Label boundaries of the surface triangulation as polylines: list of (points [n,3], closed, patch labels).
    A boundary edge separates triangles of different labels; nodes of degree != 2 split chains."""
    emap = defaultdict(set)
    for t, (a, b, c) in enumerate(tris):
        for u, v in ((a, b), (b, c), (c, a)):
            emap[(min(u, v), max(u, v))].add(int(labels[t]))
    bnd = {e for e, ls in emap.items() if len(ls) > 1}
    deg, adj = Counter(), defaultdict(list)
    for u, v in bnd:
        deg[u] += 1; deg[v] += 1; adj[u].append(v); adj[v].append(u)
    used, chains = set(), []

    def walk(a, b):
        ch = [a, b]; used.add((min(a, b), max(a, b)))
        prev, cur = a, b
        while deg[cur] == 2:
            nxt = [w for w in adj[cur] if w != prev]
            if not nxt:
                break
            w = nxt[0]; e = (min(cur, w), max(cur, w))
            if e in used:
                break
            used.add(e); ch.append(w); prev, cur = cur, w
        return ch

    for n in sorted(x for x in deg if deg[x] != 2):
        for w in adj[n]:
            if (min(n, w), max(n, w)) not in used:
                chains.append(walk(n, w))
    for e in sorted(bnd):
        if e not in used:
            chains.append(walk(*e))
    out = []
    for ch in chains:
        if len(ch) < 2:
            continue
        Q = pts[np.asarray(ch, int)]
        closed = len(ch) > 2 and ch[0] == ch[-1]
        labs = tuple(sorted(set().union(*(emap[(min(u, v), max(u, v))] for u, v in zip(ch[:-1], ch[1:])))))
        out.append((Q, closed, labs))
    return out


class Seams:
    """Seam polylines interpolated by C2 cubic splines through every node, sampled densely."""

    def __init__(self, polylines, refine=SPLINE_REFINE):
        self.pts, self.arclen, self.closed, self.trees, self.raw, self.labels = [], [], [], [], [], []
        for item in polylines:
            P, closed = item[0], item[1]
            self.labels.append(tuple(item[2]) if len(item) > 2 else ())
            if closed:
                P = P[:-1]
            P = P[np.r_[True, np.linalg.norm(np.diff(P, axis=0), axis=1) > 1e-12]]
            if len(P) >= 4:
                u = np.r_[0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
                if closed:
                    u = np.r_[u, u[-1] + np.linalg.norm(P[0] - P[-1])]
                    P = np.vstack([P, P[:1]])
                tck, _ = splprep(P.T, u=u, s=0, per=int(closed), k=3)
                Q = np.stack(splev(np.linspace(0, u[-1], refine * (len(P) - 1) + 1), tck), 1)
            else:
                Q = P if not closed else np.vstack([P, P[:1]])
            self.raw.append(P)
            self.pts.append(Q); self.closed.append(bool(closed))
            self.arclen.append(np.r_[0, np.cumsum(np.linalg.norm(np.diff(Q, axis=0), axis=1))])
            self.trees.append(cKDTree(Q))
        sp = [np.linalg.norm(np.diff(P, axis=0), axis=1) for P in self.raw if len(P) > 1]
        self.spacing = float(np.median(np.concatenate(sp))) if sp else 0.0

    def __len__(self):
        return len(self.pts)

    def on_curve(self, c, q):
        """(distance, arclength, point) of the nearest sample of curve c for points q [n,3]."""
        d, j = self.trees[c].query(np.asarray(q, float).reshape(-1, 3))
        return d, self.arclen[c][j], self.pts[c][j]

    def nearest(self, q):
        """Nearest curve of one point: (curve, distance)."""
        best = min(((float(self.on_curve(c, q)[0][0]), c) for c in range(len(self))), default=(np.inf, -1))
        return best[1], best[0]

    def sample(self, c, t):
        s, P = self.arclen[c], self.pts[c]
        return np.stack([np.interp(t, s, P[:, k]) for k in range(3)], -1)

    def turn(self, q, h):
        """Turning angle (deg) of the nearest seam over +-h at the point nearest q, and whether q sits at an
        open end (junction) of that curve. (None, None) without seams."""
        c, _ = self.nearest(q)
        if c < 0:
            return None, None
        _, t, _ = self.on_curve(c, q)
        t = float(t[0]); L = float(self.arclen[c][-1])
        if not self.closed[c] and (t < 0.5 * h or t > L - 0.5 * h):
            return None, True
        ta, tb = (t - h) % L if self.closed[c] else max(t - h, 0.0), (t + h) % L if self.closed[c] else min(t + h, L)
        A, P, B = self.sample(c, np.array([ta, t, tb]))
        u, w = P - A, B - P
        cosv = u @ w / max(np.linalg.norm(u) * np.linalg.norm(w), 1e-30)
        return float(np.degrees(np.arccos(np.clip(cosv, -1, 1)))), False

    def tangent(self, c, q):
        """Unit tangent of curve c at the point nearest q (local spline segment; same for every caller at
        the same point, so two edges meeting there get the same line), sign arbitrary."""
        _, t, _ = self.on_curve(c, q)
        s = self.arclen[c]; L = float(s[-1]); t = float(t[0])
        d = 0.5 * float(np.min(np.diff(s)[np.diff(s) > 0]))
        ta, tb = t - d, t + d
        if self.closed[c]:
            ta, tb = ta % L, tb % L
        else:
            ta, tb = max(ta, 0.0), min(tb, L)
        A, B = self.sample(c, np.array([ta, tb]))
        u = B - A
        return u / max(np.linalg.norm(u), 1e-30)

    def edge_curve(self, pa, pb, hint, tol):
        """The curve carrying both endpoints (distance <= tol); among several (junctions) the one nearest
        the edge midpoint `hint`. The midpoint itself may be far: a cubic may cut across a hairpin."""
        best = None
        for c in range(len(self)):
            d = self.on_curve(c, np.stack([pa, pb, hint]))[0]
            if max(d[0], d[1]) <= tol and (best is None or d[2] < best[0]):
                best = (float(d[2]), c)
        return None if best is None else best[1]

    def arc(self, c, pa, pb, hint, n=ARC_POINTS, others=()):
        """Curve c from pa to pb, n points by arclength. On a closed curve two arcs connect pa and pb: the one
        that contains none of the curve's other corners (arclength positions `others`) is the block edge; if
        that does not decide, the one whose midpoint is nearer `hint`. The end points are set to pa / pb."""
        s = self.arclen[c]; L = float(s[-1])
        _, t, _ = self.on_curve(c, np.stack([pa, pb]))
        ta, tb = float(t[0]), float(t[1])
        if self.closed[c]:
            fwd = (tb - ta) % L
            cands = [np.linspace(ta, ta + fwd, n), np.linspace(ta, ta - (L - fwd), n)]
        else:
            cands = [np.linspace(ta, tb, n)]

        def n_inside(q):
            sign = 1.0 if q[-1] >= q[0] else -1.0
            span = abs(q[-1] - q[0])
            cnt = 0
            for o in others:
                d = (o - ta) * sign
                if self.closed[c]:
                    d %= L
                cnt += 1e-6 < d < span - 1e-6
            return cnt
        best = min(cands, key=lambda q: (n_inside(q), np.linalg.norm(self.sample(c, q[n // 2:n // 2 + 1] % L)[0] - hint)))
        P = self.sample(c, best % L if self.closed[c] else best)
        P[0], P[-1] = pa, pb
        return P


# --------------------------------------------------------------------------
# block topology
# --------------------------------------------------------------------------

def bezier(P0, B1, B2, P1, n):
    t = np.linspace(0, 1, n)[:, None]
    return (1 - t) ** 3 * P0 + 3 * (1 - t) ** 2 * t * B1 + 3 * (1 - t) * t ** 2 * B2 + t ** 3 * P1


def undirected(sample):
    """Undirected edges (a < b) with their control points oriented a -> b, and the map from the directed
    sample edges to (undirected id, reversed)."""
    E, ctrl, ids, dmap = [], [], {}, []
    for (a, b), (c1, c2) in zip(sample["edges"], sample["edge_ctrl"]):
        a, b = int(a), int(b)
        key = (min(a, b), max(a, b))
        if key not in ids:
            ids[key] = len(E); E.append(key)
            ctrl.append((c1, c2) if a < b else (c2, c1))
        dmap.append((ids[key], a > b))
    return np.asarray(E, np.int64), np.asarray(ctrl, float), dmap


def faces_of(blocks):
    """Unique block faces as ordered quads, and their number of incident blocks."""
    cnt, quad = Counter(), {}
    for b in blocks:
        for f in FACES:
            q = tuple(int(b[i]) for i in f)
            k = tuple(sorted(q)); cnt[k] += 1; quad.setdefault(k, q)
    return [quad[k] for k in quad], [cnt[k] for k in quad]


def edge_face_sets(E, faces):
    eid = {tuple(e): k for k, e in enumerate(map(tuple, E))}
    out = [set() for _ in E]
    for fi, q in enumerate(faces):
        for i in range(4):
            k = eid.get(tuple(sorted((q[i], q[(i + 1) % 4]))))
            if k is not None:
                out[k].add(fi)
    return out


def mutual_pairs(E, cand, share, ok=None):
    """Pairs of edges at a vertex that share no element and are each other's only such partner.
    `ok(e, f)` restricts the candidates further. Returns (pairs {(v, e, f)}, ambiguous [(v, [(e, f)])])."""
    inc = defaultdict(list)
    for k in cand:
        for v in E[k]:
            inc[int(v)].append(k)
    pairs, amb = set(), []
    for v, es in inc.items():
        part = {e: [f for f in es if f != e and not (share[e] & share[f]) and (ok is None or ok(e, f))] for e in es}
        cands = sorted({tuple(sorted((e, f))) for e, fs in part.items() for f in fs})
        uniq = [(e, f) for e, f in cands if part[e] == [f] and part[f] == [e]]
        pairs |= {(v, e, f) for e, f in uniq}
        rest = [p for p in cands if p not in uniq]
        if rest:
            amb.append((v, rest))
    return pairs, amb


# --------------------------------------------------------------------------
# the post-processing
# --------------------------------------------------------------------------

def g1_cubic(P, T0, T1):
    """Cubic from P[0] to P[-1] with end tangent directions T0 (leaving P[0]) and T1 (arriving at P[-1]);
    arm lengths by least squares against the points P (chord-length parameters), clipped >= 0."""
    P0, P1 = P[0], P[-1]
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
    t = (s / max(s[-1], 1e-30))[:, None]
    b0, b1, b2, b3 = (1 - t) ** 3, 3 * (1 - t) ** 2 * t, 3 * (1 - t) * t ** 2, t ** 3
    base = (b0 + b1) * P0 + (b2 + b3) * P1                     # arms of length 0
    A = np.stack([(b1 * T0).ravel(), (-b2 * T1).ravel()], 1)
    x, *_ = np.linalg.lstsq(A, (P - base).ravel(), rcond=None)
    L = np.linalg.norm(P1 - P0)
    a, b = np.clip(x, 0.05 * L, None)
    return P0 + a * T0, P1 - b * T1


def kink(V, E, c1, c2, v, e, f):
    """Angle (deg) between the tangent arriving at v along edge e and leaving v along f (cubic form)."""
    def out_dir(k):
        a, b = E[k]
        t = c1[k] - V[a] if a == v else c2[k] - V[b]
        return t / max(np.linalg.norm(t), 1e-30)
    return float(np.degrees(np.arccos(np.clip(-out_dir(e) @ out_dir(f), -1, 1))))


def postprocess(sample, verbose=True):
    """Return the post_* arrays for one sample (dict as from export_sample.build / load)."""
    V0 = np.asarray(sample["vertices"], float)
    blocks = np.asarray(sample["blocks"], np.int64)
    pts, tris, lab = (np.asarray(sample[k]) for k in ("surface_points", "surface_tris", "surface_tri_label"))
    E, ctrl, dmap = undirected(sample)
    c1, c2 = ctrl[:, 0].copy(), ctrl[:, 1].copy()

    faces, nown = faces_of(blocks)
    bq = [q for q, n in zip(faces, nown) if n == 1]
    share3, share2 = edge_face_sets(E, faces), edge_face_sets(E, bq)
    bedge = [k for k in range(len(E)) if share2[k]]

    # patches: per vertex from the surface, per boundary quad = the patch all four corners share
    # (nearest patch to the curved quad centre when that is not unique)
    tol = VERTEX_TOL * np.linalg.norm(np.ptp(pts, 0))
    Dv, labs = patch_dists(V0, pts, tris, lab)
    vpatch = [{labs[j] for j in np.flatnonzero(Dv[i] < tol)} for i in range(len(V0))]
    eid = {tuple(e): k for k, e in enumerate(map(tuple, E))}
    qlab = []
    for q in bq:
        common = set.intersection(*(vpatch[v] for v in q))
        if len(common) == 1:
            qlab.append(common.pop()); continue
        mids = [bezier(V0[E[k][0]], c1[k], c2[k], V0[E[k][1]], 3)[1]
                for k in (eid[tuple(sorted((q[i], q[(i + 1) % 4])))] for i in range(4))]
        D, _ = patch_dists(np.mean(mids, 0)[None], pts, tris, lab)
        cand = [j for j, l in enumerate(labs) if l in common] or list(range(len(labs)))
        qlab.append(labs[min(cand, key=lambda j: D[0, j])])
    epat = [{qlab[q] for q in share2[k]} for k in range(len(E))]
    feat = {k for k in bedge if len(epat[k]) >= 2}
    qpat = defaultdict(set)
    for qi, q in enumerate(bq):
        for v in q:
            qpat[v].add(qlab[qi])

    # seams and corners
    seams = Seams(seam_polylines(pts, tris, lab))
    h = seams.spacing
    corner, vinfo = set(), {}
    for v in {int(x) for k in feat for x in E[k]}:
        turn, junction = seams.turn(V0[v], h)
        turn2, _ = seams.turn(V0[v], 2 * h)
        kinked = turn is not None and turn > KINK_DEG and turn > 0.75 * (turn2 or 0.0)
        vinfo[v] = {"patches": len(qpat[v]), "seam_turn": turn, "seam_turn_2h": turn2, "junction": junction}
        if len(qpat[v]) >= 3 or junction or turn is None or kinked:
            corner.add(v)

    # continuity pairs, priority feature > 2d > 3d
    pf, af = mutual_pairs(E, sorted(feat), share2, lambda e, f: epat[e] == epat[f])
    pf = {p for p in pf if p[0] not in corner}
    af = [x for x in af if x[0] not in corner]
    nonfeat_b = [k for k in bedge if k not in feat]
    p2, a2 = mutual_pairs(E, nonfeat_b, share2, lambda e, f: epat[e] == epat[f])
    interior = [k for k in range(len(E)) if not share2[k]]
    p3, a3 = mutual_pairs(E, interior, share3)
    pairs, used = {}, set()
    for kind, P in (("feature", pf), ("2d", p2), ("3d", p3)):
        for v, e, f in sorted(P):
            if (v, e) in used or (v, f) in used:
                continue
            pairs[(v, e, f)] = kind
            used |= {(v, e), (v, f)}
    ambiguous = {"feature": af, "2d": a2, "3d": a3}

    # feature arcs, geometry first: every feature edge lies between two patches, and the seam curve with exactly
    # that patch pair is its curve (among several such curves the one nearest to both ends). No distance
    # tolerance: an edge whose patch pair has no seam is a defect of the block structure, never silently a cubic.
    V = V0.copy()
    by_pair = defaultdict(list)
    for c, labs in enumerate(seams.labels):
        by_pair[labs].append(c)
    cur_of, failed = {}, []
    for k in sorted(feat):
        pair = tuple(sorted(epat[k]))
        cands = by_pair.get(pair, [])
        if len(pair) != 2 or not cands:
            failed.append(int(k)); continue
        a, b = E[k]
        cur_of[k] = min(cands, key=lambda c: float(seams.on_curve(c, np.stack([V[a], V[b]]))[0].max()))
    # corners onto their curves: a vertex on one feature curve onto its nearest point, a vertex where several
    # curves meet (a corner of the block structure on the geometry) onto their common end point (seam junction)
    for v in {int(x) for k in cur_of for x in E[k]}:
        cs = sorted({cur_of[k] for k in cur_of if v in E[k]})
        if len(cs) == 1 and v not in corner:
            V[v] = seams.on_curve(cs[0], V[v])[2][0]
            continue
        ends = [seams.pts[c][i] for c in cs if not seams.closed[c] for i in (0, -1)]
        if ends:
            common = [p for p in ends if sum(np.linalg.norm(q - p) < 1e-9 for q in ends) >= min(2, len(cs))]
            pool = common or ends
            V[v] = min(pool, key=lambda p: np.linalg.norm(p - V[v]))
        elif cs:
            V[v] = seams.on_curve(cs[0], V[v])[2][0]
    # corners in the order the block structure prescribes: the feature edges on one seam curve chain its corners
    # (a cycle on a closed curve, a path on an open one); where nearest-point positions contradict that order
    # (two corners a few hundredths apart, swapped -- arcs would then run over each other), the positions are
    # made monotone along the chain with a minimum gap, junction corners fixed
    order_fixed = []
    for c in sorted(set(cur_of.values())):
        ks = [k for k in cur_of if cur_of[k] == c]
        adj = defaultdict(list)
        for k in ks:
            a, b = int(E[k][0]), int(E[k][1]); adj[a].append(b); adj[b].append(a)
        if any(len(n) > 2 for n in adj.values()) or len(adj) < 3:
            continue
        ends = [v for v, n in adj.items() if len(n) == 1]
        start = ends[0] if ends else min(adj)
        chain, prev = [start], None
        while True:
            nx = [w for w in adj[chain[-1]] if w != prev and w not in chain]
            if not nx:
                break
            prev = chain[-1]; chain.append(nx[0])
        if len(chain) != len(adj):
            continue
        L = float(seams.arclen[c][-1]); closed = bool(seams.closed[c])
        t = np.array([float(seams.on_curve(c, V[v])[1][0]) for v in chain])
        fixed = np.array([(v in corner) and not closed for v in chain])
        best = None
        for sgn in (1.0, -1.0):
            if closed:
                d = (sgn * np.diff(np.r_[t, t[0]])) % L               # step to the next corner, cycle closed
                d = np.where(d > 0.5 * L, d - L, d)                   # almost a full turn = a small step back
                if abs(d.sum() - L) > 1e-6 * L:
                    continue                                          # this direction is not one revolution
                u = np.r_[0.0, np.cumsum(d[:-1])]
            else:
                u = sgn * (t - t[0])
            viol = int((np.diff(u) < 0).sum())
            if best is None or viol < best[0]:
                best = (viol, sgn, u)
        if best is None or best[0] == 0:
            continue
        viol, sgn, u = best
        gap = 0.02 * (L if closed else abs(u[-1])) / len(chain)
        w = u.copy()                                                  # pool adjacent violators (isotonic fit)
        blocks = [[i] for i in range(len(w))]
        val = [w[i] for i in range(len(w))]
        i = 0
        while i < len(val) - 1:
            if val[i + 1] < val[i] + gap and not any(fixed[j] for j in blocks[i + 1]):
                blocks[i] += blocks.pop(i + 1); val[i] = float(np.mean(u[blocks[i]])); val.pop(i + 1)
                i = max(i - 1, 0)
            else:
                i += 1
        for b, m in zip(blocks, val):                                 # spread a pooled group around its mean
            for r, j in enumerate(b):
                w[j] = m + (r - (len(b) - 1) / 2) * gap
        for v, wj, uj, f in zip(chain, w, u, fixed):
            if f or abs(wj - uj) < 1e-12:
                continue
            tw = (t[0] + sgn * wj) % L if closed else t[0] + sgn * wj
            V[v] = seams.sample(c, np.array([tw]))[0]
            order_fixed.append(int(v))
    # cubics follow moved corners (chord-relative), feature edges become arcs + G1 cubic
    dV = V - V0
    n1 = c1 + (2 * dV[E[:, 0]] + dV[E[:, 1]]) / 3
    n2 = c2 + (dV[E[:, 0]] + 2 * dV[E[:, 1]]) / 3
    arcs, fit_dev = {}, []
    corners_on = defaultdict(set)
    for k, c in cur_of.items():
        corners_on[c] |= {int(E[k][0]), int(E[k][1])}
    pos_on = {c: {v: float(seams.on_curve(c, V[v])[1][0]) for v in vs} for c, vs in corners_on.items()}
    covered = defaultdict(float)
    for k, c in cur_of.items():
        a, b = E[k]
        hint = bezier(V[a], n1[k], n2[k], V[b], 3)[1]
        others = [t for v, t in pos_on[c].items() if v not in (int(a), int(b))]
        P = seams.arc(c, V[a], V[b], hint, others=others)
        covered[c] += float(np.linalg.norm(np.diff(P, axis=0), axis=1).sum())
        arcs[k] = P
        T0, T1 = seams.tangent(c, V[a]), seams.tangent(c, V[b])     # exact seam tangents at the corners
        T0 = T0 if T0 @ (P[min(4, len(P) - 1)] - P[0]) >= 0 else -T0
        T1 = T1 if T1 @ (P[-1] - P[-min(5, len(P))]) >= 0 else -T1
        n1[k], n2[k] = g1_cubic(P, T0, T1)
        Q = bezier(V[a], n1[k], n2[k], V[b], 129)
        fit_dev.append(float(cKDTree(Q).query(P)[0].max() / max(np.linalg.norm(V[b] - V[a]), 1e-30)))

    # directed outputs, aligned with sample["edges"]
    if verbose and order_fixed:
        print(f"[curve_post] corners re-ordered along their seam: {sorted(set(order_fixed))}", flush=True)
    post_ctrl, arc_pts, arc_off, featd = [], [], [0], []
    for (k, rev) in dmap:
        post_ctrl.append([n2[k], n1[k]] if rev else [n1[k], n2[k]])
        P = arcs.get(k)
        if P is not None:
            arc_pts.append(P[::-1] if rev else P)
        arc_off.append(arc_off[-1] + (len(P) if P is not None else 0))
        featd.append(k in feat)
    role = np.zeros(len(V), np.int64)
    for k in feat:
        role[E[k]] = ROLE["feature"]
    for lst in ambiguous.values():
        for v, _ in lst:
            role[v] = ROLE["singular"]
    for v in corner:
        role[v] = ROLE["corner"]
    other = lambda k, v: int(E[k][1] if E[k][0] == v else E[k][0])  # noqa: E731
    pair_arr = np.array([[v, other(e, v), other(f, v), KIND[kd]]       # v-u continues into v-w
                         for (v, e, f), kd in sorted(pairs.items())], np.int64).reshape(-1, 4)

    kinks_before = {kd: [kink(V0, E, c1, c2, v, e, f) for (v, e, f), k2 in pairs.items() if k2 == kd]
                    for kd in KIND}
    kinks_after = {kd: [kink(V, E, n1, n2, v, e, f) for (v, e, f), k2 in pairs.items() if k2 == kd]
                   for kd in KIND}
    q = {
        "feature_edges": len(feat), "feature_arcs": len(arcs), "corners": len(corner),
        "feature_failed": failed,
        "seam_coverage": {int(c): covered.get(c, 0.0) / max(float(seams.arclen[c][-1]), 1e-30)
                          for c in range(len(seams)) if len(seams.labels[c]) == 2},
        "feature_complete": bool(not failed and all(abs(covered.get(c, 0.0) / max(float(seams.arclen[c][-1]), 1e-30) - 1) < 0.02
                                                    for c in range(len(seams)) if len(seams.labels[c]) == 2)),
        "pairs": dict(Counter(pairs.values())),
        "ambiguous": {k: sum(len(p) for _, p in v) for k, v in ambiguous.items()},
        "vertex_moved_max": float(np.linalg.norm(dV, axis=1).max()),
        "g1_fit_dev_rel_median": float(np.median(fit_dev)) if fit_dev else None,
        "g1_fit_dev_rel_max": float(np.max(fit_dev)) if fit_dev else None,
        "kink_before_median": {k: float(np.median(x)) for k, x in kinks_before.items() if x},
        "kink_before_max": {k: float(np.max(x)) for k, x in kinks_before.items() if x},
        "kink_after_max": {k: float(np.max(x)) for k, x in kinks_after.items() if x},
        "seam_curves": len(seams), "seam_spacing": seams.spacing,
    }
    if verbose:
        print(f"[curve_post] feature edges {q['feature_edges']} (arcs {q['feature_arcs']}, failed {len(q['feature_failed'])}, "
              f"complete {q['feature_complete']}), corners {q['corners']}, "
              f"pairs {q['pairs']}, ambiguous {q['ambiguous']}, kink max before {q['kink_before_max']} "
              f"after {q['kink_after_max']}, G1 fit dev/L max {q['g1_fit_dev_rel_max']}")
    return {
        "post_vertices": V,
        "post_edge_ctrl": np.asarray(post_ctrl, float),
        "post_feature_edge": np.asarray(featd, bool),
        "post_arc_polyline": np.concatenate(arc_pts) if arc_pts else np.zeros((0, 3)),
        "post_arc_offset": np.asarray(arc_off, np.int64),
        "post_vertex_role": role,
        "post_pairs": pair_arr,
        "post_quality": json.dumps(q),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("samples", nargs="+")
    ap.add_argument("--out-dir", default=None, help="write <name>_curves.npz here (default: next to the sample)")
    ap.add_argument("--inplace", action="store_true", help="add the post_* arrays to the sample file itself")
    a = ap.parse_args()
    import export_sample as es
    for path in map(Path, a.samples):
        s = es.load(path)
        post = postprocess(s, verbose=True)
        if a.inplace:
            raw = dict(np.load(path, allow_pickle=False))
            raw.update(post)
            np.savez_compressed(path, **raw)
            print(f"[curve_post] {path}: post_* added")
        else:
            name = path.parent.name if path.name == "sample.npz" else path.stem
            out = (Path(a.out_dir) if a.out_dir else path.parent) / f"{name}_curves.npz"
            out.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(out, **post)
            print(f"[curve_post] {path} -> {out}")


if __name__ == "__main__":
    main()
