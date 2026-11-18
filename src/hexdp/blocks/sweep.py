"""Building block ``section_sweep``: sweep the band-grid region inside every embedded loop (the root section, shared
nodes with the band) outward along the surface patches that start on the loop - e.g. a blade with its fillet from
the hub to the tip.

Inputs: a band (band_layout) and a surface source with ``stages`` = [[patch (nu, nw, 3), ...], ...]: per stage the
patches in cyclic order around the section (the patch joints are the loop's pins), u along a patch, w outward; stage
k + 1 starts where stage k ends. The patches belong to loop 0; loop i is loop 0 turned by 2 pi i / symmetry about z.

Every stage is cut into levels on exact cylinders r = const (the last level of the last stage is the patches' end
w = 1, e.g. a blade tip); every boundary node of the root section keeps its patch and arc fraction on all levels;
the interior of a level is TFI + Winslow on the unrolled cylinder. Blocks: the root section's grid blocks x stages.
Boundary names: <prefix> (side faces), <prefix>_tip (the end level).
"""

from __future__ import annotations

import numpy as np

from hexdp.blocks import Part
from hexdp.blocks.geom import arcfrac, at_fraction, rot_z, tfi, winslow_on_cylinder
from hexdp.contracts import ParamSpec


class SweepError(RuntimeError):
    pass


class SectionSweep:
    """Sweep the band grid inside every embedded loop outward along given surface patches (radial levels)."""

    name = "section_sweep"
    params = {
        "h": ParamSpec(0.04, "fine cell size along the sweep", unit="m", lower=1e-6),
        "grading": ParamSpec("0.5,1", "per stage: exponent of the level spacing (< 1: thick first layers, for a "
                             "fillet running out tangentially)"),
        "prefix": ParamSpec("blade", "boundary name prefix"),
        "normal_cross": ParamSpec(True, "turn the cross lines from the root grid's direction to the mid-line normal "
                                  "over the first stage"),
    }

    def run(self, inputs, params, ctx=None):
        band = next(v for v in inputs.values() if isinstance(v, dict) and v.get("kind") == "band")
        src = next(v for v in inputs.values() if isinstance(v, dict) and "stages" in v)
        grading = [float(x) for x in str(params["grading"]).split(",")]
        return section_sweep(band, src["stages"], h=params["h"], grading=grading, prefix=params["prefix"],
                             normal_cross=params["normal_cross"])


def _region(band, k):
    """Columns / rows of the band grid inside loop k (a rectangle of the tensor grid)."""
    from matplotlib.path import Path as MP

    m, C, P = band["mesh"], band["columns"], band["period"]
    oi, oj, G2 = m.oi, m.oj, band["G2"]
    NI = G2.shape[0]
    poly = MP(band["dom0"].loops[2 + k].closed())
    inside = []
    for c in range(C):
        for r in range(band["rows"]):
            I = (int(oi[c]) + int(oi[c + 1] if c + 1 < C else NI)) // 2 % NI
            J = (int(oj[r]) + int(oj[r + 1])) // 2
            x, z = m.X[G2[I, J]]
            if any(poly.contains_point((x + s * P, z)) for s in (-1, 0, 1)):
                inside.append((c, r))
    cs = sorted({c for c, _ in inside})
    rs = sorted({r for _, r in inside})
    if not inside or len(inside) != len(cs) * len(rs):
        raise SweepError(f"loop {k} is not a rectangle of the band grid: {inside}")
    if max(cs) - min(cs) + 1 != len(cs):                  # across the seam: start after the gap
        first = next(c for c in cs if (c - 1) % C not in cs)
        cs = [(first + j) % C for j in range(len(cs))]
    return cs, rs


def section_sweep(band, stages, h=0.04, grading=(0.5, 1.0), prefix="blade", normal_cross=True):
    m, C, P, R = band["mesh"], band["columns"], band["period"], band["r"]
    NB = int(band.get("symmetry", 1))
    oi, oj, G2 = m.oi, m.oj, band["G2"]
    NI = G2.shape[0]
    n_loops = len(band["dom0"].loops) - 2
    stages = [[np.asarray(f, float) for f in st] for st in stages]
    # radii of the stage ends (mean over the patch ends), levels per stage
    r_end = [float(np.mean([np.hypot(f[:, -1, 0], f[:, -1, 1]).mean() for f in st])) for st in stages]
    r_tip_min = min(float(np.hypot(f[:, -1, 0], f[:, -1, 1]).min()) for f in stages[-1])
    r_start = [R] + r_end[:-1]
    levels = []                                      # per stage: radii of the levels after the stage start
    for k, st in enumerate(stages):
        n = max(2, int(round((r_end[k] - r_start[k]) / h)))
        g = grading[k] if k < len(grading) else 1.0
        rr = r_start[k] + (r_end[k] - r_start[k]) * np.linspace(0, 1, n + 1) ** g
        if k == len(stages) - 1:                     # the end level is the real tip (w = 1), not a cylinder
            rr[-2] = min(rr[-2], r_tip_min - 1e-3 * (r_end[k] - R))
            rr = np.r_[rr[:-1], np.nan]
        levels.append(rr[1:])
    pts = list(band["points"])
    grids = []
    for k in range(n_loops):
        ang = 2 * np.pi * (k % max(NB, 1)) / max(NB, 1) if NB > 1 else 0.0
        if NB > 1 and n_loops != NB:
            raise SweepError(f"{n_loops} loops for symmetry {NB}: one loop per symmetric copy expected")
        st_k = [[rot_z(f, ang) for f in st] for st in stages]
        cs, rs = _region(band, k)
        Is = []
        for c in cs:
            a, b = int(oi[c]), int(oi[c + 1]) if c + 1 < C else NI
            Is += [i % NI for i in range(a, b)]
        Is.append(int(oi[cs[-1] + 1]) % NI if cs[-1] + 1 < C else 0)
        Js = list(range(int(oj[rs[0]]), int(oj[rs[-1] + 1]) + 1))
        S0 = G2[np.ix_(Is, Js)]
        ni, nj = S0.shape
        # root boundary node -> (patch, arc fraction) on the first stage's patch curves at w = 0
        pieces0 = []
        for f in st_k[0]:
            c3 = f[:, 0]
            pieces0.append(np.c_[R * np.unwrap(np.arctan2(c3[:, 1], c3[:, 0])), c3[:, 2]])
        bidx = ([(i, 0) for i in range(ni)] + [(ni - 1, j) for j in range(1, nj)]
                + [(i, nj - 1) for i in range(ni - 2, -1, -1)] + [(0, j) for j in range(nj - 2, 0, -1)])
        loc = {}
        for (i, j) in bidx:
            x = m.X[S0[i, j]]
            best = None
            for pi, Q in enumerate(pieces0):
                for s_ in (-1, 0, 1):
                    d = np.linalg.norm(Q + [s_ * P, 0] - x, axis=1)
                    qn = int(np.argmin(d))
                    if best is None or d[qn] < best[0]:
                        best = (d[qn], pi, arcfrac(Q)[qn])
            loc[(i, j)] = best[1:]

        def at_r(F, r):
            if np.isnan(r):
                return F[:, -1]
            Rw = np.hypot(F[..., 0], F[..., 1])
            return np.array([[np.interp(r, Rw[u], F[u, :, c]) for c in range(3)] for u in range(F.shape[0])])

        ic = [0]
        for c in cs:
            ic.append(ic[-1] + int((oi[c + 1] if c + 1 < C else NI) - oi[c]))
        jc = [int(oj[r] - oj[rs[0]]) for r in rs] + [int(oj[rs[-1] + 1] - oj[rs[0]])]

        def side_index(curves, off, ij):
            pi, t = loc[ij]
            return off[pi] + t * (len(curves[pi]) - 1)

        def reorient(L, curves, w):
            """Turn the interior cross lines (rows jc between the caps) towards the normal of the mid line by the
            weight w: their end points move along the long sides (i = 0, ni-1) to the intersections of the mid-line
            normal at the cross line's middle, blended with the patch-fraction position; the other side nodes keep
            their spacing between the new anchors. Returns the new boundary (interior left to TFI/Winslow)."""
            if w <= 0 or len(jc) < 3:
                return L
            off = np.cumsum([0] + [len(c) - 1 for c in curves])
            loop = np.vstack([c[:-1] for c in curves])
            n_lp = len(loop)
            rr = float(np.hypot(loop[:, 0], loop[:, 1]).mean())
            ref = np.arctan2(L[0, 0, 1], L[0, 0, 0])

            def u2(Y):                                   # unrolled 2D about the section
                ph = np.arctan2(Y[..., 1], Y[..., 0]) - ref
                ph = (ph + np.pi) % (2 * np.pi) - np.pi
                return np.stack([rr * ph, Y[..., 2]], -1)
            L0 = winslow_on_cylinder(tfi(L.copy()))
            mid = ni // 2
            new = L.copy()
            for i_side in (0, ni - 1):
                # the side as a loop path from the j=0 corner to the j=nj-1 corner (through the side's middle)
                k0 = side_index(curves, off, (i_side, 0)) % n_lp
                k1 = side_index(curves, off, (i_side, nj - 1)) % n_lp
                km = side_index(curves, off, (i_side, nj // 2)) % n_lp
                span_f = (k1 - k0) % n_lp
                fwd = (km - k0) % n_lp <= span_f
                steps = span_f if fwd else (k0 - k1) % n_lp
                q = np.linspace(0, steps, int(np.ceil(steps)) * 4 + 2)
                kk = (k0 + q) % n_lp if fwd else (k0 - q) % n_lp
                lo = np.floor(kk).astype(int)
                fr = (kk - lo)[:, None]
                path3 = (1 - fr) * loop[lo % n_lp] + fr * loop[(lo + 1) % n_lp]
                path2 = u2(path3)
                arc = np.r_[0, np.cumsum(np.linalg.norm(np.diff(path2, axis=0), axis=1))]
                # current anchors (arc positions of the side nodes) and the new ones of the cross lines
                node2 = u2(L[i_side])
                s_nodes = np.array([arc[np.argmin(np.linalg.norm(path2 - p_, axis=1))] for p_ in node2])
                s_new = s_nodes.copy()
                for jj in jc[1:-1]:
                    C2 = u2(L0[mid, jj])
                    tvec = u2(L0[mid, min(jj + 1, nj - 1)]) - u2(L0[mid, max(jj - 1, 0)])
                    nvec = np.array([-tvec[1], tvec[0]])
                    d = (path2 - C2) @ (tvec / (np.linalg.norm(tvec) + 1e-300))   # along the mid line
                    sgn = np.sign(d)
                    hits = np.flatnonzero(sgn[:-1] * sgn[1:] <= 0)
                    if not len(hits):
                        continue
                    h = hits[np.argmin(np.abs(arc[hits] - s_nodes[jj]))]
                    t_ = d[h] / (d[h] - d[h + 1]) if d[h] != d[h + 1] else 0.0
                    s_hit = arc[h] + t_ * (arc[h + 1] - arc[h])
                    s_new[jj] = (1 - w) * s_nodes[jj] + w * s_hit
                # the other nodes keep their relative spacing between consecutive anchors
                anchors = [0, *jc[1:-1], nj - 1]
                for a_, b_ in zip(anchors[:-1], anchors[1:]):
                    if s_nodes[b_] == s_nodes[a_]:
                        continue
                    for j in range(a_ + 1, b_):
                        f = (s_nodes[j] - s_nodes[a_]) / (s_nodes[b_] - s_nodes[a_])
                        s_new[j] = s_new[a_] + f * (s_new[b_] - s_new[a_])
                if np.any(np.diff(s_new) <= 0):            # order lost: keep the patch fractions on this side
                    continue
                for j in range(1, nj - 1):
                    new[i_side, j] = np.array([np.interp(s_new[j], arc, path3[:, c]) for c in range(3)])
            return new

        def level(faces, r, w=0.0):
            L = np.zeros((ni, nj, 3))
            curves = [at_r(F, r) for F in faces]
            for (i, j), (pi, t) in loc.items():
                L[i, j] = at_fraction(curves[pi], np.array([t]))[:, 0]
            return winslow_on_cylinder(tfi(reorient(L, curves, w)))
        n_lv = sum(len(lv) for lv in levels)
        ids = np.zeros((ni, nj, n_lv + 1), int)
        ids[:, :, 0] = S0
        lv_i, bounds = 0, [0]
        r_first_end = r_end[0]
        for sk, (faces, rr) in enumerate(zip(st_k, levels)):
            for r in rr:
                lv_i += 1
                # cross lines turn from the hub grid's direction (root) to the mid-line normal over the first stage
                w = 1.0 if (sk > 0 or np.isnan(r)) else float(np.clip((r - R) / (r_first_end - R), 0, 1))
                L = level(faces, r, w if normal_cross else 0.0)
                for i in range(ni):
                    for j in range(nj):
                        pts.append(L[i, j])
                        ids[i, j, lv_i] = len(pts) - 1
            bounds.append(lv_i)
        for a in range(len(cs)):
            for b in range(len(rs)):
                for k0, k1 in zip(bounds[:-1], bounds[1:]):
                    grids.append(ids[ic[a]:ic[a + 1] + 1, jc[b]:jc[b + 1] + 1, k0:k1 + 1])
    tip_r = 0.98 * r_tip_min

    def label(F):
        r = np.hypot(F[:, 0], F[:, 1])
        if np.abs(r - R).max() < 1e-6:
            return None                              # the root: interior once composed with the band's solid
        return f"{prefix}_tip" if r.min() >= tip_r else prefix
    return Part(prefix, np.array(pts), grids, label, {"blocks": len(grids), "stages": len(stages),
                                                      "levels": [len(lv) for lv in levels]})
