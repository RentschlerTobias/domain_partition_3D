"""Stage 7 (T5): block edges as cubic Bezier curves.

A block edge in `<name>_blocks.vtk` is already a polyline of fine-mesh
vertices -- roughly 25 points at n=60000 -- and it tracks the geometry to
0.006. Fitting one cubic through it is what makes the exported sample
higher-order: two absolute control points per directed edge instead of a
straight chord.

Why the polyline and not the frame field: decision D. The field is never
written out, it is a per-tet direction rather than a curve, and it is the stage
that produced this mesh in the first place.

Three things this module reports, and the third is the one that decides `-n`:

* the fit residual, defined exactly as the 2D study defined it, so the numbers
  are comparable: per edge, the MAXIMUM distance from the polyline to the
  fitted curve divided by the chord length, then median / p95 / max over
  edges. 2D reference: 0.74 % median over 330 164 edges
  (`meshtron/docs/ho_quad_transformer/06_edge_geometry_study.md`).
* the inflection count, the 3D analogue of the 2D sign-change test.
* `--subsample K`, which is what makes the residual honest at low `-n`. See
  the warning on `fit_edges`.
"""

import json
from collections import defaultdict
from pathlib import Path

import numpy as np

import tfi

REPO = Path(__file__).resolve().parent.parent.parent
DEFAULT_BLOCKS = (REPO / "output" / "hex3d_algohex" / "deliverable"
                  / "T1_9_blocks_v11.vtk")


# --------------------------------------------------------------------------
# the polylines
# --------------------------------------------------------------------------

def edge_chains(lat):
    """Every block edge of the complex as a chain of mesh vertex indices.

    A block edge is a lattice axis line: for one of the three axes, the line
    of vertices at a fixed corner of the other two. Twelve per block, and a
    shared edge yields the SAME chain from both blocks, so the dedup key is
    the chain itself.

    Keyed on the whole chain, not on its endpoints. Two distinct block edges
    can join the same pair of corners along different curves -- decision D2
    keeps them apart with half-edges for exactly this reason -- and any
    endpoint-based key would silently merge them.

    Returns {canonical chain: [(block, axis), ...]} so a caller can see how
    many blocks claim each edge: 1 on the boundary of the complex, more
    inside.
    """
    owners = defaultdict(list)
    for r, (_dims, vert) in lat.items():
        for axis in (0, 1, 2):
            other = [a for a in (0, 1, 2) if a != axis]
            for c0 in (0, -1):
                for c1 in (0, -1):
                    sl = [slice(None)] * 3
                    sl[other[0]] = c0
                    sl[other[1]] = c1
                    chain = tuple(int(v) for v in vert[tuple(sl)])
                    owners[_canonical(chain)].append((int(r), axis))
    return dict(owners)


def _canonical(chain):
    """One of the two traversals, picked so both blocks agree on the key."""
    return chain if chain[0] <= chain[-1] else chain[::-1]


# --------------------------------------------------------------------------
# the fit
# --------------------------------------------------------------------------

def fit_cubic_bezier(P0, P1, pts):
    """Best-fit cubic Bezier with fixed endpoints. Returns (B1, B2, mode).

    Lifted from `meshtron/prototype_twostage.py::TwoStageTokenizer.
    _fit_cubic_bezier`, and the lift itself was free -- the 2D version is
    already dimension-agnostic, because the only shape-bearing operation is an
    lstsq whose right-hand side carries the coordinate axis. Kept as a separate
    function rather than imported so this repository does not grow a dependency
    on meshtron; the seam runs the other way (decision C).

    The parameter is CHORD-LENGTH along the polyline, not index. Uniform
    parameterisation on a polyline whose segments differ in length pulls the
    curve towards the densely sampled end.

    **What had to be added for 3D, and it is not a 3D issue at all.** Two free
    control vectors need two interior points. In 2D that was never binding --
    the streamlines carry many samples -- but a block edge at n=2000 carries
    about four points total, so the underdetermined case is the normal case
    here, and the 2D code fails silently in it: with two points both control
    points come back as the ORIGIN (`b1 = b2 = 0`, so `A` is the zero matrix
    and lstsq returns the zero vector), and the residual check cannot see it,
    because the only two points being measured are the fixed endpoints, which
    lie on any such curve by construction. A perfect score for a curve through
    the origin.

    So the model is chosen by how much data there is, and the choice is
    reported:

    * ``cubic`` -- >= 2 interior points, both control vectors free.
    * ``quadratic`` -- exactly 1 interior point. A quadratic Bezier through it,
      degree-elevated to cubic. One free control vector, exactly determined,
      and a real arc rather than lstsq's minimum-norm artifact.
    * ``chord`` -- no interior point. The straight line, honestly labelled.
    """
    pts = np.asarray(pts, float)
    P0 = np.asarray(P0, float)
    P1 = np.asarray(P1, float)
    n_int = len(pts) - 2

    if n_int <= 0:
        return P0 + (P1 - P0) / 3.0, P0 + 2.0 * (P1 - P0) / 3.0, "chord"

    d = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    cum = np.concatenate([[0.0], np.cumsum(d)])
    t = cum / cum[-1] if cum[-1] > 1e-12 else np.linspace(0, 1, len(pts))

    if n_int == 1:
        # quadratic in Bernstein form, then elevated: B1 = P0 + 2/3 (Q - P0),
        # B2 = P1 + 2/3 (Q - P1)
        q = (2 * (1 - t) * t)[:, None]
        rhs = pts - ((1 - t) ** 2)[:, None] * P0 - (t ** 2)[:, None] * P1
        Q, *_ = np.linalg.lstsq(q, rhs, rcond=None)
        Q = Q[0]
        return (P0 + 2.0 / 3.0 * (Q - P0), P1 + 2.0 / 3.0 * (Q - P1),
                "quadratic")

    b1 = 3 * (1 - t) ** 2 * t
    b2 = 3 * (1 - t) * t ** 2
    rhs = pts - ((1 - t) ** 3)[:, None] * P0 - (t ** 3)[:, None] * P1
    A = np.stack([b1, b2], axis=1)
    sol, *_ = np.linalg.lstsq(A, rhs, rcond=None)
    return sol[0], sol[1], "cubic"


def bezier_points(P0, B1, B2, P1, n=200):
    t = np.linspace(0, 1, n)[:, None]
    return ((1 - t) ** 3 * P0 + 3 * (1 - t) ** 2 * t * B1
            + 3 * (1 - t) * t ** 2 * B2 + t ** 3 * P1)


def max_dist_to_curve(pts, curve):
    """Largest distance from any polyline vertex to the fitted curve.

    Distance to the curve as a curve -- point to nearest SEGMENT of the dense
    sampling, not point to nearest sample. Nearest-sample would report a
    spurious error of up to half a sample spacing and would improve when the
    sampling is refined, which is not a property a residual may have.
    """
    A = curve[:-1]
    d = curve[1:] - A
    L2 = np.einsum("ij,ij->i", d, d)
    L2 = np.where(L2 > 0, L2, 1.0)
    w = pts[:, None, :] - A[None, :, :]
    s = np.clip(np.einsum("pij,ij->pi", w, d) / L2, 0.0, 1.0)
    foot = A[None, :, :] + s[:, :, None] * d[None, :, :]
    return float(np.linalg.norm(pts[:, None, :] - foot, axis=2).min(1).max())


# --------------------------------------------------------------------------
# shape descriptors: is one cubic even the right model?
# --------------------------------------------------------------------------

def _in_plane(pts):
    """The polyline in the 2D coordinates of its own best-fit plane."""
    Q = pts - pts.mean(0)
    V = np.linalg.svd(Q, full_matrices=False)[2]
    return Q @ V[0], Q @ V[1]


def _sign_runs(v, mag, floor):
    """Sign changes between consecutive runs that clear `floor`.

    The threshold is applied per RUN, not per sample, which is the difference
    between counting shape and counting jitter: a stretch that wanders across
    zero at the 1e-6 level contributes nothing, while a genuine excursion to
    one side and back counts once.
    """
    runs = []
    for s, m in zip(np.sign(v), mag):
        if runs and runs[-1][0] == s:
            runs[-1][1] = max(runs[-1][1], m)
        else:
            runs.append([s, m])
    keep = [s for s, m in runs if s != 0 and m > floor]
    return int(sum(1 for a, b in zip(keep[:-1], keep[1:]) if a != b))


def inflection_count(pts, tol=0.01):
    """How often the polyline genuinely crosses its own chord.

    The operational meaning of "single arc", and the metric this module
    reports: a monotone bow stays on one side of its chord, an S-curve does
    not. Runs in the polyline's own best-fit plane -- justified by
    `planarity`, which reports what that plane leaves out.

    `tol` is a fraction of chord length: an excursion counts only if the
    polyline leaves the chord by more than 1 % of it, the same order as the
    residual the cubic achieves, and therefore the scale below which "which
    side" is not a meaningful question.

    **This is NOT the 2D study's test**, and the substitution is measured, not
    preferred -- see `curvature_sign_changes`.
    """
    pts = np.asarray(pts, float)
    if len(pts) < 4:
        return 0
    chord = float(np.linalg.norm(pts[-1] - pts[0]))
    if chord <= 0:
        return 0
    x, y = _in_plane(pts)
    axis = np.array([x[-1] - x[0], y[-1] - y[0]])
    if np.linalg.norm(axis) <= 0:
        return 0
    perp = np.array([-axis[1], axis[0]]) / np.linalg.norm(axis)
    h = np.stack([x[1:-1] - x[0], y[1:-1] - y[0]], axis=1) @ perp
    return _sign_runs(h, np.abs(h) / chord, tol)


def curvature_sign_changes(pts, tol=0.08, smooth=3):
    """The 2D study's inflection test, ported and then REJECTED. Do not use
    this as the reported metric; `inflection_count` replaced it.

    `06_edge_geometry_study.md` counts inflections as sign changes of the
    smoothed curvature, and finds 0 over 330 164 2D streamline edges. Ported
    here it reports 120 inflections on 37 of v11's 107 edges -- and that
    number is noise, which took three attempts to establish:

    * unsmoothed, 315 on 47 edges;
    * smoothed but thresholded against each edge's OWN peak curvature, 132 --
      circular, since on a straight edge the peak is itself noise;
    * smoothed and thresholded against ``|k| * chord**2``, still 120.

    So it was measured directly against a ground truth instead of tuned.
    Splitting v11's edges by `inflection_count`, on the **85 edges that
    provably do not cross their chord**, ``|k| * chord**2`` has median 0.275,
    p95 1.354 and reaches **5.995** -- while the 4 genuine S-curves only reach
    **2.740**. The noise does not merely overlap the signal, it exceeds it, so
    no threshold separates them and the statistic has no discriminating power
    here.

    The cause is not 3D and not the mesh being bad. A block edge is ~17
    fine-mesh vertices whose deviation from their own 3-point smoothing is
    7.25e-3 median against a chord of 0.72 -- about 1 %. Curvature
    differentiates twice, amplifying that by 1/h^2 with h ~ chord/17, i.e. by
    ~290. The chord-excursion test integrates instead, which is why it
    survives. The 2D study was measuring smooth streamline samples, not a
    mesh's own vertices; that is the difference, and it is why the 2D method
    does not transfer.

    Kept because a later session with a smoother edge source (say, projected
    onto the input triangulation -- Option C) may reasonably want it back.
    """
    pts = np.asarray(pts, float)
    if len(pts) < 5:
        return 0
    chord = float(np.linalg.norm(pts[-1] - pts[0]))
    if chord <= 0:
        return 0
    x, y = _in_plane(pts)
    u = np.stack([np.diff(x), np.diff(y)], axis=1)
    a, b = u[:-1], u[1:]
    cross = a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]
    den = (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1)
           * np.linalg.norm(a + b, axis=1))
    k = np.where(den > 0, 2.0 * cross / np.where(den > 0, den, 1.0), 0.0)
    if smooth > 1 and len(k) >= smooth:
        k = np.convolve(k, np.ones(smooth) / smooth, mode="valid")
    return _sign_runs(k, np.abs(k) * chord ** 2, tol)


def planarity(pts):
    """Deviation from the best-fit plane, as a fraction of chord length.

    New information in 3D that the 2D study could not have had: a cubic
    Bezier is planar iff its four control points are coplanar, so a
    non-planar edge is one that a single cubic can only approximate. If this
    stays small, the 2D conclusions transfer; if not, that is the finding.
    """
    pts = np.asarray(pts, float)
    chord = float(np.linalg.norm(pts[-1] - pts[0]))
    if len(pts) < 4 or chord <= 0:
        return 0.0
    Q = pts - pts.mean(0)
    n = np.linalg.svd(Q, full_matrices=False)[2][-1]
    return float(np.abs(Q @ n).max() / chord)


# --------------------------------------------------------------------------
# per-edge records
# --------------------------------------------------------------------------

def fit_edges(P, chains, subsample=0, n_curve=200):
    """One record per UNDIRECTED edge; the exporter makes them directed.

    `subsample` is the honest-residual switch, and T6 is why it exists. With
    4 polyline points, two of them the fixed endpoints, a cubic with two free
    control points has exactly as many unknowns as data: the residual is zero
    BY CONSTRUCTION, not because the curve is right. So a residual measured at
    n=2000 cannot be compared with one measured at n=60000 -- the low-n number
    is a statement about degrees of freedom, not about geometry.

    `--subsample K` fits from K points evenly spaced along the polyline and
    then measures against ALL of them. That asks the question T6 actually
    means: what does a coarse polyline cost, holding the geometry fixed?
    """
    out = []
    for chain in chains:
        pts = P[list(chain)]
        chord = float(np.linalg.norm(pts[-1] - pts[0]))
        arc = float(np.linalg.norm(np.diff(pts, axis=0), axis=1).sum())
        if subsample and len(pts) > subsample:
            idx = np.unique(np.linspace(0, len(pts) - 1, subsample).astype(int))
            src = pts[idx]
        else:
            src = pts
        B1, B2, mode = fit_cubic_bezier(pts[0], pts[-1], src)
        curve = bezier_points(pts[0], B1, B2, pts[-1], n_curve)
        err = max_dist_to_curve(pts, curve)
        # what the straight chord would have cost on the same edge: without
        # it, a small residual cannot be distinguished from a straight edge
        chord_err = max_dist_to_curve(
            pts, bezier_points(pts[0], pts[0] + (pts[-1] - pts[0]) / 3.0,
                               pts[0] + 2.0 * (pts[-1] - pts[0]) / 3.0,
                               pts[-1], n_curve))
        out.append({
            "chain": [int(v) for v in chain],
            "n_pts": len(pts),
            "n_fit_pts": len(src),
            "fit_mode": mode,
            "chord": chord,
            "arc": arc,
            "sagitta_over_chord": (arc / chord - 1.0) if chord > 0 else 0.0,
            "ctrl": [B1.tolist(), B2.tolist()],
            "err_abs": err,
            "err_rel": (err / chord) if chord > 0 else float("nan"),
            "chord_err_rel": (chord_err / chord) if chord > 0 else float("nan"),
            "inflections": inflection_count(pts),
            "planarity": planarity(pts),
        })
    return out


def _pct(v, q):
    return float(np.percentile(v, q)) if len(v) else float("nan")


def summary(recs):
    e = np.array([r["err_rel"] for r in recs], float)
    e = e[np.isfinite(e)]
    ce = np.array([r["chord_err_rel"] for r in recs], float)
    ce = ce[np.isfinite(ce)]
    a = np.array([r["err_abs"] for r in recs], float)
    npts = np.array([r["n_pts"] for r in recs], int)
    pl = np.array([r["planarity"] for r in recs], float)
    infl = np.array([r["inflections"] for r in recs], int)
    ch = np.array([r["chord"] for r in recs], float)
    # arc / chord, the predictor of fit failure in 3D -- see print_summary
    ac = np.array([(r["arc"] / r["chord"]) if r["chord"] > 0 else np.nan
                   for r in recs], float)
    ac = ac[np.isfinite(ac)]
    modes = defaultdict(int)
    for r in recs:
        modes[r["fit_mode"]] += 1
    return {
        "fit_mode_cubic": modes["cubic"],
        "fit_mode_quadratic": modes["quadratic"],
        "fit_mode_chord": modes["chord"],
        "chord_err_rel_median": _pct(ce, 50),
        "chord_err_rel_p95": _pct(ce, 95),
        "arc_over_chord_median": _pct(ac, 50),
        "arc_over_chord_p95": _pct(ac, 95),
        "arc_over_chord_max": float(ac.max()) if len(ac) else float("nan"),
        "edges_over_5pct": int((e > 0.05).sum()),
        "n_edges": len(recs),
        "pts_per_edge_min": int(npts.min()) if len(npts) else 0,
        "pts_per_edge_median": float(np.median(npts)) if len(npts) else 0.0,
        "pts_per_edge_max": int(npts.max()) if len(npts) else 0,
        "fit_pts_median": float(np.median([r["n_fit_pts"] for r in recs]))
        if recs else 0.0,
        "err_rel_median": _pct(e, 50),
        "err_rel_p95": _pct(e, 95),
        "err_rel_max": float(e.max()) if len(e) else float("nan"),
        "err_abs_median": _pct(a, 50),
        "err_abs_max": float(a.max()) if len(a) else float("nan"),
        "chord_min": float(ch.min()) if len(ch) else float("nan"),
        "chord_median": _pct(ch, 50),
        "planarity_median": _pct(pl, 50),
        "planarity_p95": _pct(pl, 95),
        "planarity_max": float(pl.max()) if len(pl) else float("nan"),
        "edges_with_inflection": int((infl > 0).sum()),
        "inflections_total": int(infl.sum()),
        # edges too short to say anything about: the 2D study's worst case was
        # exactly this, a chord of 0.016, and it dominated the max
        "degenerate_edges": int((ch < 0.02 * float(np.median(ch))).sum())
        if len(ch) else 0,
    }


def print_summary(name, s):
    print(f"\n[block_edges] {name}")
    print(f"  edges                {s['n_edges']}")
    print(f"  points per edge      min {s['pts_per_edge_min']}, "
          f"median {s['pts_per_edge_median']:.0f}, max {s['pts_per_edge_max']}"
          f"   (fitted from {s['fit_pts_median']:.0f})")
    print(f"  fit model            {s['fit_mode_cubic']} cubic, "
          f"{s['fit_mode_quadratic']} quadratic, "
          f"{s['fit_mode_chord']} straight (too few points)")
    print(f"  residual / chord     median {s['err_rel_median'] * 100:.3f} %, "
          f"p95 {s['err_rel_p95'] * 100:.3f} %, "
          f"max {s['err_rel_max'] * 100:.2f} %")
    print(f"  straight would cost  median "
          f"{s['chord_err_rel_median'] * 100:.3f} %, "
          f"p95 {s['chord_err_rel_p95'] * 100:.3f} %"
          f"   <- what the cubic buys")
    print(f"  residual absolute    median {s['err_abs_median']:.3e}, "
          f"max {s['err_abs_max']:.3e}")
    print(f"  chord                min {s['chord_min']:.4f}, "
          f"median {s['chord_median']:.4f}"
          f"   ({s['degenerate_edges']} degenerate)")
    # An edge whose arc is several times its chord is wound around something,
    # and one cubic cannot represent it at any residual -- a hard geometric
    # limit, not a fit failure. Measured on cand_001: the four edges above 5 %
    # residual are exactly the four with the highest arc/chord (4.51, 3.06,
    # 1.66, 1.55) and the lowest planarity, while the median edge sits at
    # 1.001. Chord LENGTH does not predict it -- the worst is a well-sampled
    # 26-point edge at a third of the median chord. That is the difference
    # from 2D, where the worst case was a degenerate mini-edge.
    print(f"  arc / chord          median {s['arc_over_chord_median']:.3f}, "
          f"p95 {s['arc_over_chord_p95']:.3f}, "
          f"max {s['arc_over_chord_max']:.3f}"
          f"   ({s['edges_over_5pct']} edges over 5 % residual)")
    print(f"  planarity / chord    median {s['planarity_median'] * 100:.3f} %, "
          f"p95 {s['planarity_p95'] * 100:.3f} %, "
          f"max {s['planarity_max'] * 100:.2f} %")
    print(f"  inflections          {s['inflections_total']} on "
          f"{s['edges_with_inflection']} of {s['n_edges']} edges "
          f"(chord-excursion test; see curvature_sign_changes)")


# --------------------------------------------------------------------------
# a VTK of the fitted curves, because a number is not a picture
# --------------------------------------------------------------------------

def write_curves_vtk(path, recs, n=24):
    """The fitted curves as polylines, coloured by relative residual x 1e6.

    The branch convention: VTK is authoritative, so a fit that looks right in
    the table and wrong in ParaView is caught here rather than at T9.
    """
    pts, cells, val = [], [], []
    for r in recs:
        B1, B2 = (np.array(c, float) for c in r["ctrl"])
        p0 = np.array(r["_p0"], float)
        p1 = np.array(r["_p1"], float)
        c = bezier_points(p0, B1, B2, p1, n)
        base = len(pts)
        pts.extend(c.tolist())
        cells.append(list(range(base, base + n)))
        val.append(int(round(r["err_rel"] * 1e6))
                   if np.isfinite(r["err_rel"]) else -1)
    with open(path, "w") as f:
        f.write("# vtk DataFile Version 2.0\nfitted block-edge cubics\n"
                "ASCII\nDATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {len(pts)} float\n")
        for p in pts:
            f.write(f"{p[0]} {p[1]} {p[2]}\n")
        tot = sum(len(c) + 1 for c in cells)
        f.write(f"CELLS {len(cells)} {tot}\n")
        for c in cells:
            f.write(f"{len(c)} " + " ".join(str(i) for i in c) + "\n")
        f.write(f"CELL_TYPES {len(cells)}\n")
        for _ in cells:
            f.write("4\n")                                  # VTK_POLY_LINE
        f.write(f"CELL_DATA {len(cells)}\n"
                "SCALARS err_rel_x1e6 int 1\nLOOKUP_TABLE default\n")
        for v in val:
            f.write(f"{v}\n")


# --------------------------------------------------------------------------

def run(blocks_vtk, subsample=0, verbose=True):
    """(records, summary, owners) for one written block structure."""
    P, H, B, f2h = tfi.load_blocks(str(blocks_vtk))
    lat, missing = tfi.lattices(P, H, B, f2h, verbose=verbose)
    owners = edge_chains(lat)
    chains = sorted(owners)
    recs = fit_edges(P, chains, subsample=subsample)
    for r in recs:
        r["_p0"] = P[r["chain"][0]].tolist()
        r["_p1"] = P[r["chain"][-1]].tolist()
        r["n_blocks"] = len(owners[tuple(r["chain"])])
    s = summary(recs)
    s["blocks"] = len(lat)
    s["blocks_without_lattice"] = len(missing)
    if verbose:
        print(f"[block_edges] {len(lat)} lattice blocks, "
              f"{len(chains)} distinct block edges "
              f"(expected 12 per block before dedup: "
              f"{12 * len(lat)} -> {len(chains)})")
    return recs, s, owners


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("blocks", nargs="?", default=str(DEFAULT_BLOCKS))
    ap.add_argument("--subsample", type=int, default=0,
                    help="fit from this many points, measure against all of "
                         "them. The control that makes low-n residuals "
                         "comparable -- see fit_edges")
    ap.add_argument("--json", default=None, help="write the per-edge records")
    ap.add_argument("--vtk", default=None, help="write the fitted curves")
    a = ap.parse_args()

    recs, s, _own = run(a.blocks, subsample=a.subsample)
    print_summary(Path(a.blocks).name
                  + (f"  [subsampled to {a.subsample}]" if a.subsample else ""),
                  s)
    if a.json:
        Path(a.json).parent.mkdir(parents=True, exist_ok=True)
        with open(a.json, "w") as f:
            json.dump({"source": str(a.blocks), "subsample": a.subsample,
                       "summary": s, "edges": recs}, f)
        print(f"wrote {a.json}")
    if a.vtk:
        Path(a.vtk).parent.mkdir(parents=True, exist_ok=True)
        write_curves_vtk(a.vtk, recs)
        print(f"wrote {a.vtk}")
