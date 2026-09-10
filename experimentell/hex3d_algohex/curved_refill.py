"""How much does a cubic block edge buy over a linear one, and when?

`tfi.refill_block` resamples a block's boundary faces with `_lerp_axis`:
piecewise-linear interpolation through the existing grid points. At h = 0.05
that is measurably fine -- the chordal term is about 0.0005 of a total
boundary error of 0.034, so the sampling is not what limits fidelity there.

That measurement does NOT carry over to the regime this is built for. The
plan-B route wants FEW, LARGE blocks and few cells, and a coarse sampling of a
curved edge cuts corners in proportion to how coarse it is. The fix is the one
Meshtron's plan-B branch already uses in 2D (`geom_head_prototype.py`, stage
3): represent the edges at higher order -- there as cubic Bezier control
points per half-edge, regressed and scored against the true streamline -- and
sample the curve rather than the chords.

Higher order lives in the CONSTRUCTION, not in the output: the cells written
are still linear hexes, so a finite-volume solver is unaffected.

Measured on v11, boundary error against the input surface, p95 over points
sampled INSIDE the new quads:

    divisions   chord p95   coons p95   gain
      100 %      0.08172     0.01059    7.7x
       50 %      0.08009     0.01103    7.3x
       25 %      0.08056     0.01533    5.2x
     12.5 %      0.07895     0.04065    1.9x
      6.2 %      0.07263     0.07022    1.0x

So in the GENERATIVE setting the curved edges are worth a factor of 5 to 8,
and 0.011 is close to what the full AlgoHex mesh itself achieves (0.006).
But the gain needs cells to show: below about five cells per block edge it
vanishes, because the SURFACE between samples is flat again however good the
curve is. Few blocks, yes; arbitrarily few cells, no.

The four modes, and which question each answers:

    linear  what refill_block does today: piecewise-linear through the fine
            grid. Already tracks the surface, because the fine grid does.
    cubic   the same with a cubic through the fine grid. Worth 1.00-1.01x --
            measured, and the reason is that both place the new vertices on
            the same fine polyline.
    chord   the generative case WITHOUT edge geometry: straight lines between
            the four block corners.
    coons   the generative case WITH it: cubic edges, Coons patch inside.

    PY curved_refill.py BLOCKS.vtk --input-vtk INPUT.vtk
"""

import argparse
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "experimentell" / "hex3d_algohex"))

import clean_blocks as cb                                              # noqa: E402
import tfi                                                             # noqa: E402


def _cubic_axis(Q, frac, axis):
    """Cubic-spline resampling of `Q` along `axis` at normalised positions.

    Same parameterisation as `tfi._lerp_axis` -- normalised index, not arc
    length -- so that the two differ only in interpolation order. Falls back
    to linear where there are too few points for a cubic (n < 4), which is
    exactly the coarse case where the two agree anyway."""
    from scipy.interpolate import CubicSpline
    n = Q.shape[axis]
    if n < 4:
        return tfi._lerp_axis(Q, frac, axis)
    t = np.linspace(0.0, 1.0, n)
    cs = CubicSpline(t, Q, axis=axis)
    return cs(np.asarray(frac, float))


def _coons(e0, e1, e2, e3):
    """Coons patch from four boundary curves: bilinearly blended edges minus
    the bilinear corner term. The 2D case of the Gordon-Hall map `tfi.tfi`
    already applies in 3D."""
    n, m = len(e0), len(e2)
    u = np.linspace(0, 1, n)[:, None, None]
    v = np.linspace(0, 1, m)[None, :, None]
    F = ((1 - v) * e0[:, None, :] + v * e1[:, None, :]
         + (1 - u) * e2[None, :, :] + u * e3[None, :, :])
    # e0 = S(u,0) and e1 = S(u,1), so e0[0]=S(0,0), e0[-1]=S(1,0),
    # e1[0]=S(0,1), e1[-1]=S(1,1). Pairing those with the wrong blend term --
    # v against S(1,0) instead of S(0,1) -- is a transposed corner grid, and
    # it does not fail loudly: the patch still interpolates its four corners
    # and is merely wrong in between, which measured as a boundary error of
    # 0.77 where a flat bilinear patch gets 0.08.
    p00, p10, p01, p11 = e0[0], e0[-1], e1[0], e1[-1]
    B = ((1 - u) * (1 - v) * p00 + (1 - u) * v * p01
         + u * (1 - v) * p10 + u * v * p11)
    return F - B


def resample_face(Q, fu, fv, mode):
    """`linear` and `cubic` resample the FINE grid, which is what a refill of
    an existing mesh does. `chord` and `coons` model the GENERATIVE case of
    plan B, where no fine mesh exists: a model emits the block corners, and
    the only question is whether the edges between them carry geometry.

        chord   straight lines between the four block corners, bilinear inside
        coons   cubic edges through the true edge points, Coons patch inside
    """
    if mode == "linear":
        return tfi._lerp_axis(tfi._lerp_axis(Q, fu, 0), fv, 1)
    if mode == "cubic":
        return _cubic_axis(_cubic_axis(Q, fu, 0), fv, 1)
    n, m = len(fu), len(fv)
    if mode == "chord":
        c = np.array([Q[0, 0], Q[0, -1], Q[-1, 0], Q[-1, -1]])
        u = np.asarray(fu)[:, None, None]
        v = np.asarray(fv)[None, :, None]
        return ((1 - u) * (1 - v) * c[0] + (1 - u) * v * c[1]
                + u * (1 - v) * c[2] + u * v * c[3])
    # coons: the four edges resampled as cubics, interior by Coons blend
    e0 = _cubic_axis(Q[:, 0], fu, 0)      # u-edge at v=0
    e1 = _cubic_axis(Q[:, -1], fu, 0)     # u-edge at v=1
    e2 = _cubic_axis(Q[0, :], fv, 0)      # v-edge at u=0
    e3 = _cubic_axis(Q[-1, :], fv, 0)     # v-edge at u=1
    return _coons(e0, e1, e2, e3)


def boundary_face_grids(S, lat):
    """Every block face that lies on the domain boundary, as a point grid."""
    out = []
    for r, (_dims, vert) in lat.items():
        for ax in (0, 1, 2):
            for side in (0, 1):
                G = tfi.face_grid(vert, ax, side)
                if G.shape[0] < 2 or G.shape[1] < 2:
                    continue
                on = 0
                for i in range(G.shape[0] - 1):
                    for j in range(G.shape[1] - 1):
                        fk = frozenset((int(G[i, j]), int(G[i + 1, j]),
                                        int(G[i + 1, j + 1]), int(G[i, j + 1])))
                        if len(S.f2h.get(fk, ())) == 1:
                            on += 1
                if on:
                    out.append(S.P[G])
    return out


def measure(S, lab, grids, keep, mode, samples=3):
    """Resample every boundary face to `keep` of its divisions and report how
    far the result lands from the input surface.

    Points are taken INSIDE the new quads, not at their corners: corners of a
    coarse sampling still sit on original mesh vertices and would hide exactly
    the chordal error being looked for."""
    pts = []
    for Q in grids:
        n, m = Q.shape[0] - 1, Q.shape[1] - 1
        nn, mm = max(1, int(round(n * keep))), max(1, int(round(m * keep)))
        F = resample_face(Q, np.linspace(0, 1, nn + 1),
                          np.linspace(0, 1, mm + 1), mode)
        for a in np.linspace(0.2, 0.8, samples):
            for b in np.linspace(0.2, 0.8, samples):
                pts.append(((1 - a) * (1 - b) * F[:-1, :-1]
                            + a * (1 - b) * F[1:, :-1]
                            + a * b * F[1:, 1:]
                            + (1 - a) * b * F[:-1, 1:]).reshape(-1, 3))
    X = np.vstack(pts)
    d = np.linalg.norm(cb.project_to_surface(lab, X) - X, axis=1)
    return len(X), float(np.percentile(d, 95)), float(d.max())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("blocks")
    ap.add_argument("--input-vtk", required=True)
    ap.add_argument("--keep", type=float, nargs="+",
                    default=[1.0, 0.5, 0.25, 0.125, 0.0625])
    a = ap.parse_args()

    P, H, B, f2h = tfi.load_blocks(a.blocks)
    S, _bid, lab = cb.read_blocks_vtk(a.blocks, a.input_vtk)
    lat, _missing = tfi.lattices(P, H, B, f2h, verbose=False)
    grids = boundary_face_grids(S, lat)
    print(f"[curved] {len(lat)} blocks, {len(grids)} boundary face grids")
    print(f"\n{'divisions kept':>14} {'points':>8} "
          f"{'chord p95':>11} {'chord max':>11} "
          f"{'coons p95':>11} {'coons max':>11} {'gain p95':>9}")
    for keep in a.keep:
        nl, s95, smx = measure(S, lab, grids, keep, "chord")
        _n, k95, kmx = measure(S, lab, grids, keep, "coons")
        gain = s95 / k95 if k95 > 0 else float("inf")
        print(f"{keep * 100:13.1f}% {nl:8d} {s95:11.5f} {smx:11.5f} "
              f"{k95:11.5f} {kmx:11.5f} {gain:8.2f}x")


if __name__ == "__main__":
    main()
