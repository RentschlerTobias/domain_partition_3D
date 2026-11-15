"""Simple CAD envelope around the diffuser lip = mantle (gv_ru_shroudOut, axisymmetric, round nose at
its upstream tip) + ring 1 (dtIn / dtOut / bevelBack_dt) as ONE body.

outer skin: mantle meridian offset by delta outward (in the meridian plane, nose cut back where the
            curvature radius < R_TRIM delta, continued straight upstream to ext before the tip), revolved
            over the quarter; continued G1 by the offset of dtOut (same normal at z of the junction)
inner skin: offset of dtIn into the diffuser; upstream of it (inside the gv/ru hole, cut from the
            domain anyway) a surface of revolution leaving dtIn tangentially (G1) and turning axis-parallel
            to the front level, so it cannot drift into the mantle nose
ends:       ruled caps front / back, planar caps on the symmetry planes.
The interface planes between diffuser and gv/ru domain lose a rim strip of ~delta under the envelope.

run with LD_LIBRARY_PATH=/dtOO-install/lib /opt/venv-dtoo/bin/python lip_envelope.py <delta> [ext] [R_TRIM]
"""
import numpy as np
import simple_envelope as se
from simple_envelope import (faces, delta, ext, R_TRIM, offset_pair, aligned, close_and_check, BRep_Tool,
                             gp_Pnt, BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeEdge, brepfill)
from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnSurf, GeomAPI_PointsToBSpline
from OCC.Core.GeomLProp import GeomLProp_SLProps
from OCC.Core.TColgp import TColgp_Array1OfPnt
from OCC.Core.GeomAbs import GeomAbs_C2
from OCC.Core.Geom import Geom_SurfaceOfRevolution, Geom_BezierCurve
from OCC.Core.gp import gp_Ax1, gp_Dir

IN, OUT, BEVEL, MANTLE, TIP = 0, 3, 8, 11, 12          # face indices in cad_blk.stp (export order)


def mantle_meridian_offset(face, z_front, n=400):
    """(r, z) polyline of the mantle meridian offset by delta away from the axis, nose cut back to the
    first point with curvature radius >= R_TRIM delta, plus a straight upstream continuation to z_front"""
    srf = BRep_Tool.Surface(face)
    u, v = GeomAPI_ProjectPointOnSurf(gp_Pnt(2.15, 0, -1.0), srf).LowerDistanceParameters()
    u0, u1, v0, v1 = srf.Bounds()
    rows = []
    for t in np.linspace(0, 1, n) ** 1.5:                  # denser at the nose (low v)
        q = v0 + t * (v1 - v0)
        pr = GeomLProp_SLProps(srf, u, q, 2, 1e-9)
        p = np.array(pr.Value().Coord()); nn = np.array(pr.Normal().Coord())
        r = np.hypot(p[0], p[1]); er = p[:2] / r
        n2 = np.array([nn[:2] @ er, nn[2]]); n2 /= np.linalg.norm(n2)
        k = max(abs(pr.MinCurvature()), abs(pr.MaxCurvature())) if pr.IsCurvatureDefined() else 0
        rows.append((r, p[2], *n2, 1 / k if k > 0 else 1e9))
    A = np.array(rows)
    if A[len(A) // 2, 2] < 0:                              # normal away from the axis
        A[:, 2:4] *= -1
    if A[0, 1] > A[-1, 1]:                                 # upstream first
        A = A[::-1]
    i0 = int(np.argmax(A[:, 4] >= R_TRIM * delta))
    print(f"   mantle nose cut back to z {A[i0, 1]:.3f} (curvature radius < {R_TRIM * delta:.3f} before)")
    Q = A[i0:, :2] + delta * A[i0:, 2:4]
    t = Q[0] - Q[1]; t /= np.linalg.norm(t)                # upstream tangent
    L = (Q[0, 1] - z_front) / -t[1]
    m = max(4, int(L / np.linalg.norm(Q[1] - Q[0])))
    pre = [Q[0] + t * L * k / m for k in range(m, 0, -1)]
    return np.concatenate([np.array(pre), Q]), A


CIRC = 0.15        # skins reach this far across the symmetry planes, trimmed to the quarter afterwards


def revolve_quarter(RZ, bezier=False):
    pts = TColgp_Array1OfPnt(1, len(RZ))
    for k, (r, z) in enumerate(RZ):
        pts.SetValue(k + 1, gp_Pnt(r, 0, z))
    crv = Geom_BezierCurve(pts) if bezier else GeomAPI_PointsToBSpline(pts, 3, 3, GeomAbs_C2, 1e-6).Curve()
    rev = Geom_SurfaceOfRevolution(crv, gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)))
    c0, c1 = crv.FirstParameter(), crv.LastParameter()
    a = CIRC / RZ[:, 0].min()
    face = BRepBuilderAPI_MakeFace(rev, -a, np.pi / 2 + a, c0, c1, 1e-6).Face()
    return (face, BRepBuilderAPI_MakeEdge(rev.VIso(c0), -a, np.pi / 2 + a).Edge(),
            BRepBuilderAPI_MakeEdge(rev.VIso(c1), -a, np.pi / 2 + a).Edge())


def edge_rz(edge, n=41):
    """(r, z) samples of an edge inside the quarter"""
    c, f0, f1 = BRep_Tool.Curve(edge)
    P = [c.Value(f0 + s * (f1 - f0)) for s in np.linspace(0, 1, n)]
    return np.array([[np.hypot(p.X(), p.Y()), p.Z()] for p in P if min(p.X(), p.Y()) >= 0])


def main():
    z_tip = min(BRep_Tool.Surface(faces[TIP]).Value(*uv).Z() for uv in
                [(a, b) for a in np.linspace(*BRep_Tool.Surface(faces[TIP]).Bounds()[:2], 9)
                 for b in np.linspace(*BRep_Tool.Surface(faces[TIP]).Bounds()[2:], 9)])
    z_front = z_tip - ext
    # ring 1: back end like rings 2/3; upstream: dtIn straight on to z_front, dtOut attached to the mantle
    srf_in = BRep_Tool.Surface(faces[IN]); u0, u1, v0, v1 = srf_in.Bounds()
    z_in0 = min(srf_in.Value(u, v).Z() for u in (u0, u1) for v in (v0, v1))
    (Fin, fin_front, fin_back), (Fout, fout_front, fout_back) = offset_pair(
        faces[IN], faces[OUT], None, front=(None, None), circ_ext=CIRC)
    # inner skin in the hole: cubic meridian, tangent to the dtIn offset at its upstream edge, axis-parallel
    # at the front level
    J = edge_rz(fin_front); r0, z0 = J.mean(0)
    G, _ = se.skin_grid(faces[IN])
    d = np.array([np.hypot(*G[0, :, :2].T).mean() - np.hypot(*G[3, :, :2].T).mean(), G[0, :, 2].mean() - G[3, :, 2].mean()])
    d /= np.linalg.norm(d)                                  # upstream meridian tangent of dtIn
    h = (z0 - z_front) / 3
    P1 = np.array([r0, z0]) + d * h / max(-d[1], 1e-3)
    Fhole, fh_front, fh_back = revolve_quarter(np.array([[r0, z0], P1, [P1[0], z_front + h], [P1[0], z_front]])[::-1], bezier=True)
    print(f"   inner skin in the hole: r {r0:.4f} at z {z0:.3f} -> r {P1[0]:.4f} at z {z_front:.3f} "
          f"(dtIn offset edge roundness {np.ptp(J[:, 0]):.1e})")
    RZ, A = mantle_meridian_offset(faces[MANTLE], z_front)
    # junction mantle offset / dtOut offset (same normal there): end the meridian exactly on the dtOut
    # offset edge so both faces share it
    c, f0, f1 = BRep_Tool.Curve(fout_front)
    J = np.array([[np.hypot(*c.Value(f0 + s * (f1 - f0)).Coord()[:2]), c.Value(f0 + s * (f1 - f0)).Z()]
                  for s in np.linspace(0, 1, 41)
                  if min(c.Value(f0 + s * (f1 - f0)).X(), c.Value(f0 + s * (f1 - f0)).Y()) >= 0])
    print(f"   junction mantle/dtOut offsets: mismatch r {abs(J[:, 0].mean() - RZ[-1, 0]):.1e} "
          f"z {abs(J[:, 1].mean() - RZ[-1, 1]):.1e}, dtOut edge roundness {np.ptp(J[:, 0]):.1e} / {np.ptp(J[:, 1]):.1e}")
    RZ[-1] = J.mean(0)
    Fm, fm_front, _ = revolve_quarter(RZ)
    cap_front = brepfill.Face(fh_front, aligned(fm_front, fh_front))
    cap_back = brepfill.Face(fin_back, aligned(fout_back, fin_back))
    quarter = lambda p: p.X() > 1e-6 and p.Y() > 1e-6
    close_and_check("lip", (Fin, Fhole, Fout, Fm, cap_front, cap_back), (IN, OUT, BEVEL, MANTLE, TIP), quarter, tol=1e-3, quarter=True)


if __name__ == "__main__":
    main()
