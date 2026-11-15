"""Simple CAD envelope around a thin ring: its inner and outer skin (single smooth dtOO faces) extended
streamwise beyond the ring tips and offset by delta away from the ring; the ends closed with ruled
surfaces between the two offset skins, the symmetry planes with planar caps. Bevels / knife edges of
the ring play no role, the envelope has 6 faces and ~90 deg edges.

run with LD_LIBRARY_PATH=/dtOO-install/lib /opt/venv-dtoo/bin/python simple_envelope.py <delta> [ext]
"""
import sys
import numpy as np
from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_EDGE, TopAbs_SHELL, TopAbs_SOLID, TopAbs_WIRE
from OCC.Core.TopoDS import topods
from OCC.Core.BRep import BRep_Tool
from OCC.Core.Geom import Geom_OffsetSurface, Geom_BSplineSurface, Geom_Plane
from OCC.Core.BRepLib import BRepLib_FindSurface
from OCC.Core.GeomLib import geomlib
from OCC.Core.GeomConvert import geomconvert
from OCC.Core.BRepBuilderAPI import (BRepBuilderAPI_MakeFace, BRepBuilderAPI_Sewing, BRepBuilderAPI_MakeSolid,
                                     BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeEdge)
from OCC.Core.BRepFill import brepfill
from OCC.Core.ShapeAnalysis import ShapeAnalysis_FreeBounds
from OCC.Core.ShapeFix import ShapeFix_Solid, ShapeFix_Shape
from OCC.Core.BRepCheck import BRepCheck_Analyzer
from OCC.Core.GProp import GProp_GProps
from OCC.Core.BRepGProp import brepgprop
from OCC.Core.BRepTools import breptools
from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnSurf
from OCC.Core.BRepClass3d import BRepClass3d_SolidClassifier
from OCC.Core.TopAbs import TopAbs_IN
from OCC.Core.gp import gp_Pnt, gp_Vec
from OCC.Core.BRepExtrema import BRepExtrema_DistShapeShape
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
from OCC.Core.TColStd import TColStd_Array1OfReal, TColStd_Array1OfInteger, TColStd_Array2OfReal
from OCC.Core.GeomLProp import GeomLProp_SLProps
from OCC.Core.TColgp import TColgp_Array2OfPnt
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface
from OCC.Core.GeomAbs import GeomAbs_C2, GeomAbs_C1
from OCC.Core.GeomConvert import GeomConvert_ApproxSurface
from OCC.Core.Approx import Approx_ChordLength

delta = float(sys.argv[1]) if len(sys.argv) > 1 else 0.0354
ext = float(sys.argv[2]) if len(sys.argv) > 2 else 2 * delta
R_TRIM = float(sys.argv[3]) if len(sys.argv) > 3 else 4.0   # trim skin ends with curvature radius < R_TRIM * delta

r = STEPControl_Reader(); r.ReadFile("/work/canada/case/cad_blk.stp"); r.TransferRoots(); shape = r.OneShape()
faces = []
ex = TopExp_Explorer(shape, TopAbs_FACE)
while ex.More():
    faces.append(topods.Face(ex.Current())); ex.Next()
RINGS = {"ring2": (1, 4, [1, 4, 6, 7]), "ring3": (2, 5, [2, 5, 9, 10])}   # (inner skin, outer skin, all faces)


def volume(s):
    g = GProp_GProps(); brepgprop.VolumeProperties(s, g); return g.Mass()


def count(s, t):
    e = TopExp_Explorer(s, t); n = 0
    while e.More():
        n += 1; e.Next()
    return n


def ring_solid(idx):
    sew = BRepBuilderAPI_Sewing(1e-5)
    for i in idx:
        sew.Add(faces[i])
    sew.Perform()
    sh = sew.SewedShape()
    fb = ShapeAnalysis_FreeBounds(sh, 1e-6, False, False)
    s2 = BRepBuilderAPI_Sewing(1e-5); s2.Add(sh)
    e = TopExp_Explorer(fb.GetClosedWires(), TopAbs_WIRE)
    while e.More():
        s2.Add(BRepBuilderAPI_MakeFace(topods.Wire(e.Current()), True).Face()); e.Next()
    s2.Perform()
    ms = BRepBuilderAPI_MakeSolid()
    e = TopExp_Explorer(s2.SewedShape(), TopAbs_SHELL)
    while e.More():
        ms.Add(topods.Shell(e.Current())); e.Next()
    f = ShapeFix_Solid(ms.Solid()); f.Perform(); return f.Solid()


def skin_grid(face, ns=41, nc=61):
    """points of a skin on a (streamwise x circumferential) grid, circumferential index by increasing phi"""
    srf = BRep_Tool.Surface(face); u0, u1, v0, v1 = srf.Bounds()
    um, vm = 0.5 * (u0 + u1), 0.5 * (v0 + v1)
    dz_u = abs(srf.Value(u1, vm).Z() - srf.Value(u0, vm).Z())
    dz_v = abs(srf.Value(um, v1).Z() - srf.Value(um, v0).Z())
    U, V = np.linspace(u0, u1, ns if dz_u > dz_v else nc), np.linspace(v0, v1, nc if dz_u > dz_v else ns)
    G = np.array([[[*srf.Value(u, v).Coord()] for v in V] for u in U])
    N = np.zeros_like(G)
    for i, u in enumerate(U):
        for j, v in enumerate(V):
            pr = GeomLProp_SLProps(srf, u, v, 1, 1e-9)
            N[i, j] = [*pr.Normal().Coord()] if pr.IsNormalDefined() else np.nan
    if dz_u <= dz_v:
        G, N = G.transpose(1, 0, 2), N.transpose(1, 0, 2)
    if G[-1, :, 2].mean() < G[0, :, 2].mean():          # streamwise index by increasing z
        G, N = G[::-1], N[::-1]
    phi = np.arctan2(G[..., 1], G[..., 0])
    if phi[:, -1].mean() < phi[:, 0].mean():
        G, N = G[:, ::-1], N[:, ::-1]
    return G, N


def _poles(srf):
    nu, nv = srf.NbUPoles(), srf.NbVPoles()
    P = np.array([[[*srf.Pole(i, j).Coord(), 1.0] for j in range(1, nv + 1)] for i in range(1, nu + 1)])
    for i in range(nu):
        for j in range(nv):
            w = srf.Weight(i + 1, j + 1); P[i, j, :3] *= w; P[i, j, 3] = w
    return P


def _knots(srf, in_u):
    n = srf.NbUKnots() if in_u else srf.NbVKnots()
    return ([srf.UKnot(k) if in_u else srf.VKnot(k) for k in range(1, n + 1)],
            [srf.UMultiplicity(k) if in_u else srf.VMultiplicity(k) for k in range(1, n + 1)])


def _arr(vals, typ):
    A = typ(1, len(vals))
    for k, v in enumerate(vals):
        A.SetValue(k + 1, v)
    return A


def extend_linear(srf, in_u, L_lo, L_hi):
    """continue a (non-periodic) B-spline surface beyond both ends of one parameter direction by
    straight tangent extensions, chord length per pole column from L_lo / L_hi = (phi, length) samples,
    G1 at the junction (tangent direction kept, length per column free). Implemented as an
    extra span with poles on the tangent lines (pythonocc's GeomLib.ExtendSurfByLength loses its result).
    Rational surfaces: the span is linear in homogeneous coordinates, weights held constant."""
    deg = srf.UDegree() if in_u else srf.VDegree()
    P = _poles(srf)
    if not in_u:
        P = P.transpose(1, 0, 2)                    # first index = extension direction
    P0 = P.copy()
    K, M = _knots(srf, in_u)
    a, b = K[0], K[-1]
    # end derivative of the homogeneous curve per column: d (P_end - P_prev) / (span param)
    # derivative at the ends of a clamped B-spline: deg / (t_{n+deg} - t_n) * (P_n - P_{n-1})
    flat = np.repeat(K, M)
    n = len(P) - 1
    D_hi = deg / (flat[n + deg] - flat[n]) * (P[-1] - P[-2])
    D_lo = deg / (flat[deg + 1] - flat[1]) * (P[1] - P[0])
    # per pole column chord length, interpolated over the circumferential angle of the end pole;
    # unit parameter span for the extension, so the junction is G1 (tangent planes match)
    def col_len(Lphi, Pend, D):
        phi_p = np.arctan2(Pend[:, 1] / Pend[:, 3], Pend[:, 0] / Pend[:, 3])
        o = np.argsort(Lphi[0])
        L = np.interp(phi_p, Lphi[0][o], Lphi[1][o])
        return (L / (np.linalg.norm(D[:, :3], axis=1) / Pend[:, 3]))[:, None]
    e_par = 0.05 * (b - a)                              # parameter length of the extension spans
    if L_hi is not None:
        e = col_len(L_hi, P[-1], D_hi)
        P = np.concatenate([P, np.array([P[-1] + D_hi * e * k / deg for k in range(1, deg + 1)])])
        K = K + [b + e_par]; M = M[:-1] + [deg, deg + 1]
    if L_lo is not None:
        e = col_len(L_lo, P[0], D_lo)
        P = np.concatenate([np.array([P[0] - D_lo * e * k / deg for k in range(deg, 0, -1)]), P])
        K = [a - e_par] + K; M = [deg + 1, deg] + M[1:]
    # boundary columns lying in a symmetry plane (x = 0 or y = 0) keep their extension in that plane
    for j in (0, -1):
        for c in (0, 1):
            if np.all(np.abs(P0[:, j, c] / P0[:, j, 3]) < 1e-6):
                P[:, j, c] = 0.0
    if not in_u:
        P = P.transpose(1, 0, 2)
    nu, nv = P.shape[:2]
    poles = TColgp_Array2OfPnt(1, nu, 1, nv); wts = TColStd_Array2OfReal(1, nu, 1, nv)
    for i in range(nu):
        for j in range(nv):
            w = P[i, j, 3]; poles.SetValue(i + 1, j + 1, gp_Pnt(*(P[i, j, :3] / w))); wts.SetValue(i + 1, j + 1, w)
    uk, um = _knots(srf, True); vk, vm = _knots(srf, False)
    if in_u: uk, um = K, M
    else: vk, vm = K, M
    return Geom_BSplineSurface(poles, wts, _arr(uk, TColStd_Array1OfReal), _arr(vk, TColStd_Array1OfReal),
                               _arr(um, TColStd_Array1OfInteger), _arr(vm, TColStd_Array1OfInteger),
                               srf.UDegree(), srf.VDegree(), False, False)


def min_radius_row(srf, in_u, par, n=201):
    o = srf.Bounds(); other = np.linspace(o[2], o[3], n) if in_u else np.linspace(o[0], o[1], n)
    r = 1e9
    for q in other:
        pr = GeomLProp_SLProps(srf, par, q, 2, 1e-9) if in_u else GeomLProp_SLProps(srf, q, par, 2, 1e-9)
        if pr.IsCurvatureDefined():
            k = max(abs(pr.MinCurvature()), abs(pr.MaxCurvature()))
            if k > 0:
                r = min(r, 1 / k)
    return r


def trim_tight_ends(srf, in_u, r_min, ends=(True, True)):
    """cut each end of the streamwise direction back to where the surface curvature radius is >= r_min
    (an offset by delta folds where the radius on its concave side is < delta)"""
    u0, u1, v0, v1 = srf.Bounds()
    a, b = (u0, u1) if in_u else (v0, v1)
    cut = []
    for end, inward, on in ((a, 1, ends[0]), (b, -1, ends[1])):
        p = end
        while on and min_radius_row(srf, in_u, p) < r_min and abs(p - end) < 0.3 * (b - a):
            p += inward * 0.005 * (b - a)
        cut.append(p)
    srf = Geom_BSplineSurface.DownCast(srf.Copy())
    if cut != [a, b]:
        if in_u: srf.Segment(cut[0], cut[1], v0, v1)
        else: srf.Segment(u0, u1, cut[0], cut[1])
    return srf, (cut[0] - a) / (b - a), (b - cut[1]) / (b - a)


def offset_pair(face_a, face_b, solid=None, front=None, circ_ext=0.0):
    """both skins: the dtOO surface extended streamwise (C2, GeomLib) at each end far enough to reach ext
    beyond the farther ring tip (oblique bevels put the tip of one skin past the end of the other), then
    offset by delta away from the ring (exact offset surface of the smooth skin, no fitting).
    Offset side: outside `solid`, or (solid None) away from the opposite skin.
    front: None = upstream ends like the downstream ends; else (L_a, L_b), per skin a fixed straight
    extension length or None for no extension (and no trimming) at the upstream end.
    circ_ext: straight extension across the symmetry-plane edges (then trim with close_and_check(quarter=True);
    needed where the skin does not meet the symmetry planes at right angles).
    Returns per skin (face, front edge, back edge)."""
    grids = []
    for face, other in ((face_a, face_b), (face_b, face_a)):
        G, N = skin_grid(face)
        mid = G[len(G) // 2, G.shape[1] // 2]; nm = N[len(G) // 2, G.shape[1] // 2]
        if solid is not None:
            sgn = 1.0 if BRepClass3d_SolidClassifier(solid, gp_Pnt(*(mid + 0.5 * delta * nm)), 1e-7).State() != TopAbs_IN else -1.0
        else:
            dist = lambda p: BRepExtrema_DistShapeShape(BRepBuilderAPI_MakeVertex(gp_Pnt(*p)).Vertex(), other).Value()
            sgn = 1.0 if dist(mid + 0.5 * delta * nm) > dist(mid - 0.5 * delta * nm) else -1.0
        grids.append((G, G + sgn * delta * N, sgn))
    unit = lambda v: v / np.linalg.norm(v, axis=-1)[..., None]
    out = []
    for si, (face, (G, O, sgn)) in enumerate(zip((face_a, face_b), grids)):
        srf = geomconvert.SurfaceToBSplineSurface(BRep_Tool.Surface(face))
        u0, u1, v0, v1 = srf.Bounds(); um, vm = 0.5 * (u0 + u1), 0.5 * (v0 + v1)
        in_u = abs(srf.Value(u1, vm).Z() - srf.Value(u0, vm).Z()) > abs(srf.Value(um, v1).Z() - srf.Value(um, v0).Z())
        # which parameter end is upstream (low z)
        lo_first = (srf.Value(u0, vm).Z() < srf.Value(u1, vm).Z()) if in_u else (srf.Value(um, v0).Z() < srf.Value(um, v1).Z())
        Ls = {}
        for end, nxt in ((0, 1), (-1, -2)):
            ts = unit(O[end] - O[nxt])
            beyond = np.max([((G_[end] - O[end]) * ts).sum(1) for G_, _, _ in grids], axis=0)
            # both skins should end on a common level along their mean tangent t (end cap ~normal to them)
            t = unit(sum(unit(O_[end] - O_[nxt]) for _, O_, _ in grids))
            level = np.max([(P[end] * t).sum(1) for G_, O_, _ in grids for P in (G_, O_)], axis=0) + ext
            match = (level - (O[end] * t).sum(1)) / np.maximum((ts * t).sum(1), 0.3)
            need = ext + np.maximum(0.0, beyond)                        # per column, never short of a tip
            L = np.maximum(match, need)
            phi = np.arctan2(O[end][:, 1], O[end][:, 0])
            # smooth along the circumference (abrupt length changes make neighbouring extension lines
            # cross); upward-only so containment holds
            w = np.exp(-0.5 * ((phi[:, None] - phi[None, :]) / np.radians(8)) ** 2); w /= w.sum(1)[:, None]
            for _ in range(20):
                L = np.maximum(w @ L, need)
            Ls[end] = (phi, L)
        if front is not None:                                 # fixed upstream end
            Ls[0] = None if front[si] is None else (Ls[0][0], np.full_like(Ls[0][1], front[si]))
        trim = (True, True)
        if front is not None and front[si] is None:
            trim = (not lo_first, lo_first)                 # never trim the attached upstream end
        if circ_ext:
            # edges reaching across the symmetry planes: cut back tiny folds / tight curvature at the
            # symmetry-plane edges as well (the straight extension starts from the edge tangent)
            srf, c_lo, c_hi = trim_tight_ends(srf, not in_u, R_TRIM * delta)
            if c_lo or c_hi:
                print(f"   skin trimmed at the symmetry-plane edges by {c_lo:.3f} / {c_hi:.3f} of its width")
        srf, f_lo, f_hi = trim_tight_ends(srf, in_u, R_TRIM * delta, trim)
        # trimmed-off part of the ring end has to be covered by the extension too (arc length of the cut)
        u0, u1, v0, v1 = srf.Bounds()
        span = np.linalg.norm(G[-1] - G[0], axis=1).mean()
        # grid end 0 is upstream (low z); it is the parameter min end iff lo_first
        add = lambda pl, f: None if pl is None else (pl[0], pl[1] + f * span)
        L_lo = add(Ls[0] if lo_first else Ls[-1], f_lo)
        L_hi = add(Ls[-1] if lo_first else Ls[0], f_hi)
        if f_lo or f_hi:
            print(f"   skin trimmed at the ends by {f_lo:.3f} / {f_hi:.3f} of its length (curvature radius < {R_TRIM * delta:.3f})")
        # offset the trimmed skin (radius >= 1.5 delta: no fold), as B-spline, then extend it straight;
        # extending before offsetting folds where the extension lines converge
        ap = GeomConvert_ApproxSurface(Geom_OffsetSurface(srf, sgn * delta, True), 1e-6, GeomAbs_C1, GeomAbs_C1,
                                       5, 5, 200, 1)
        off = extend_linear(ap.Surface(), in_u, L_lo, L_hi)
        if circ_ext:
            const = (np.array([-1e3, 1e3]), np.array([circ_ext, circ_ext]))
            off = extend_linear(off, not in_u, const, const)
        u0, u1, v0, v1 = off.Bounds()
        isos = [off.UIso(u0), off.UIso(u1)] if in_u else [off.VIso(v0), off.VIso(v1)]
        if not lo_first:
            isos = isos[::-1]
        out.append((BRepBuilderAPI_MakeFace(off, 1e-6).Face(),
                    BRepBuilderAPI_MakeEdge(isos[0]).Edge(), BRepBuilderAPI_MakeEdge(isos[1]).Edge()))
    return out


def aligned(e, ref):
    """edge e with the same running direction as ref (ruled caps must not twist)"""
    c, f0, f1 = BRep_Tool.Curve(e); r, g0, g1 = BRep_Tool.Curve(ref)
    if c.Value(f0).Distance(r.Value(g0)) > c.Value(f1).Distance(r.Value(g0)):
        return BRepBuilderAPI_MakeEdge(c.Reversed()).Edge()
    return e


def close_and_check(name, skin_faces, check_faces, check_filter=lambda p: True, tol=1e-4, quarter=False):
    """sew skins + end caps, close the symmetry planes with planar caps, solid, containment check.
    quarter: skins reach across the symmetry planes and are trimmed to x >= 0, y >= 0 first"""
    sew = BRepBuilderAPI_Sewing(tol)
    box = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, -1e3), gp_Pnt(1e3, 1e3, 1e3)).Shape()
    for f in skin_faces:
        if quarter:
            ex_ = TopExp_Explorer(BRepAlgoAPI_Common(f, box).Shape(), TopAbs_FACE)
            while ex_.More():
                sew.Add(ex_.Current()); ex_.Next()
        else:
            sew.Add(f)
    sew.Perform()
    sh = sew.SewedShape()
    fb = ShapeAnalysis_FreeBounds(sh, 0.1 * tol, False, False)
    s2 = BRepBuilderAPI_Sewing(tol); s2.Add(sh)
    nc = 0
    e = TopExp_Explorer(fb.GetClosedWires(), TopAbs_WIRE)
    while e.More():
        w = topods.Wire(e.Current())
        fs = BRepLib_FindSurface(w, tol, True)          # symmetry-plane wire, planar within 1e-4
        if fs.Found():
            mf = BRepBuilderAPI_MakeFace(Geom_Plane.DownCast(fs.Surface()).Pln(), w, True)
            if mf.IsDone():
                s2.Add(mf.Face()); nc += 1
        e.Next()
    s2.Perform()
    print(f"   {name}: free wires {count(fb.GetClosedWires(), TopAbs_WIRE)} -> caps {nc}, "
          f"shells {count(s2.SewedShape(), TopAbs_SHELL)}", flush=True)
    ms = BRepBuilderAPI_MakeSolid()
    e = TopExp_Explorer(s2.SewedShape(), TopAbs_SHELL)
    while e.More():
        ms.Add(topods.Shell(e.Current())); e.Next()
    fix = ShapeFix_Solid(ms.Solid()); fix.Perform(); env = fix.Solid()
    if not BRepCheck_Analyzer(env).IsValid():
        sf = ShapeFix_Shape(env); sf.SetPrecision(1e-6); sf.Perform()
        e = TopExp_Explorer(sf.Shape(), TopAbs_SOLID)
        env = topods.Solid(e.Current())
    inside = True
    for i in check_faces:
        srf = BRep_Tool.Surface(faces[i]); u0, u1, v0, v1 = srf.Bounds()
        for uu in np.linspace(u0, u1, 25):
            for vv in np.linspace(v0, v1, 25):
                p = srf.Value(uu, vv)
                if min(p.X(), p.Y()) < 1e-6 or not check_filter(p):   # symmetry planes / outside the quarter
                    continue
                st = BRepClass3d_SolidClassifier(env, p, 1e-7).State()
                if st != TopAbs_IN:
                    inside = False; print(f"   not inside: face {i} state {st} at {p.X():.3f} {p.Y():.3f} {p.Z():.3f}")
    print(f"[env] {name}: faces {count(env, TopAbs_FACE)} (caps {nc}), edges {count(env, TopAbs_EDGE)}, "
          f"valid {BRepCheck_Analyzer(env).IsValid()}, volume {volume(env):.4f}, wall strictly inside {inside}", flush=True)
    breptools.Write(env, f"/work/canada/cad/{name}_simple_envelope.brep")
    return env


def main():
    for name, (fi, fo, allf) in RINGS.items():
        solid = ring_solid(allf)
        (Fin, fin_front, fin_back), (Fout, fout_front, fout_back) = offset_pair(faces[fi], faces[fo], solid)
        cap_front = brepfill.Face(fin_front, aligned(fout_front, fin_front))
        cap_back = brepfill.Face(fin_back, aligned(fout_back, fin_back))
        close_and_check(name, (Fin, Fout, cap_front, cap_back), allf)


if __name__ == "__main__":
    main()
