"""Draft tube (blk domain, quarter) with the thin walls as zero-thickness mid-surfaces ("baffles").
Stage A (pythonOCC): writes the CAD parts for the gmsh assembly (assemble_mesh.py) into parts/.

  - thin-wall skin pairs are detected generically among the curved CAD faces: two faces whose mutual
    distance is small and nearly constant compared with their size (no face indices / labels used);
  - mid-surface of a pair: points of skin A projected onto skin B, midpoints, fitted B-spline surface.
    The grid is mirrored across the symmetry planes (x=0, y=0) before fitting, so the mid-surface crosses
    them smoothly; the quarter is cut out later by the boolean fragment;
  - the pair attached to the gv/ru shroud (the lip = mantle + ring 1) is continued upstream along the
    mid-line between the shroud (straight line between the outer rims of the two interface faces) and the
    mantle, so the hub-chain obstacle (bulb + gv/ru hole + dt hub) is closed by the lip mid-surface and the
    ring-1 part of the lip continues it as an interior baffle;
  - outer body (far field + inlet/outlet + symmetry planes) from domain_blk.outer_body(); hub-chain faces and
    the two interface faces (extended radially beyond the lip mid-surface, trimmed in the gmsh stage).
Hub-chain faces are taken from domain_blk.F (same face naming as the cad1 domain).

run: LD_LIBRARY_PATH=/dtOO-install/lib /opt/venv-dtoo/bin/python build_parts.py
"""
import sys, os, json
sys.path.insert(0, "/work/canada/cad")
import numpy as np
import domain_blk as db
from simple_envelope import faces, BRep_Tool, skin_grid
from OCC.Core.gp import gp_Pnt, gp_Ax1, gp_Dir
from OCC.Core.TColgp import TColgp_Array2OfPnt, TColgp_Array1OfPnt
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface, GeomAPI_ProjectPointOnSurf, GeomAPI_PointsToBSpline
from OCC.Core.GeomAbs import GeomAbs_C2, GeomAbs_Plane
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeVertex, BRepBuilderAPI_MakeEdge
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeRevol
from OCC.Core.BRepExtrema import BRepExtrema_DistShapeShape
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.BRepTools import breptools
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib

OUT = "/work/canada/blocking_study/dt_baffle/parts"
os.makedirs(OUT, exist_ok=True)
info = {}


def bbox(s):
    b = Bnd_Box(); brepbndlib.Add(s, b); return np.array(b.Get())


def dist_pt(p, shape):
    return BRepExtrema_DistShapeShape(BRepBuilderAPI_MakeVertex(gp_Pnt(*p)).Vertex(), shape).Value()


def face_samples(f, n=7):
    u0, u1, v0, v1 = breptools.UVBounds(f); s = BRep_Tool.Surface(f)
    return np.array([[*s.Value(u, v).Coord()] for u in np.linspace(u0, u1, n + 2)[1:-1]
                     for v in np.linspace(v0, v1, n + 2)[1:-1]])


# ---------------------------------------------------------------- 1. thin-wall pairs (generic)
def thin_wall_pairs(cands):
    """pairs (a, b): every interior sample of a lies within a near-constant small distance of b and vice versa"""
    S = {i: face_samples(faces[i]) for i in cands}
    B = {i: bbox(faces[i]) for i in cands}
    diag = {i: np.linalg.norm(B[i][3:] - B[i][:3]) for i in cands}
    pairs = []
    for ia, a in enumerate(cands):
        for b in cands[ia + 1:]:
            lo = np.maximum(B[a][:3], B[b][:3]); hi = np.minimum(B[a][3:], B[b][3:])
            if np.any(lo > hi + 0.2):
                continue
            dab = np.array([dist_pt(p, faces[b]) for p in S[a]])
            dba = np.array([dist_pt(p, faces[a]) for p in S[b]])
            d = np.r_[dab, dba]; size = min(diag[a], diag[b])
            ok = np.median(d) > 1e-6 and np.median(d) < 0.05 * size and np.percentile(d, 95) < 0.1 * size
            if ok:
                pairs.append((a, b, float(np.median(d)), float(d.max())))
    return pairs


curved = [i for i, f in enumerate(faces) if BRepAdaptor_Surface(f).GetType() != GeomAbs_Plane]
pairs = thin_wall_pairs(curved)
print("[parts] thin-wall pairs (face a, face b, median gap, max gap):", [(a, b, round(m, 4), round(x, 4)) for a, b, m, x in pairs])


# ---------------------------------------------------------------- 2. mid-surfaces
def mid_grid(a, b, ns=41, nc=61):
    G, _ = skin_grid(faces[a], ns, nc)                      # streamwise (increasing z) x phi (increasing)
    sb = BRep_Tool.Surface(faces[b])
    M = np.zeros_like(G)
    for i in range(G.shape[0]):
        for j in range(G.shape[1]):
            pr = GeomAPI_ProjectPointOnSurf(gp_Pnt(*G[i, j]), sb)
            q = np.array(pr.NearestPoint().Coord())
            M[i, j] = 0.5 * (G[i, j] + q)
    # exact symmetry-plane columns (projection may drift by the fit tolerance)
    M[:, 0, 1] = 0.0; M[:, -1, 0] = 0.0
    return M


def resample_columns(M, ns=None):
    """each streamwise column resampled to uniform arc length (projection clamping at knife edges bunches rows)"""
    ns = ns or M.shape[0]
    out = np.zeros((ns, M.shape[1], 3))
    for j in range(M.shape[1]):
        c = M[:, j]; s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(c, axis=0), axis=1))]
        t = np.linspace(0, s[-1], ns)
        out[:, j] = np.c_[[np.interp(t, s, c[:, k]) for k in range(3)]].T
    return out


def mirror_ext(M, k=6):
    """add k mirrored columns beyond y=0 (first column) and beyond x=0 (last column)"""
    lo = M[:, 1:k + 1][:, ::-1] * np.array([1, -1, 1])
    hi = M[:, -k - 1:-1][:, ::-1] * np.array([-1, 1, 1])
    return np.concatenate([lo, M, hi], axis=1)


def fit_surface(M, name):
    ns, nc = M.shape[:2]
    arr = TColgp_Array2OfPnt(1, ns, 1, nc)
    for i in range(ns):
        for j in range(nc):
            arr.SetValue(i + 1, j + 1, gp_Pnt(*M[i, j]))
    srf = GeomAPI_PointsToBSplineSurface(arr, 3, 8, GeomAbs_C2, 1e-5).Surface()
    face = BRepBuilderAPI_MakeFace(srf, 1e-6).Face()
    dev = max(GeomAPI_ProjectPointOnSurf(gp_Pnt(*M[i, j]), srf).LowerDistance()
              for i in range(0, ns, 2) for j in range(0, nc, 2))
    print(f"   {name}: grid {ns}x{nc}, fit max deviation {dev:.2e}, poles {srf.NbUPoles()}x{srf.NbVPoles()}")
    return face, dev


def meridian_pts(face, n=60):
    """(r, z) samples of a surface-of-revolution face, from a dense parametric grid (axisymmetric)"""
    s = BRep_Tool.Surface(face); u0, u1, v0, v1 = breptools.UVBounds(face)
    P = np.array([[*s.Value(u, v).Coord()] for u in np.linspace(u0, u1, 40) for v in np.linspace(v0, v1, 400)])
    RZ = np.c_[np.hypot(P[:, 0], P[:, 1]), P[:, 2]]
    zb = np.linspace(RZ[:, 1].min(), RZ[:, 1].max(), n)
    r = np.array([RZ[np.argsort(np.abs(RZ[:, 1] - z))[:40], 0].mean() for z in zb])
    return np.c_[r, zb]


def rim(face):
    """outer (max r) end point (r, z) of the meridian of an interface face"""
    m = db.meridian(face); pts = []
    for e in m:
        p, q = db.ends(e); pts += [p, q]
    p = max(pts, key=lambda p: p[0])
    return np.array([p[0], p[2]])


mid_faces = {}
F = db.F
mantle = faces[11]          # gv/ru mantle (outer skin of the lip upstream of the dt interface), as in lip_envelope.py
for a, b, gap, gmax in pairs:
    # orient: a = skin nearer to the axis (inner)
    ra = np.hypot(*face_samples(faces[a])[:, :2].T).mean(); rb = np.hypot(*face_samples(faces[b])[:, :2].T).mean()
    if ra > rb:
        a, b = b, a
    M = mid_grid(a, b)
    is_lip = min(dist_pt(p, mantle) for p in M[0, ::10]) < 0.2     # upstream end of this wall sits on the mantle
    name = f"mid_{a}_{b}"
    if is_lip:
        name = "mid_lip"
        B_ = rim(faces[F["if_gv"]]); C_ = rim(faces[F["if_dt"]])     # shroud line B -> C
        RZm = meridian_pts(mantle)
        z0 = B_[1] - 0.03                                       # a little upstream of the interface (trimmed later)
        zs = np.linspace(z0, C_[1], 26)[:-1]
        r_sh = np.interp(zs, [B_[1], C_[1]], [B_[0], C_[0]], left=None, right=None)
        r_sh = B_[0] + (zs - B_[1]) * (C_[0] - B_[0]) / (C_[1] - B_[1])
        r_ma = np.interp(zs, RZm[:, 1], RZm[:, 0])             # constant extrapolation upstream of the nose
        r_mid = 0.5 * (r_sh + r_ma)
        phi = np.arctan2(M[0, :, 1], M[0, :, 0])
        U = np.stack([np.c_[r_mid * np.cos(p), r_mid * np.sin(p), zs] for p in phi], axis=1)
        U[:, 0, 1] = 0.0; U[:, -1, 0] = 0.0
        print(f"   lip: shroud B {B_.round(4)} C {C_.round(4)}, mantle-mid r {r_mid.min():.4f}..{r_mid.max():.4f}, "
              f"ring-1 mid start r {np.hypot(*M[0, :, :2].T).mean():.4f} z {M[0, :, 2].mean():.4f}")
        M = np.concatenate([U, M], axis=0)
        info["lip_mid_r_at_dt_interface"] = float(0.5 * (C_[0] + np.interp(C_[1], RZm[:, 1], RZm[:, 0])))
        info["lip_z_gv"] = float(B_[1]); info["lip_z_dt"] = float(C_[1])
    M = resample_columns(M)
    face, dev = fit_surface(mirror_ext(M), name)
    breptools.Write(face, f"{OUT}/{name}.brep")
    np.save(f"{OUT}/{name}_grid.npy", M)
    mid_faces[name] = dict(skins=[a, b], median_gap=gap, max_gap=gmax, fit_dev=dev, lip=bool(is_lip))


# ---------------------------------------------------------------- 3. body, hub chain, interfaces
body, caps = db.outer_body()
breptools.Write(body, f"{OUT}/body.brep")
for k in ("bulb", "hub0", "hub1", "hubEnd"):
    breptools.Write(faces[F[k]], f"{OUT}/{k}.brep")


def revolve_meridian(RZ, name):
    pts = TColgp_Array1OfPnt(1, len(RZ))
    for k, (r, z) in enumerate(RZ):
        pts.SetValue(k + 1, gp_Pnt(r, 0, z))
    crv = GeomAPI_PointsToBSpline(pts, 1, 3, GeomAbs_C2, 1e-7).Curve()
    e = BRepBuilderAPI_MakeEdge(crv).Edge()
    f = BRepPrimAPI_MakeRevol(e, gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), np.pi / 2).Shape()
    breptools.Write(f, f"{OUT}/{name}.brep")


R_EXT = 2.3
for k in ("if_gv", "if_dt"):
    m = db.meridian(faces[F[k]])
    P = np.array([[np.hypot(p[0], p[1]), p[2]] for e in m for p in db.ends(e)])
    P = P[np.argsort(P[:, 0])]
    lo, hi = P[0], P[-1]
    d = (hi - lo) / np.linalg.norm(hi - lo)
    dev = max(abs(np.cross(d, q - lo)) for q in P)
    assert dev < 1e-6, f"{k}: interface meridian not straight ({dev})"   # straight: extend along the line
    revolve_meridian(np.array([lo, lo + d * (R_EXT - lo[0]) / d[0]]), f"{k}_ext")
    print(f"   {k}: meridian {lo.round(4)} -> {hi.round(4)}, extended to r {R_EXT}")

info.update(pairs=[list(p) for p in pairs], mid=mid_faces, r_ext=R_EXT)
json.dump(info, open(f"{OUT}/parts.json", "w"), indent=1)
print("[parts] written to", OUT)
