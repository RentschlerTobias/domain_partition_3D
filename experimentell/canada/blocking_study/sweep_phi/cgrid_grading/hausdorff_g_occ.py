"""Hausdorff of the C-grid mesh boundary to the REAL CAD: case/cad_blk.stp faces (walls with thickness, far field,
hub chain, interfaces) + the planar faces of the domain that are symmetry planes or inlet/outlet (domain bbox
extremes). Forward: mesh boundary nodes / quad centres -> CAD (with the nearest face, grouped wall / other).
Reverse samples: wall faces + far field + hub faces + planes (interface discs extend into the hub / lip and are
used forward only)."""
import sys, json, numpy as np
sys.path.insert(0, "/work/canada/blocking_study/sweep_phi")
from section_rz import load
from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_IN, TopAbs_ON
from OCC.Core.TopoDS import topods
from OCC.Core.BRepExtrema import BRepExtrema_DistShapeShape
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.GeomAbs import GeomAbs_Plane
from OCC.Core.gp import gp_Pnt, gp_Pnt2d
from OCC.Core.BRepTools import breptools
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepClass import BRepClass_FaceClassifier
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib
r = STEPControl_Reader(); r.ReadFile("/work/canada/case/cad_blk.stp"); r.TransferRoots(); sh = r.OneShape()
cad = []; ex = TopExp_Explorer(sh, TopAbs_FACE)
while ex.More(): cad.append(topods.Face(ex.Current())); ex.Next()
W = json.load(open("/work/canada/blocking_study/sweep_phi/cgrid/walls.json")); wall_faces = sorted({f for c in W["components"] for f in c["faces"]})
solid, sfaces, _ = load()
bb = Bnd_Box(); brepbndlib.Add(solid, bb); B = np.array(bb.Get())
planes = []
for f in sfaces:
    A = BRepAdaptor_Surface(topods.Face(f))
    if A.GetType() != GeomAbs_Plane: continue
    n = np.array(A.Plane().Axis().Direction().Coord()); p = np.array(A.Plane().Location().Coord())
    if abs(n[2]) < 1e-9 and abs(p @ n) < 1e-9:
        planes.append(topods.Face(f))                                   # plane through the axis (symmetry)
    elif abs(abs(n[2]) - 1) < 1e-9 and (abs(p[2] - B[2]) < 1e-4 or abs(p[2] - B[5]) < 1e-4):
        planes.append(topods.Face(f))                                   # inlet / outlet
print('planes found', len(planes)); faces = cad + planes
group = ["wall" if i in wall_faces else "cad" for i in range(len(cad))] + ["plane"] * len(planes)
boxes = []
for f in faces:
    b = Bnd_Box(); brepbndlib.Add(f, b); boxes.append(np.array(b.Get()))
def dist(p):
    best = (1e9, -1)
    bd = [np.linalg.norm(np.maximum(0, np.maximum(bx[:3] - p, p - bx[3:]))) for bx in boxes]
    for i in np.argsort(bd):
        if bd[i] >= best[0]: break
        d = BRepExtrema_DistShapeShape(BRepBuilderAPI_MakeVertex(gp_Pnt(*p)).Vertex(), faces[i]); d.Perform()
        if d.Value() < best[0]: best = (d.Value(), i)
    return best
D = np.load("boundary_pts.npz"); rng = np.random.default_rng(0); out = {}
NS = int(sys.argv[1]) if len(sys.argv) > 1 else 6000
for name in ("nodes", "centres"):
    P = D[name]; sel = rng.choice(len(P), min(len(P), NS), replace=False)
    res = np.array([dist(P[i]) for i in sel]); d, fi = res[:, 0], res[:, 1].astype(int)
    g = np.array([group[i] for i in fi])
    out[name] = d; out[name + "_pts"] = P[sel]; out[name + "_grp"] = g
    print(f"mesh->CAD {name}: max {d.max():.3e} (wall {d[g == 'wall'].max():.3e}, n_wall {int((g == 'wall').sum())}), p99 {np.percentile(d, 99):.3e}", flush=True)
# reverse samples: wall faces (walls.json) and planes everywhere; every other CAD face only where it lies on the
# boundary of the baffle domain (drops e.g. the parts of the interface discs that are inside the hub chain)
def on_domain(q):
    v = BRepBuilderAPI_MakeVertex(gp_Pnt(*q)).Vertex()
    for f in sfaces:
        d = BRepExtrema_DistShapeShape(v, f); d.Perform()
        if d.Value() < 1e-4:
            return True
    return False
S, F = [], []
for fi in range(len(faces)):
    f = faces[fi]; u0, u1, v0, v1 = breptools.UVBounds(f); s = BRep_Tool.Surface(f)
    for u in np.linspace(u0, u1, 70):
        for v in np.linspace(v0, v1, 70):
            if BRepClass_FaceClassifier(f, gp_Pnt2d(u, v), 1e-7).State() in (TopAbs_IN, TopAbs_ON):
                q = np.array(s.Value(u, v).Coord())
                if q[0] >= -1e-9 and q[1] >= -1e-9 and (group[fi] != "cad" or on_domain(q)):
                    S.append(q); F.append(fi)
np.savez("hausdorff_g_occ.npz", cad=np.array(S), cad_face=np.array(F), cad_grp=np.array([group[i] for i in F]), **out)
print("reverse samples", len(S))
