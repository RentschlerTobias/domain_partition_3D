"""Hausdorff part 1 (pythonOCC): (a) distance of mesh boundary nodes / quad centres to the CAD (domain.brep:
solid faces + embedded baffle faces), (b) sample points on the CAD faces for the reverse direction."""
import sys, numpy as np
sys.path.insert(0, "/work/canada/blocking_study/sweep_phi")
from section_rz import load, items
from OCC.Core.BRepExtrema import BRepExtrema_DistShapeShape
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
from OCC.Core.gp import gp_Pnt, gp_Pnt2d
from OCC.Core.BRepTools import breptools
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepClass import BRepClass_FaceClassifier
from OCC.Core.TopAbs import TopAbs_IN, TopAbs_ON
from OCC.Core.TopoDS import topods
from OCC.Core.BRepBndLib import brepbndlib
from OCC.Core.Bnd import Bnd_Box
solid, sfaces, free = load()
faces = [topods.Face(f) for f in sfaces + free]
D = np.load("boundary_pts.npz")
boxes = []
for f in faces:
    b = Bnd_Box(); brepbndlib.Add(f, b); boxes.append(np.array(b.Get()))
def dist(p):
    best = 1e9
    # faces sorted by bbox distance, early exit
    bd = [np.linalg.norm(np.maximum(0, np.maximum(bx[:3] - p, p - bx[3:]))) for bx in boxes]
    for i in np.argsort(bd):
        if bd[i] >= best: break
        d = BRepExtrema_DistShapeShape(BRepBuilderAPI_MakeVertex(gp_Pnt(*p)).Vertex(), faces[i]); d.Perform()
        best = min(best, d.Value())
    return best
rng = np.random.default_rng(0)
out = {}
for name in ("nodes", "centres"):
    P = D[name]; sel = rng.choice(len(P), min(len(P), int(sys.argv[1]) if len(sys.argv) > 1 else 8000), replace=False)
    d = np.array([dist(P[i]) for i in sel])
    out[name] = d; out[name + "_pts"] = P[sel]
    print(f"mesh->CAD {name}: {len(sel)} samples, max {d.max():.3e}, p99 {np.percentile(d,99):.3e}, mean {d.mean():.3e}", flush=True)
# CAD samples (reverse direction)
S = []; F = []
for fi, f in enumerate(faces):
    u0, u1, v0, v1 = breptools.UVBounds(f); s = BRep_Tool.Surface(f)
    for u in np.linspace(u0, u1, 60):
        for v in np.linspace(v0, v1, 60):
            if BRepClass_FaceClassifier(f, gp_Pnt2d(u, v), 1e-7).State() in (TopAbs_IN, TopAbs_ON):
                S.append(s.Value(u, v).Coord()); F.append(fi)
np.savez("hausdorff_occ.npz", cad=np.array(S), cad_face=np.array(F), **out)
print("CAD samples", len(S))
