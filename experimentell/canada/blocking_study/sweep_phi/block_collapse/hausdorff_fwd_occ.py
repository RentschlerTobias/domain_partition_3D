"""mesh -> exact CAD distance (pythonOCC) for boundary nodes / quad centres of a replay output (generic).
usage: LD_LIBRARY_PATH=/dtOO-install/lib /opt/venv-dtoo/bin/python hausdorff_fwd_occ.py shape.brep pts.npz out.json [n]"""
import sys, json, numpy as np
from OCC.Core.BRepTools import breptools
from OCC.Core.BRep import BRep_Builder
from OCC.Core.TopoDS import TopoDS_Shape, topods
from OCC.Core.TopExp import topexp
from OCC.Core.TopTools import TopTools_IndexedMapOfShape
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.BRepExtrema import BRepExtrema_DistShapeShape
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
from OCC.Core.gp import gp_Pnt
from OCC.Core.BRepBndLib import brepbndlib
from OCC.Core.Bnd import Bnd_Box
sh = TopoDS_Shape(); breptools.Read(sh, sys.argv[1], BRep_Builder())
m = TopTools_IndexedMapOfShape(); topexp.MapShapes(sh, TopAbs_FACE, m)
faces = [topods.Face(m.FindKey(i)) for i in range(1, m.Size() + 1)]
boxes = []
for f in faces:
    b = Bnd_Box(); brepbndlib.Add(f, b); boxes.append(np.array(b.Get()))
def dist(p):
    best = 1e9
    bd = [np.linalg.norm(np.maximum(0, np.maximum(bx[:3] - p, p - bx[3:]))) for bx in boxes]
    for i in np.argsort(bd):
        if bd[i] >= best: break
        d = BRepExtrema_DistShapeShape(BRepBuilderAPI_MakeVertex(gp_Pnt(*p)).Vertex(), faces[i]); d.Perform()
        best = min(best, d.Value())
    return best
D = np.load(sys.argv[2]); n = int(sys.argv[4]) if len(sys.argv) > 4 else 6000
rng = np.random.default_rng(0); out = {}
for name in ("nodes", "centres"):
    P = D[name]; sel = rng.choice(len(P), min(n, len(P)), replace=False)
    d = np.array([dist(P[i]) for i in sel]); j = d.argmax()
    out[name] = dict(n=int(len(sel)), max=float(d.max()), p99=float(np.percentile(d, 99)), mean=float(d.mean()), argmax=P[sel][j].tolist())
    print(name, out[name], flush=True)
json.dump(out, open(sys.argv[3], "w"), indent=1)
