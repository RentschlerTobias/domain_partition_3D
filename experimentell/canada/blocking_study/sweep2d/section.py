"""Meridional section (half plane phi) of the draft-tube CAD domain -> planar BREP face in local (r, z)
coordinates (rotated into the x-z plane) + label per boundary edge."""
import sys, math, json
sys.path.insert(0, "/work/canada/blocking_study")
from OCC.Core.BRepTools import breptools
from OCC.Core.BRep import BRep_Builder
from OCC.Core.TopoDS import TopoDS_Shape, topods
from OCC.Core.gp import gp_Pnt, gp_Dir, gp_Ax3, gp_Pln, gp_Trsf, gp_Ax1
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_Transform, BRepBuilderAPI_MakeVertex
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_EDGE
from OCC.Core.BRepAdaptor import BRepAdaptor_Curve
from OCC.Core.BRepExtrema import BRepExtrema_DistShapeShape
src, deg, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
s = TopoDS_Shape(); breptools.Read(s, src, BRep_Builder())
def items(sh, t):
    e = TopExp_Explorer(sh, t); o = []
    while e.More(): o.append(e.Current()); e.Next()
    return o
p = math.radians(deg)
ax = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(-math.sin(p), math.cos(p), 0), gp_Dir(math.cos(p), math.sin(p), 0))
f = BRepBuilderAPI_MakeFace(gp_Pln(ax), 0.0, 20.0, -10.0, 20.0).Face()
c = BRepAlgoAPI_Common(s, f); c.Build(); fc = items(c.Shape(), TopAbs_FACE)[0]
tr = gp_Trsf(); tr.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), -p)     # rotate into the x-z plane
fc2 = topods.Face(BRepBuilderAPI_Transform(fc, tr, True).Shape())
breptools.Write(fc2, out)
lab = {d['face']: d['label'] for d in json.load(open('/work/canada/cad/domain_blk_faces.json'))}
dfaces = items(s, TopAbs_FACE)
labels = []
for e in items(fc, TopAbs_EDGE):
    a = BRepAdaptor_Curve(topods.Edge(e)); m = a.Value(0.5 * (a.FirstParameter() + a.LastParameter()))
    v = BRepBuilderAPI_MakeVertex(m).Vertex(); best = (1e9, None)
    for i, F in enumerate(dfaces):
        d = BRepExtrema_DistShapeShape(v, F); d.Perform()
        if d.Value() < best[0]: best = (d.Value(), i)
    labels.append(lab.get(best[1], str(best[1])))
json.dump(labels, open(out.replace(".brep", "_edges.json"), "w"))
print(f"phi {deg}: section face with {len(labels)} edges written to {out}")
