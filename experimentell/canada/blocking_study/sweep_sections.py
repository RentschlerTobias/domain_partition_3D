"""Is the draft-tube quarter domain a rotational sweep?  Cut the CAD solid with meridional half planes
phi = const and report the topology of each section (faces, holes per face, Euler characteristic)."""
import sys, math
from OCC.Core.BRepTools import breptools
from OCC.Core.BRep import BRep_Builder
from OCC.Core.TopoDS import TopoDS_Shape, topods
from OCC.Core.gp import gp_Pnt, gp_Dir, gp_Ax3, gp_Pln
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_WIRE, TopAbs_EDGE
from OCC.Core.GProp import GProp_GProps
from OCC.Core.BRepGProp import brepgprop
s = TopoDS_Shape(); breptools.Read(s, sys.argv[1], BRep_Builder())
def items(sh, t):
    e = TopExp_Explorer(sh, t); out = []
    while e.More(): out.append(e.Current()); e.Next()
    return out
for deg in [0.5, 10, 20, 30, 45, 60, 70, 80, 89.5]:
    p = math.radians(deg)
    # half plane: spanned by radial dir (cos p, sin p, 0) and z, radial coordinate u in [0, 20]
    ax = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(-math.sin(p), math.cos(p), 0), gp_Dir(math.cos(p), math.sin(p), 0))
    f = BRepBuilderAPI_MakeFace(gp_Pln(ax), 0.0, 20.0, -10.0, 20.0).Face()
    c = BRepAlgoAPI_Common(s, f); c.Build()
    r = c.Shape(); fs = items(r, TopAbs_FACE)
    desc = []
    for fc in fs:
        g = GProp_GProps(); brepgprop.SurfaceProperties(fc, g)
        nw = len(items(fc, TopAbs_WIRE)); ne = len(items(fc, TopAbs_EDGE))
        desc.append(f"area {g.Mass():.3f} wires {nw} (holes {nw-1}) edges {ne}")
    print(f"phi {deg:5.1f}: {len(fs)} section faces; " + " | ".join(desc), flush=True)
