import sys
from OCC.Core.BRepTools import breptools
from OCC.Core.BRep import BRep_Builder
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Iterator
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_SOLID, TopAbs_EDGE
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
s = TopoDS_Shape(); breptools.Read(s, sys.argv[1], BRep_Builder())
print("type", s.ShapeType())
it = TopoDS_Iterator(s)
while it.More():
    c = it.Value(); n=0; e=TopExp_Explorer(c, TopAbs_FACE)
    while e.More(): n+=1; e.Next()
    print(" child type", c.ShapeType(), "faces", n); it.Next()
