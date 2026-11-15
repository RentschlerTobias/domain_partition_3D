import numpy as np
from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.TopExp import TopExp_Explorer, topexp
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_EDGE
from OCC.Core.TopoDS import topods
from OCC.Core.TopTools import TopTools_IndexedDataMapOfShapeListOfShape
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib
from OCC.Core.GProp import GProp_GProps
from OCC.Core.BRepGProp import brepgprop
r = STEPControl_Reader(); r.ReadFile("/work/canada/case/cad_blk.stp"); r.TransferRoots(); shape = r.OneShape()
faces = []; ex = TopExp_Explorer(shape, TopAbs_FACE)
while ex.More(): faces.append(topods.Face(ex.Current())); ex.Next()
m = TopTools_IndexedDataMapOfShapeListOfShape(); topexp.MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, m)
def fid(f):
    for i, g in enumerate(faces):
        if g.IsSame(f): return i
adj = {i: set() for i in range(len(faces))}
for k in range(1, m.Size() + 1):
    L = [fid(x) for x in m.FindFromIndex(k)]
    for a in L:
        for b in L:
            if a != b: adj[a].add(b)
for i, f in enumerate(faces):
    b = Bnd_Box(); brepbndlib.Add(f, b); bb = np.array(b.Get())
    g = GProp_GProps(); brepgprop.SurfaceProperties(f, g)
    print(i, BRepAdaptor_Surface(f).GetType(), "area %.3f" % g.Mass(), "bbox", bb.round(2).tolist(), "adj", sorted(adj[i]))
