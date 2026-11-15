import numpy as np
from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.TopoDS import topods
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib
r = STEPControl_Reader(); r.ReadFile("/work/canada/case/cad_blk.stp"); r.TransferRoots(); s = r.OneShape()
ex = TopExp_Explorer(s, TopAbs_FACE); i=0
while ex.More():
    f = topods.Face(ex.Current()); b = Bnd_Box(); brepbndlib.Add(f, b)
    x0,y0,z0,x1,y1,z1 = b.Get()
    srf = BRep_Tool.Surface(f); u0,u1,v0,v1 = srf.Bounds()
    pts = np.array([[*srf.Value(u,v).Coord()] for u in np.linspace(u0,u1,9) for v in np.linspace(v0,v1,9)])
    rr = np.hypot(pts[:,0],pts[:,1])
    print(i, BRepAdaptor_Surface(f).GetType(), f"bb x[{x0:.2f},{x1:.2f}] y[{y0:.2f},{y1:.2f}] z[{z0:.2f},{z1:.2f}] r[{rr.min():.3f},{rr.max():.3f}]")
    i+=1; ex.Next()
