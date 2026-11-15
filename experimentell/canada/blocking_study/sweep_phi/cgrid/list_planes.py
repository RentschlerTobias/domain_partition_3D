import sys; sys.path.insert(0, "/work/canada/blocking_study/sweep_phi")
import numpy as np
from section_rz import load
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.GeomAbs import GeomAbs_Plane
from OCC.Core.TopoDS import topods
s, sf, _ = load()
for i, f in enumerate(sf):
    A = BRepAdaptor_Surface(topods.Face(f))
    if A.GetType() == GeomAbs_Plane:
        print(i, np.round(A.Plane().Axis().Direction().Coord(), 4), np.round(A.Plane().Location().Coord(), 4))
