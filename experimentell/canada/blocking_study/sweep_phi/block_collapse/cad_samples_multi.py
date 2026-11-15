"""CAD face samples from several shape files (no shared topology needed), restricted to a box/half-spaces.
For CAD that comes as unconnected faces (e.g. case/cad_blk.stp) the block-collapse code derives the CAD edge /
vertex topology from the block-face labels instead (block_collapse.py --cad-topology labels).
Here (C-grid): all faces of cad_blk.stp + the planar faces of the baffle domain that are symmetry planes (through
the axis) or inlet/outlet (domain bbox z extremes) -- the same selection as cgrid/hausdorff_c_occ.py.
Samples with x < 0 or y < 0 are dropped (quarter model).
usage: LD_LIBRARY_PATH=/dtOO-install/lib /opt/venv-dtoo/bin/python cad_samples_multi.py step brep out.npz [ds]"""
import sys, json, numpy as np
from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_IN, TopAbs_ON, TopAbs_SOLID
from OCC.Core.TopoDS import topods, TopoDS_Shape, TopoDS_Iterator
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.GeomAbs import GeomAbs_Plane
from OCC.Core.gp import gp_Pnt2d
from OCC.Core.BRepTools import breptools
from OCC.Core.BRep import BRep_Tool, BRep_Builder
from OCC.Core.BRepClass import BRepClass_FaceClassifier
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib
step, brep, out = sys.argv[1:4]; DS = float(sys.argv[4]) if len(sys.argv) > 4 else 0.03
r = STEPControl_Reader(); r.ReadFile(step); r.TransferRoots(); sh = r.OneShape()
faces = []; ex = TopExp_Explorer(sh, TopAbs_FACE)
while ex.More(): faces.append(topods.Face(ex.Current())); ex.Next()
ncad = len(faces)
s = TopoDS_Shape(); breptools.Read(s, brep, BRep_Builder()); solid = None
it = TopoDS_Iterator(s)
while it.More():
    if it.Value().ShapeType() == TopAbs_SOLID: solid = it.Value()
    it.Next()
bb = Bnd_Box(); brepbndlib.Add(solid, bb); B = np.array(bb.Get())
ex = TopExp_Explorer(solid, TopAbs_FACE)
while ex.More():
    f = topods.Face(ex.Current()); A = BRepAdaptor_Surface(f)
    if A.GetType() == GeomAbs_Plane:
        n = np.array(A.Plane().Axis().Direction().Coord()); p = np.array(A.Plane().Location().Coord())
        if (abs(n[2]) < 1e-9 and abs(p @ n) < 1e-9) or (abs(abs(n[2]) - 1) < 1e-9 and (abs(p[2] - B[2]) < 1e-4 or abs(p[2] - B[5]) < 1e-4)):
            faces.append(f)
    ex.Next()
print(f"{ncad} CAD faces + {len(faces) - ncad} planes", flush=True)
P, I = [], []
for fi, f in enumerate(faces):
    u0, u1, v0, v1 = breptools.UVBounds(f); S = BRep_Tool.Surface(f)
    G = np.array([[S.Value(x, y).Coord() for y in np.linspace(v0, v1, 9)] for x in np.linspace(u0, u1, 9)])
    Lu = np.linalg.norm(np.diff(G, axis=0), axis=-1).sum(0).max(); Lv = np.linalg.norm(np.diff(G, axis=1), axis=-1).sum(1).max()
    nu = int(min(1500, max(8, Lu / DS))); nv = int(min(1500, max(8, Lv / DS))); c = 0
    for x in np.linspace(u0, u1, nu):
        for y in np.linspace(v0, v1, nv):
            if BRepClass_FaceClassifier(f, gp_Pnt2d(x, y), 1e-7).State() in (TopAbs_IN, TopAbs_ON):
                q = S.Value(x, y).Coord()
                if q[0] >= -1e-9 and q[1] >= -1e-9:
                    P.append(q); I.append(fi); c += 1
    print(f"face {fi}: {c} samples", flush=True)
np.savez(out, face_pts=np.array(P), face_id=np.array(I), n_cad=ncad, src=json.dumps([step, brep]))
