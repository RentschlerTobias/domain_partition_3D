"""Stage 1 of the CAD envelope: ring bodies of the canada diffuser as solids, offset by delta.

Reads the dtOO STEP export (cad_blk.stp, quarter model x>=0, y>=0), sews the faces of each ring
(inner skin, outer skin, front bevel, back bevel) into a shell, closes it with planar caps on the two
symmetry planes, makes a solid, offsets it by delta with rounded joins (BRepOffsetAPI_MakeOffsetShape,
GeomAbs_Arc) and trims the result back to the quarter. Writes one STEP per envelope plus a summary of
its faces and edges.

run with LD_LIBRARY_PATH=/dtOO-install/lib /opt/venv-dtoo/bin/python cad_envelope.py <delta>
"""
import sys
from OCC.Core.STEPControl import STEPControl_Reader, STEPControl_Writer, STEPControl_AsIs
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_EDGE, TopAbs_SHELL, TopAbs_SOLID, TopAbs_WIRE
from OCC.Core.TopoDS import topods
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Sewing, BRepBuilderAPI_MakeSolid
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_MakeOffsetShape
from OCC.Core.BRepOffset import BRepOffset_Skin
from OCC.Core.GeomAbs import GeomAbs_Arc
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Core.BRepCheck import BRepCheck_Analyzer
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib
from OCC.Core.GProp import GProp_GProps
from OCC.Core.BRepGProp import brepgprop
from OCC.Core.gp import gp_Pnt
from OCC.Core.ShapeFix import ShapeFix_Solid, ShapeFix_Shell
from OCC.Core.BRepFill import BRepFill_Filling
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeWire
from OCC.Core.ShapeAnalysis import ShapeAnalysis_FreeBounds
from OCC.Core.TopTools import TopTools_HSequenceOfShape

delta = float(sys.argv[1]) if len(sys.argv) > 1 else 0.035

r = STEPControl_Reader(); r.ReadFile("/work/canada/case/cad_blk.stp"); r.TransferRoots(); shape = r.OneShape()
faces = []
ex = TopExp_Explorer(shape, TopAbs_FACE)
while ex.More():
    faces.append(topods.Face(ex.Current())); ex.Next()
# face indices from the export order (see export_step.py / inspection): ring 2 = dt2In, dt2Out,
# bevelFront_dt2, bevelBack_dt2; ring 3 = dt3In, dt3Out, bevelFront_dt3, bevelBack_dt3
RINGS = {"ring2": [1, 4, 6, 7], "ring3": [2, 5, 9, 10]}


def bbox(s):
    b = Bnd_Box(); brepbndlib.Add(s, b); return b.Get()


def count(s, t):
    ex = TopExp_Explorer(s, t); n = 0
    while ex.More():
        n += 1; ex.Next()
    return n


def volume(s):
    g = GProp_GProps(); brepgprop.VolumeProperties(s, g); return g.Mass()


def close_with_caps(shell):
    """free boundary wires of the sewn shell (on the symmetry planes) -> planar cap faces"""
    fb = ShapeAnalysis_FreeBounds(shell, 1e-6, False, False)
    wires = fb.GetClosedWires()
    sew = BRepBuilderAPI_Sewing(1e-5)
    sew.Add(shell)
    ex = TopExp_Explorer(wires, TopAbs_WIRE)
    n = 0
    while ex.More():
        w = topods.Wire(ex.Current())
        mf = BRepBuilderAPI_MakeFace(w, True)       # planar if possible
        if mf.IsDone():
            sew.Add(mf.Face()); n += 1
        ex.Next()
    sew.Perform()
    return sew.SewedShape(), n


for name, idx in RINGS.items():
    sew = BRepBuilderAPI_Sewing(1e-5)
    for i in idx:
        sew.Add(faces[i])
    sew.Perform()
    sh = sew.SewedShape()
    closed, ncap = close_with_caps(sh)
    shells = []
    ex = TopExp_Explorer(closed, TopAbs_SHELL)
    while ex.More():
        shells.append(topods.Shell(ex.Current())); ex.Next()
    ms = BRepBuilderAPI_MakeSolid()
    for s_ in shells:
        ms.Add(s_)
    solid = ms.Solid()
    fix = ShapeFix_Solid(solid); fix.Perform(); solid = fix.Solid()
    ok = BRepCheck_Analyzer(solid).IsValid()
    print(f"[cad] {name}: {len(idx)} faces + {ncap} caps -> {len(shells)} shell(s), solid valid {ok}, "
          f"volume {volume(solid):.5f}, bbox {[round(v, 2) for v in bbox(solid)]}", flush=True)
    off = BRepOffsetAPI_MakeOffsetShape()
    off.PerformByJoin(solid, delta, 1e-5, BRepOffset_Skin, False, False, GeomAbs_Arc, False)
    if not off.IsDone():
        print(f"[cad] {name}: offset FAILED", flush=True)
        continue
    env = off.Shape()
    x0, y0, z0, x1, y1, z1 = bbox(env)
    box = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, z0 - 1), gp_Pnt(x1 + 1, y1 + 1, z1 + 1)).Shape()
    env_q = BRepAlgoAPI_Common(env, box).Shape()
    ok2 = BRepCheck_Analyzer(env_q).IsValid()
    print(f"[cad] {name}: envelope offset {delta}: valid {ok2}, faces {count(env_q, TopAbs_FACE)}, "
          f"edges {count(env_q, TopAbs_EDGE)}, volume {volume(env_q):.5f} (ring {volume(solid):.5f})", flush=True)
    w = STEPControl_Writer(); w.Transfer(env_q, STEPControl_AsIs); w.Write(f"/work/canada/cad/{name}_envelope.stp")
    w = STEPControl_Writer(); w.Transfer(solid, STEPControl_AsIs); w.Write(f"/work/canada/cad/{name}_solid.stp")
