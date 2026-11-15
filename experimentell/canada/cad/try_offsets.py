import sys
from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_MakeOffsetShape
from OCC.Core.BRepOffset import BRepOffset_Skin
from OCC.Core.GeomAbs import GeomAbs_Arc, GeomAbs_Intersection
from OCC.Core.BRepCheck import BRepCheck_Analyzer
from OCC.Core.GProp import GProp_GProps
from OCC.Core.BRepGProp import brepgprop
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE
sys.argv=['x','0.0354']; exec(open('/work/canada/cad_envelope.py').read().split("for name, idx in RINGS.items():")[0])
from OCC.Core.TopoDS import topods
def build(idx):
    sew=BRepBuilderAPI_Sewing(1e-5)
    for i in idx: sew.Add(faces[i])
    sew.Perform(); closed,_=close_with_caps(sew.SewedShape())
    ms=BRepBuilderAPI_MakeSolid(); ex=TopExp_Explorer(closed,TopAbs_SHELL)
    while ex.More(): ms.Add(topods.Shell(ex.Current())); ex.Next()
    s=ms.Solid(); f=ShapeFix_Solid(s); f.Perform(); return f.Solid()
solid=build(RINGS["ring2"])
for label,kw in [("arc tol1e-3",dict(tol=1e-3,join=GeomAbs_Arc,rm=False)),("arc tol1e-3 removeIntEdges",dict(tol=1e-3,join=GeomAbs_Arc,rm=True)),
                 ("intersection",dict(tol=1e-4,join=GeomAbs_Intersection,rm=False)),("arc tol1e-2",dict(tol=1e-2,join=GeomAbs_Arc,rm=True))]:
    try:
        off=BRepOffsetAPI_MakeOffsetShape(); off.PerformByJoin(solid,0.0354,kw["tol"],BRepOffset_Skin,False,False,kw["join"],kw["rm"])
        ok=off.IsDone(); v=volume(off.Shape()) if ok else 0; val=BRepCheck_Analyzer(off.Shape()).IsValid() if ok else False
        print(f"[try] {label}: done {ok} valid {val} volume {v:.4f}",flush=True)
    except Exception as e: print(f"[try] {label}: exception {e}",flush=True)
try:
    off=BRepOffsetAPI_MakeOffsetShape(); off.PerformBySimple(solid,0.0354)
    ok=off.IsDone(); print(f"[try] simple: done {ok} valid {BRepCheck_Analyzer(off.Shape()).IsValid() if ok else False} volume {volume(off.Shape()) if ok else 0:.4f}",flush=True)
except Exception as e: print("[try] simple: exception",e,flush=True)
