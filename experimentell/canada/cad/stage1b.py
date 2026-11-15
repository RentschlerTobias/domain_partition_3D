import sys
sys.argv=['x','0.0354']; exec(open('/work/canada/cad_envelope.py').read().split("for name, idx in RINGS.items():")[0])
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_MakeOffsetShape
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.BOPAlgo import BOPAlgo_ArgumentAnalyzer
from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh
from OCC.Core.StlAPI import StlAPI_Writer
from OCC.Core.TopExp import topexp
from OCC.Core.TopTools import TopTools_IndexedDataMapOfShapeListOfShape
from OCC.Core.BRepLProp import BRepLProp_SLProps
from OCC.Core.BRep import BRep_Tool
from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnSurf
import numpy as np
def build(idx):
    sew=BRepBuilderAPI_Sewing(1e-5)
    for i in idx: sew.Add(faces[i])
    sew.Perform(); closed,_=close_with_caps(sew.SewedShape())
    ms=BRepBuilderAPI_MakeSolid(); ex=TopExp_Explorer(closed,TopAbs_SHELL)
    while ex.More(): ms.Add(topods.Shell(ex.Current())); ex.Next()
    s=ms.Solid(); f=ShapeFix_Solid(s); f.Perform(); return f.Solid()
for name,idx in RINGS.items():
    solid=build(idx)
    off=BRepOffsetAPI_MakeOffsetShape(); off.PerformBySimple(solid,delta); env=off.Shape()
    x0,y0,z0,x1,y1,z1=bbox(env)
    box=BRepPrimAPI_MakeBox(gp_Pnt(0,0,z0-1),gp_Pnt(x1+1,y1+1,z1+1)).Shape()
    envq=BRepAlgoAPI_Common(env,box).Shape()
    an=BOPAlgo_ArgumentAnalyzer(); an.SetShape1(envq); an.SetSelfInterMode(True) if hasattr(an,'SetSelfInterMode') else None; an.Perform()
    print(f"[cad] {name}: simple offset valid {BRepCheck_Analyzer(envq).IsValid()}, self-intersections {an.HasFaulty()}, faces {count(envq,TopAbs_FACE)}, edges {count(envq,TopAbs_EDGE)}, volume {volume(envq):.4f} (ring {volume(solid):.4f})",flush=True)
    # dihedral along every edge between two faces (sampled at the edge middle)
    m=TopTools_IndexedDataMapOfShapeListOfShape(); topexp.MapShapesAndAncestors(envq,TopAbs_EDGE,TopAbs_FACE,m)
    for i in range(1,m.Size()+1):
        e=topods.Edge(m.FindKey(i)); fl=m.FindFromIndex(i)
        if fl.Size()!=2: continue
        it=iter([topods.Face(f) for f in fl])
        fa,fb=next(it),next(it)
        c,f0,f1=BRep_Tool.Curve(e)
        if c is None: continue
        p=c.Value(0.5*(f0+f1)); ns=[]
        for f in (fa,fb):
            srf=BRep_Tool.Surface(f); pr=GeomAPI_ProjectPointOnSurf(p,srf); u,v=pr.LowerDistanceParameters()
            pr_=BRepLProp_SLProps(BRepAdaptor_Surface(f),u,v,1,1e-6); n=pr_.Normal(); n=np.array([n.X(),n.Y(),n.Z()])
            if f.Orientation()==1: n=-n
            ns.append(n)
        ang=np.degrees(np.arccos(np.clip(ns[0]@ns[1],-1,1)))
        if ang>5: print(f"   edge at ({p.X():.2f},{p.Y():.2f},{p.Z():.2f}): dihedral {ang:.0f}")
    BRepMesh_IncrementalMesh(envq,0.005,False,0.2,True)
    w=StlAPI_Writer(); w.Write(envq,f"/work/canada/cad/{name}_envelope.stl")
    from OCC.Core.STEPControl import STEPControl_Writer, STEPControl_AsIs
    w=STEPControl_Writer(); w.Transfer(envq,STEPControl_AsIs); w.Write(f"/work/canada/cad/{name}_envelope.stp")
    BRepMesh_IncrementalMesh(solid,0.005,False,0.2,True); w=StlAPI_Writer(); w.Write(solid,f"/work/canada/cad/{name}_solid.stl")
