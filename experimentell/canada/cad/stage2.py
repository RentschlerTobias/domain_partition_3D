"""Blunt the knife edges of each ring solid (pipe of radius rho along every edge with interior angle
< 60 deg subtracted), then offset by delta (arc, else simple), trim to the quarter, check, export."""
import sys
sys.argv=['x','0.0354']; exec(open('/work/canada/cad_envelope.py').read().split("for name, idx in RINGS.items():")[0])
import numpy as np
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_MakeOffsetShape, BRepOffsetAPI_MakePipe
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface, BRepAdaptor_Curve
from OCC.Core.BOPAlgo import BOPAlgo_ArgumentAnalyzer
from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh
from OCC.Core.StlAPI import StlAPI_Writer
from OCC.Core.TopExp import topexp
from OCC.Core.TopTools import TopTools_IndexedDataMapOfShapeListOfShape
from OCC.Core.BRepLProp import BRepLProp_SLProps
from OCC.Core.BRep import BRep_Tool
from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnSurf
from OCC.Core.gp import gp_Circ, gp_Ax2, gp_Dir, gp_Vec
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeFace
from OCC.Core.GeomAbs import GeomAbs_Arc
from OCC.Core.BRepOffset import BRepOffset_Skin
rho = 2 * delta

def build(idx):
    sew=BRepBuilderAPI_Sewing(1e-5)
    for i in idx: sew.Add(faces[i])
    sew.Perform(); closed,_=close_with_caps(sew.SewedShape())
    ms=BRepBuilderAPI_MakeSolid(); ex=TopExp_Explorer(closed,TopAbs_SHELL)
    while ex.More(): ms.Add(topods.Shell(ex.Current())); ex.Next()
    s=ms.Solid(); f=ShapeFix_Solid(s); f.Perform(); return f.Solid()

def edge_dihedral(shape):
    m=TopTools_IndexedDataMapOfShapeListOfShape(); topexp.MapShapesAndAncestors(shape,TopAbs_EDGE,TopAbs_FACE,m)
    out=[]
    for i in range(1,m.Size()+1):
        e=topods.Edge(m.FindKey(i)); fl=[topods.Face(f) for f in m.FindFromIndex(i)]
        if len(fl)!=2: continue
        c,f0,f1=BRep_Tool.Curve(e)
        if c is None: continue
        angs=[]
        for t in np.linspace(0.1,0.9,5):
            p=c.Value(f0+t*(f1-f0)); ns=[]
            for f in fl:
                srf=BRep_Tool.Surface(f); pr=GeomAPI_ProjectPointOnSurf(p,srf); u,v=pr.LowerDistanceParameters()
                n=BRepLProp_SLProps(BRepAdaptor_Surface(f),u,v,1,1e-6).Normal(); n=np.array([n.X(),n.Y(),n.Z()])
                if f.Orientation()==1: n=-n
                ns.append(n)
            angs.append(np.degrees(np.arccos(np.clip(ns[0]@ns[1],-1,1))))
        out.append((e,max(angs)))
    return out

for name,idx in RINGS.items():
    solid=build(idx)
    knives=[(e,a) for e,a in edge_dihedral(solid) if a>120]
    print(f"[cad] {name}: knife edges (dihedral > 120): {[round(a) for _,a in knives]}",flush=True)
    blunt=solid
    for e,a in knives:
        c=BRepAdaptor_Curve(e); p=c.Value(c.FirstParameter()); d=gp_Vec(); pp=gp_Pnt(); c.D1(c.FirstParameter(),pp,d)
        circ=gp_Circ(gp_Ax2(p,gp_Dir(d)),rho)
        prof=BRepBuilderAPI_MakeFace(BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(circ).Edge()).Wire()).Face()
        pipe=BRepOffsetAPI_MakePipe(BRepBuilderAPI_MakeWire(e).Wire(),prof); pipe.Build()
        cutter=pipe.Shape()
        blunt=BRepAlgoAPI_Cut(blunt,cutter).Shape()
    print(f"[cad] {name}: blunted valid {BRepCheck_Analyzer(blunt).IsValid()}, volume {volume(blunt):.4f} (ring {volume(solid):.4f}), "
          f"remaining edges > 120: {[round(a) for _,a in edge_dihedral(blunt) if a>120]}",flush=True)
    env=None
    for mode in ("arc","simple"):
        off=BRepOffsetAPI_MakeOffsetShape()
        try:
            if mode=="arc": off.PerformByJoin(blunt,delta,1e-4,BRepOffset_Skin,False,False,GeomAbs_Arc,False)
            else: off.PerformBySimple(blunt,delta)
        except Exception as ex_: print("   ",mode,"exception",ex_); continue
        if off.IsDone():
            env=off.Shape(); print(f"[cad] {name}: offset by {mode} done",flush=True); break
        print(f"[cad] {name}: offset by {mode} failed",flush=True)
    if env is None: continue
    x0,y0,z0,x1,y1,z1=bbox(env)
    envq=BRepAlgoAPI_Common(env,BRepPrimAPI_MakeBox(gp_Pnt(0,0,z0-1),gp_Pnt(x1+1,y1+1,z1+1)).Shape()).Shape()
    # make sure the result is ONE solid (Common returns a compound); gmsh read only faces before
    sols=[]; ex=TopExp_Explorer(envq,TopAbs_SOLID)
    while ex.More(): sols.append(topods.Solid(ex.Current())); ex.Next()
    print(f"[cad] {name}: solids in the trimmed envelope: {len(sols)}",flush=True)
    if len(sols)==1: envq=sols[0]
    an=BOPAlgo_ArgumentAnalyzer(); an.SetShape1(envq); an.SetSelfInterMode(True); an.Perform()
    print(f"[cad] {name}: envelope valid {BRepCheck_Analyzer(envq).IsValid()}, self-intersections {an.HasFaulty()}, faces {count(envq,TopAbs_FACE)}, "
          f"volume {volume(envq):.4f}; edge dihedrals {sorted(round(a) for _,a in edge_dihedral(envq) if a>5)}",flush=True)
    from OCC.Core.STEPControl import STEPControl_Writer, STEPControl_AsIs
    w=STEPControl_Writer(); w.Transfer(envq,STEPControl_AsIs); w.Write(f"/work/canada/cad/{name}_envelope_b.stp")
    from OCC.Core.BRepTools import breptools
    breptools.Write(envq, f"/work/canada/cad/{name}_envelope_b.brep")
    breptools.Write(blunt, f"/work/canada/cad/{name}_blunt.brep")
    BRepMesh_IncrementalMesh(envq,0.003,False,0.2,True); StlAPI_Writer().Write(envq,f"/work/canada/cad/{name}_envelope_b.stl")
    BRepMesh_IncrementalMesh(blunt,0.003,False,0.2,True); StlAPI_Writer().Write(blunt,f"/work/canada/cad/{name}_blunt.stl")
