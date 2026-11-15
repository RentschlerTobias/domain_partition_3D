"""Meridional half-plane sections of the baffle domain (dt_baffle/in_a/domain.brep: one solid + embedded
zero-thickness baffle faces). For each phi: the section face of the solid (outer wire, (r,z)) and the
section curves of the free (embedded) faces = slits. Each section edge is tagged with the index of the 3D
CAD face it lies on (geometric id, used only for correspondence between phi, never as a decision).
Writes sections/sec_<phi>.json : {"loops": [[{"face": i, "pts": [[r,z],...]}, ...]], "slits": [...]}"""
import sys, math, json, os
import numpy as np
from OCC.Core.BRepTools import breptools
from OCC.Core.BRep import BRep_Builder, BRep_Tool
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Iterator, topods
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_EDGE, TopAbs_WIRE, TopAbs_SOLID, TopAbs_REVERSED
from OCC.Core.gp import gp_Pnt, gp_Dir, gp_Ax3, gp_Pln
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeVertex
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Section
from OCC.Core.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface
from OCC.Core.GeomAbs import GeomAbs_Plane
from OCC.Core.GCPnts import GCPnts_UniformAbscissa
from OCC.Core.BRepExtrema import BRepExtrema_DistShapeShape
from OCC.Core.BRepTools import BRepTools_WireExplorer

SRC = os.environ.get("HEXDP_DOMAIN_BREP", "domain.brep")   # the domain solid; --brep overrides
def items(sh, t):
    e = TopExp_Explorer(sh, t); o = []
    while e.More(): o.append(e.Current()); e.Next()
    return o
def load():
    s = TopoDS_Shape(); breptools.Read(s, SRC, BRep_Builder())
    solid, free = None, []
    it = TopoDS_Iterator(s)
    while it.More():
        c = it.Value()
        if c.ShapeType() == TopAbs_SOLID: solid = c
        else: free += items(c, TopAbs_FACE)
        it.Next()
    return solid, items(solid, TopAbs_FACE), free

def sample(e, ds=0.01):
    a = BRepAdaptor_Curve(topods.Edge(e)); f0, f1 = a.FirstParameter(), a.LastParameter()
    t = np.linspace(f0, f1, 4001)
    P = np.array([a.Value(x).Coord() for x in t])
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
    n = max(2, int(math.ceil(s[-1] / ds)) + 1)
    u = np.linspace(0, s[-1], n)
    Q = np.c_[[np.interp(u, s, P[:, k]) for k in range(3)]].T
    if e.Orientation() == TopAbs_REVERSED: Q = Q[::-1]
    return Q

def nearest_face(p, fl, nrm=None):
    """index of the CAD face nearest to p; planar faces parallel to the section plane are skipped"""
    v = BRepBuilderAPI_MakeVertex(gp_Pnt(*p)).Vertex(); best = (1e9, -1)
    for i, F in enumerate(fl):
        A = BRepAdaptor_Surface(topods.Face(F))
        if nrm is not None and A.GetType() == GeomAbs_Plane and abs(np.dot(A.Plane().Axis().Direction().Coord(), nrm)) > 0.999:
            continue
        d = BRepExtrema_DistShapeShape(v, F); d.Perform()
        if d.Value() < best[0]: best = (d.Value(), i)
    return best[1]

def section(deg, solid, sfaces, free):
    p = math.radians(deg)
    er = np.array([math.cos(p), math.sin(p), 0.0])
    ax = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(-math.sin(p), math.cos(p), 0), gp_Dir(*er))
    pl = BRepBuilderAPI_MakeFace(gp_Pln(ax), 0.0, 30.0, -20.0, 30.0).Face()
    c = BRepAlgoAPI_Common(solid, pl); c.Build(); fcs = items(c.Shape(), TopAbs_FACE)
    assert len(fcs) == 1, len(fcs)
    fc = topods.Face(fcs[0])
    nrm = np.array([-math.sin(p), math.cos(p), 0.0])
    torz = lambda Q: np.c_[Q @ er, Q[:, 2]]
    loops = []
    for w in items(fc, TopAbs_WIRE):
        L = []
        we = BRepTools_WireExplorer(topods.Wire(w), fc)
        while we.More():
            e = we.Current(); Q = sample(e)
            # orientation from the wire explorer
            if we.Orientation() == TopAbs_REVERSED and e.Orientation() != TopAbs_REVERSED: pass
            rz = torz(Q)
            fi = -1 if np.abs(rz[:, 0]).max() < 1e-6 else nearest_face(Q[len(Q) // 2], sfaces, nrm)  # -1 = axis
            L.append({"face": int(fi), "pts": torz(Q).tolist()}); we.Next()
        for k in range(len(L)):                       # chain check / orientation fix
            nxt = np.array(L[(k + 1) % len(L)]["pts"]); cur = np.array(L[k]["pts"])
            if min(np.linalg.norm(cur[-1] - nxt[0]), np.linalg.norm(cur[-1] - nxt[-1])) > 1e-5 and \
               min(np.linalg.norm(cur[0] - nxt[0]), np.linalg.norm(cur[0] - nxt[-1])) < 1e-5:
                L[k]["pts"] = cur[::-1].tolist()
        if L: loops.append(L)
    slits = []
    for j, F in enumerate(free):
        sc = BRepAlgoAPI_Section(F, pl); sc.Build()
        for e in items(sc.Shape(), TopAbs_EDGE):
            slits.append({"face": 100 + j, "pts": torz(sample(e)).tolist()})
    return {"phi": deg, "loops": loops, "slits": slits}

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="r-z sections of the domain solid at the given phi (deg) -> sections/")
    ap.add_argument("phi", type=float, nargs="+")
    ap.add_argument("--brep", default=SRC)
    a = ap.parse_args()
    SRC = a.brep
    os.makedirs("sections", exist_ok=True)
    solid, sf, fr = load()
    for deg in a.phi:
        degc = min(max(deg, 1e-4), 90 - 1e-4)   # exactly 0/90 coincides with the symmetry planes
        S = section(degc, solid, sf, fr); S["phi"] = deg
        json.dump(S, open(f"sections/sec_{deg:05.1f}.json", "w"))
        print(deg, "loops", [len(l) for l in S["loops"]], [[d["face"] for d in l] for l in S["loops"]],
              "slits", [(d["face"], len(d["pts"])) for d in S["slits"]], flush=True)
