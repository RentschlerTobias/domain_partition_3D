"""Generic CAD export for the block-level collapse (pythonOCC).

Reads a .brep / .step shape (solid(s) and/or free faces), enumerates its unique faces, edges and vertices
(TopExp maps, so the ids are those of the shape, not of any mesh), and writes
  faces : dense surface samples (point, face id)  -> boundary classification + reverse Hausdorff
  edges : dense polylines per edge (ds)            -> feature-curve classification / projection
  verts : vertex coordinates
  topology : edge->faces, vertex->edges
usage: LD_LIBRARY_PATH=/dtOO-install/lib /opt/venv-dtoo/bin/python cad_export.py shape.brep out.npz [ds_face] [ds_edge]
"""
import sys, json, math
import numpy as np
from OCC.Core.BRepTools import breptools
from OCC.Core.BRep import BRep_Builder, BRep_Tool
from OCC.Core.TopoDS import TopoDS_Shape, topods
from OCC.Core.TopExp import topexp
from OCC.Core.TopTools import TopTools_IndexedMapOfShape, TopTools_IndexedDataMapOfShapeListOfShape
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX, TopAbs_IN, TopAbs_ON
from OCC.Core.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface
from OCC.Core.BRepClass import BRepClass_FaceClassifier
from OCC.Core.gp import gp_Pnt2d
from OCC.Core.GProp import GProp_GProps
from OCC.Core.BRepGProp import brepgprop

src, out = sys.argv[1], sys.argv[2]
DSF = float(sys.argv[3]) if len(sys.argv) > 3 else 0.04
DSE = float(sys.argv[4]) if len(sys.argv) > 4 else 0.004
shape = TopoDS_Shape()
if src.lower().endswith((".stp", ".step")):
    from OCC.Core.STEPControl import STEPControl_Reader
    r = STEPControl_Reader(); r.ReadFile(src); r.TransferRoots(); shape = r.OneShape()
else:
    breptools.Read(shape, src, BRep_Builder())


def imap(t):
    m = TopTools_IndexedMapOfShape(); topexp.MapShapes(shape, t, m)
    return [m.FindKey(i) for i in range(1, m.Size() + 1)], m


faces, fmap = imap(TopAbs_FACE); edges, emap = imap(TopAbs_EDGE); verts, vmap = imap(TopAbs_VERTEX)
print(f"{len(faces)} faces, {len(edges)} edges, {len(verts)} vertices", flush=True)
e2f = TopTools_IndexedDataMapOfShapeListOfShape(); topexp.MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, e2f)
v2e = TopTools_IndexedDataMapOfShapeListOfShape(); topexp.MapShapesAndAncestors(shape, TopAbs_VERTEX, TopAbs_EDGE, v2e)


def anc(dm, sh, m):
    lst = dm.FindFromKey(sh); o = set()
    for x in lst:
        o.add(m.FindIndex(x) - 1)
    return sorted(o)


edge_faces = [anc(e2f, e, fmap) for e in edges]
vert_edges = [anc(v2e, v, emap) for v in verts]
V = np.array([BRep_Tool.Pnt(topods.Vertex(v)).Coord() for v in verts])

EP, EI = [], []
for ei, e in enumerate(edges):
    if BRep_Tool.Degenerated(topods.Edge(e)):
        continue
    a = BRepAdaptor_Curve(topods.Edge(e)); t = np.linspace(a.FirstParameter(), a.LastParameter(), 4001)
    P = np.array([a.Value(x).Coord() for x in t])
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
    n = max(2, int(math.ceil(s[-1] / DSE)) + 1); u = np.linspace(0, s[-1], n)
    EP.append(np.c_[[np.interp(u, s, P[:, k]) for k in range(3)]].T); EI.append(np.full(n, ei))

FP, FI = [], []
for fi, f in enumerate(faces):
    f = topods.Face(f); g = GProp_GProps(); brepgprop.SurfaceProperties(f, g); A = g.Mass()
    u0, u1, v0, v1 = breptools.UVBounds(f); s = BRep_Tool.Surface(f)
    # grid size from the physical extent of the uv box (sampled), target spacing DSF
    uu = np.linspace(u0, u1, 9); vv = np.linspace(v0, v1, 9)
    G = np.array([[s.Value(x, y).Coord() for y in vv] for x in uu])
    Lu = np.linalg.norm(np.diff(G, axis=0), axis=-1).sum(0).max(); Lv = np.linalg.norm(np.diff(G, axis=1), axis=-1).sum(1).max()
    nu = int(min(1500, max(8, Lu / DSF))); nv = int(min(1500, max(8, Lv / DSF)))
    cnt = 0
    for x in np.linspace(u0, u1, nu):
        for y in np.linspace(v0, v1, nv):
            if BRepClass_FaceClassifier(f, gp_Pnt2d(x, y), 1e-7).State() in (TopAbs_IN, TopAbs_ON):
                FP.append(s.Value(x, y).Coord()); FI.append(fi); cnt += 1
    print(f"face {fi}: area {A:.3f}, {nu}x{nv} uv, {cnt} samples, edges {[i for i, ef in enumerate(edge_faces) if fi in ef]}", flush=True)

np.savez(out, face_pts=np.array(FP), face_id=np.array(FI), edge_pts=np.vstack(EP), edge_id=np.concatenate(EI),
         verts=V, topo=json.dumps(dict(edge_faces=edge_faces, vert_edges=vert_edges, src=src)))
print("written", out, len(FP), "face samples", sum(len(p) for p in EP), "edge samples")
