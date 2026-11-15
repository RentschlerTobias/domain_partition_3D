"""Diffuser (blk) domain of canada as one CAD solid, quarter model, with the thin walls replaced by their
envelopes:  outer body (farMid + inlet/outlet planes + symmetry planes)
            minus hub chain (gv hub bulb, gv/ru hole between the two interfaces, dt hub; bulb and hole
                  are solids of revolution from the exact dtOO meridians (section with the plane y = 0),
                  the dt hub is closed from its own dtOO faces since its end is not axisymmetric)
            minus envelopes (lip = mantle + ring 1, ring 2, ring 3).
Every face of the result gets a label from the source face it lies on (names as in the dtOO physical
export, envelopes ENV_*). Writes domain_blk.brep + domain_blk_faces.json (label + interior point per face).

run with LD_LIBRARY_PATH=/dtOO-install/lib /opt/venv-dtoo/bin/python domain_blk.py
"""
import json, sys
import numpy as np
sys.argv = sys.argv[:1]
import simple_envelope as se
from simple_envelope import faces, BRep_Tool, gp_Pnt, count, volume
from OCC.Core.gp import gp_Pln, gp_Dir, gp_Ax1
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Section, BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCC.Core.BRepBuilderAPI import (BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace, BRepBuilderAPI_Sewing,
                                     BRepBuilderAPI_MakeSolid)
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeRevol
from OCC.Core.ShapeAnalysis import ShapeAnalysis_FreeBounds
from OCC.Core.TopTools import TopTools_HSequenceOfShape, TopTools_ListOfShape
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_EDGE, TopAbs_WIRE, TopAbs_FACE, TopAbs_SHELL, TopAbs_SOLID, TopAbs_IN, TopAbs_ON
from OCC.Core.TopoDS import topods, TopoDS_Shape
from OCC.Core.BRepTools import breptools
from OCC.Core.BRep import BRep_Builder
from OCC.Core.BRepCheck import BRepCheck_Analyzer
from OCC.Core.BRepLib import BRepLib_FindSurface
from OCC.Core.Geom import Geom_Plane
from OCC.Core.ShapeFix import ShapeFix_Solid
from OCC.Core.BRepExtrema import BRepExtrema_DistShapeShape
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeVertex, BRepBuilderAPI_MakeWire
from OCC.Core.TopTools import TopTools_IndexedDataMapOfShapeListOfShape
from OCC.Core.TopExp import topexp
from OCC.Core.BRepClass import BRepClass_FaceClassifier
from OCC.Core.gp import gp_Pnt2d

F = dict(bulb=13, hub0=14, hub1=15, hubEnd=16, far=17, if_gv=18, if_dt=19)
LABEL = {"bulb": "IN_HUB", "hub0": "DT_HUB", "hub1": "DT_HUB", "hubEnd": "DT_HUB", "far": "INGVRUDTOUT_FAR",
         "if_gv": "IN_OUTLET", "if_dt": "DT_INLET", "inlet": "IN_INLET", "outlet": "OUT_OUTLET",
         "lip": "ENV_LIP", "ring2": "ENV_RING2", "ring3": "ENV_RING3"}


def edges_of(shape):
    out = []; e = TopExp_Explorer(shape, TopAbs_EDGE)
    while e.More():
        out.append(topods.Edge(e.Current())); e.Next()
    return out


def ends(e):
    c, f0, f1 = BRep_Tool.Curve(e)
    return np.array(c.Value(f0).Coord()), np.array(c.Value(f1).Coord())


def meridian(face):
    """edges of the section of a (revolved) face with the half plane y = 0, x >= 0"""
    sec = BRepAlgoAPI_Section(face, gp_Pln(gp_Pnt(0, 0, 0), gp_Dir(0, 1, 0)))
    sec.Build()
    out = []
    for e in edges_of(sec.Shape()):
        c, f0, f1 = BRep_Tool.Curve(e)
        if c.Value(0.5 * (f0 + f1)).X() > 1e-9:
            out.append(e)
    return out


def chain_ends(es):
    """the two free end points (in the y = 0 plane, as (x, z)) of a chain of edges"""
    P = [p for e in es for p in ends(e)]
    free = [p for p in P if sum(np.linalg.norm(p - q) < 1e-6 for q in P) == 1]
    return free


def connect(seq):
    """edges (any order / orientation) -> one closed wire, chained by end points (tolerance 1e-6)"""
    es = [topods.Edge(seq.Value(k)) for k in range(1, seq.Length() + 1)]
    order = [es.pop(0)]; tip = ends(order[0])[1]; start = ends(order[0])[0]
    while es:
        k = next((k for k, e in enumerate(es) if min(np.linalg.norm(p - tip) for p in ends(e)) < 1e-6), None)
        if k is None:
            raise RuntimeError(f"profile not connected at {tip}")
        e = es.pop(k); a, b = ends(e)
        tip = b if np.linalg.norm(a - tip) < 1e-6 else a
        order.append(e)
    mw = BRepBuilderAPI_MakeWire()
    for e in order:
        mw.Add(e)
    out = TopTools_HSequenceOfShape(); out.Append(mw.Wire())
    return out


def line(a, b):
    return BRepBuilderAPI_MakeEdge(gp_Pnt(*a), gp_Pnt(*b)).Edge()


def revolve_profile(es, name):
    seq = TopTools_HSequenceOfShape()
    for e in es:
        seq.Append(e)
    wires = connect(seq)
    assert wires.Length() == 1, f"{name}: profile not one wire ({wires.Length()})"
    w = topods.Wire(wires.Value(1))
    face = BRepBuilderAPI_MakeFace(w, True).Face()
    rev = BRepPrimAPI_MakeRevol(face, gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), np.pi / 2).Shape()
    print(f"   {name}: profile {len(es)} edges, solid volume {volume(rev):.4f}", flush=True)
    return rev


def axis_body(fs, name):
    """solid from quarter faces around the z axis: sew them, close every free chain that is not on a
    symmetry plane with a planar cap through the axis, close the symmetry planes (they share the axis)"""
    sew = BRepBuilderAPI_Sewing(1e-6)
    for f in fs:
        sew.Add(f)
    sew.Perform(); sh = sew.SewedShape()
    m = TopTools_IndexedDataMapOfShapeListOfShape(); topexp.MapShapesAndAncestors(sh, TopAbs_EDGE, TopAbs_FACE, m)
    free = [topods.Edge(m.FindKey(i)) for i in range(1, m.Size() + 1) if m.FindFromIndex(i).Size() == 1]
    def in_plane(e, k):
        c, f0, f1 = BRep_Tool.Curve(e)
        return all(abs(c.Value(f0 + t * (f1 - f0)).Coord()[k]) < 1e-6 for t in np.linspace(0, 1, 7))
    sym = {1: [e for e in free if in_plane(e, 1)], 0: [e for e in free if in_plane(e, 0)]}
    other = [e for e in free if not in_plane(e, 0) and not in_plane(e, 1)]
    caps = []
    for e in other:                                   # e.g. the arc at the interface plane
        p, q = ends(e)
        a = np.array([0, 0, p[2]])
        lp, lq = line(p, a), line(a, q)
        seq = TopTools_HSequenceOfShape()
        for x in (e, lp, lq):
            seq.Append(x)
        caps.append(BRepBuilderAPI_MakeFace(topods.Wire(connect(seq).Value(1)), True).Face())
        for x in (lp, lq):
            sym[1 if in_plane(x, 1) else 0].append(x)
    axis_z = sorted({round(p[2], 9) for k in (0, 1) for e in sym[k] for p in ends(e) if np.hypot(p[0], p[1]) < 1e-6})
    axis = line(np.array([0, 0, axis_z[0]]), np.array([0, 0, axis_z[-1]]))
    sew2 = BRepBuilderAPI_Sewing(1e-6); sew2.Add(sh)
    for c in caps:
        sew2.Add(c)
    for k in (0, 1):
        seq = TopTools_HSequenceOfShape()
        for x in sym[k] + [axis]:
            seq.Append(x)
        sew2.Add(BRepBuilderAPI_MakeFace(topods.Wire(connect(seq).Value(1)), True).Face())
    sew2.Perform()
    ms = BRepBuilderAPI_MakeSolid()
    e = TopExp_Explorer(sew2.SewedShape(), TopAbs_SHELL)
    while e.More():
        ms.Add(topods.Shell(e.Current())); e.Next()
    fx = ShapeFix_Solid(ms.Solid()); fx.Perform(); body = fx.Solid()
    print(f"   {name}: {len(fs)} faces + {len(caps)} caps + 2 symmetry planes, valid {BRepCheck_Analyzer(body).IsValid()}, "
          f"volume {volume(body):.4f}", flush=True)
    return body


def hub_chain():
    O = lambda z: np.array([0.0, 0.0, z])
    # gv hub bulb (tip on the axis) closed downstream by a radial line
    mb = meridian(faces[F["bulb"]]); a, b = sorted(chain_ends(mb), key=lambda p: p[2])
    bulb = revolve_profile(mb + [line(b, O(b[2])), line(O(b[2]), O(a[2]))] + ([line(O(a[2]), a)] if a[0] > 1e-6 else []), "bulb")
    # gv/ru hole: interface gv (upstream), cylinder at the interface's outer radius (lies inside the lip
    # envelope), interface dt (downstream), hub core on the axis
    m18 = meridian(faces[F["if_gv"]]); A, B = sorted(chain_ends(m18), key=lambda p: p[0])
    m19 = meridian(faces[F["if_dt"]]); D, C = sorted(chain_ends(m19), key=lambda p: p[0])
    Bz = np.array([B[0], 0, C[2]])
    hole = revolve_profile(m18 + m19 + [line(B, Bz), line(Bz, C), line(D, O(D[2])), line(O(D[2]), O(A[2])),
                                        line(O(A[2]), A)], "hole")
    # dt hub: hub_0, hub_1, hubEnd from the dtOO faces themselves (the hub end is not axisymmetric)
    hub = axis_body([faces[F[k]] for k in ("hub0", "hub1", "hubEnd")], "dt hub")
    s = BRepAlgoAPI_Fuse(BRepAlgoAPI_Fuse(bulb, hole).Shape(), hub).Shape()
    return s


def outer_body():
    """farMid closed by the inlet / outlet planes and the two symmetry planes (which share the axis)"""
    far = faces[F["far"]]
    es = edges_of(far)
    zs = [p[2] for e in es for p in ends(e)]; zmin, zmax = min(zs), max(zs)
    on = lambda e, k, v: all(abs(p[k] - v) < 1e-6 for p in ends(e)) and \
        abs(BRep_Tool.Curve(e)[0].Value(0.5 * sum(BRep_Tool.Curve(e)[1:])).Coord()[k] - v) < 1e-6
    e_in = next(e for e in es if on(e, 2, zmin)); e_out = next(e for e in es if on(e, 2, zmax))
    e_y0 = next(e for e in es if on(e, 1, 0.0)); e_x0 = next(e for e in es if on(e, 0, 0.0))
    o_in, o_out = np.array([0, 0, zmin]), np.array([0, 0, zmax])
    px = lambda e: next(p for p in ends(e) if abs(p[1]) < 1e-6)       # end point on y = 0
    py = lambda e: next(p for p in ends(e) if abs(p[0]) < 1e-6)       # end point on x = 0
    Lix, Liy = line(o_in, px(e_in)), line(o_in, py(e_in))
    Lox, Loy = line(o_out, px(e_out)), line(o_out, py(e_out))
    axis = line(o_in, o_out)
    def planar(edges):
        seq = TopTools_HSequenceOfShape()
        for x in edges:
            seq.Append(x)
        return BRepBuilderAPI_MakeFace(topods.Wire(connect(seq).Value(1)), True).Face()
    caps = {"inlet": planar([e_in, Lix, Liy]), "outlet": planar([e_out, Lox, Loy])}
    sym = [planar([e_y0, Lix, axis, Lox]), planar([e_x0, Liy, axis, Loy])]
    sew = BRepBuilderAPI_Sewing(1e-6)
    for f in [far, caps["inlet"], caps["outlet"]] + sym:
        sew.Add(f)
    sew.Perform()
    ms = BRepBuilderAPI_MakeSolid()
    e = TopExp_Explorer(sew.SewedShape(), TopAbs_SHELL)
    while e.More():
        ms.Add(topods.Shell(e.Current())); e.Next()
    fx = ShapeFix_Solid(ms.Solid()); fx.Perform()
    body = fx.Solid()
    print(f"   outer body: faces {count(body, TopAbs_FACE)}, valid {BRepCheck_Analyzer(body).IsValid()}, "
          f"volume {volume(body):.3f}", flush=True)
    return body, caps


def interior_point(face):
    """a point inside the face (not on its boundary)"""
    u0, u1, v0, v1 = breptools.UVBounds(face)
    srf = BRep_Tool.Surface(face)
    rng = np.random.default_rng(0)
    for k in range(400):
        u, v = (u0 + (u1 - u0) * (0.5 if k == 0 else rng.random()), v0 + (v1 - v0) * (0.5 if k == 0 else rng.random()))
        if BRepClass_FaceClassifier(face, gp_Pnt2d(u, v), 1e-7).State() == TopAbs_IN:
            return np.array(srf.Value(u, v).Coord())
    raise RuntimeError("no interior point")


def main():
    body, caps = outer_body()
    hubs = hub_chain()
    envs = {}
    for n in ("lip", "ring2", "ring3"):
        s = TopoDS_Shape(); breptools.Read(s, f"/work/canada/cad/{n}_simple_envelope.brep", BRep_Builder()); envs[n] = s
    args = TopTools_ListOfShape(); args.Append(body)
    tools = TopTools_ListOfShape()
    for s in [hubs] + list(envs.values()):
        tools.Append(s)
    cut = BRepAlgoAPI_Cut(); cut.SetArguments(args); cut.SetTools(tools); cut.SetFuzzyValue(1e-6); cut.Build()
    dom = cut.Shape()
    sols = []; e = TopExp_Explorer(dom, TopAbs_SOLID)
    while e.More():
        sols.append(topods.Solid(e.Current())); e.Next()
    print(f"[domain] solids {len(sols)}, faces {count(dom, TopAbs_FACE)}, valid {BRepCheck_Analyzer(dom).IsValid()}, "
          f"volume {volume(dom):.3f} (outer {volume(body):.3f})", flush=True)
    dom = sols[0] if len(sols) == 1 else dom
    # labels: source face nearest to an interior point of each result face
    src = {k: faces[i] for k, i in F.items()}
    src.update(caps)
    for n, s in envs.items():
        src[n] = s
    out = []
    fe = TopExp_Explorer(dom, TopAbs_FACE); k = 0
    from collections import Counter
    cnt = Counter()
    while fe.More():
        f = topods.Face(fe.Current()); p = interior_point(f)
        if abs(p[1]) < 1e-7:
            lab, d = "INGVRUDTOUT_SYMNY", 0.0
        elif abs(p[0]) < 1e-7:
            lab, d = "INGVRUDTOUT_SYMNX", 0.0
        else:
            v = BRepBuilderAPI_MakeVertex(gp_Pnt(*p)).Vertex()
            dd = {n: BRepExtrema_DistShapeShape(v, s).Value() for n, s in src.items()}
            n = min(dd, key=dd.get); lab, d = LABEL[n], dd[n]
        if d > 1e-5:
            print(f"   WARNING face {k}: nearest source {lab} at distance {d:.2e}")
        out.append({"face": k, "label": lab, "point": p.tolist(), "dist": d}); cnt[lab] += 1
        k += 1; fe.Next()
    print("[domain] faces per label:", dict(cnt), flush=True)
    breptools.Write(dom, "/work/canada/cad/domain_blk.brep")
    json.dump(out, open("/work/canada/cad/domain_blk_faces.json", "w"), indent=1)


if __name__ == "__main__":
    main()
