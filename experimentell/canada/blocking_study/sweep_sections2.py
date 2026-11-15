import sys, math, json
exec(open('sweep_sections.py').read().split('for deg in')[0])
from OCC.Core.BRepAdaptor import BRepAdaptor_Curve
from OCC.Core.BRepExtrema import BRepExtrema_DistShapeShape
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
lab = {d['face']: d['label'] for d in json.load(open('/work/canada/cad/domain_blk_faces.json'))}
dfaces = items(s, TopAbs_FACE)
def on_face(pt):
    v = BRepBuilderAPI_MakeVertex(pt).Vertex(); best = (1e9, None)
    for i, fc in enumerate(dfaces):
        d = BRepExtrema_DistShapeShape(v, fc); d.Perform()
        if d.Value() < best[0]: best = (d.Value(), i)
    return lab.get(best[1], best[1])
for deg in [float(a) for a in sys.argv[2:]]:
    p = math.radians(deg)
    ax = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(-math.sin(p), math.cos(p), 0), gp_Dir(math.cos(p), math.sin(p), 0))
    f = BRepBuilderAPI_MakeFace(gp_Pln(ax), 0.0, 20.0, -10.0, 20.0).Face()
    c = BRepAlgoAPI_Common(s, f); c.Build()
    fc = items(c.Shape(), TopAbs_FACE)[0]
    out = []
    for w in items(fc, TopAbs_WIRE):
        labs = []
        for e in items(w, TopAbs_EDGE):
            a = BRepAdaptor_Curve(topods.Edge(e)); m = a.Value(0.5 * (a.FirstParameter() + a.LastParameter()))
            labs.append(on_face(m))
        out.append(labs)
    print(f"phi {deg}: " + " || ".join(",".join(x) for x in out), flush=True)
