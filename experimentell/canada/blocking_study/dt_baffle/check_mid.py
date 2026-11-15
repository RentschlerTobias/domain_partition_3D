"""distance of baffle mesh nodes to the two skins of their wall: (d_a - d_b) should be ~0, d_a+d_b ~ thickness"""
import sys, json; sys.path.insert(0, "/work/canada/cad")
import numpy as np
from simple_envelope import faces
from OCC.Core.BRepExtrema import BRepExtrema_DistShapeShape
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
from OCC.Core.gp import gp_Pnt
A = np.load("in_a/tet_arrays.npz"); P, T, L = A["P"], A["tris"], A["tlab"]
lab = {l.split()[1]: int(l.split()[0]) for l in open("in_a/labels.txt")}
info = json.load(open("parts/parts.json"))
d = lambda p, f: BRepExtrema_DistShapeShape(BRepBuilderAPI_MakeVertex(gp_Pnt(*p)).Vertex(), f).Value()
rng = np.random.default_rng(0); out = {}
for name, m in info["mid"].items():
    key = "BAFFLE_LIP" if m["lip"] else "BAFFLE_" + name.upper()
    V = np.unique(T[L == lab[key]]); V = rng.choice(V, min(150, len(V)), replace=False)
    a, b = m["skins"]
    da = np.array([d(P[v], faces[a]) for v in V]); db = np.array([d(P[v], faces[b]) for v in V])
    out[key] = dict(asym_max=float(np.abs(da - db).max()), asym_median=float(np.median(np.abs(da - db))),
                    thickness_median=float(np.median(da + db)))
    print(key, {k: round(v, 4) for k, v in out[key].items()})
json.dump(out, open("in_a/mid_check.json", "w"), indent=1)
