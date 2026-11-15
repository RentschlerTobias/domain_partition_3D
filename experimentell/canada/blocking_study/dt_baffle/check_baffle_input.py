"""Check of an AlgoHex input with interior baffles (tet_arrays.npz from assemble_mesh.py).
 - every boundary face of the tet mesh (1 tet) is a labelled triangle with a non-baffle label
 - every baffle triangle is an INTERIOR face (2 tets) -> baffle kept interior, not a crack
 - no duplicated nodes (min node distance), every tet positively oriented / non-degenerate
 - boundary triangles form a closed 2-manifold (edge multiplicity 2 among non-baffle triangles)
 - feature curves: every curve end (valence != 2) is a feature vertex; kinks > 30 deg at valence-2
   vertices are feature vertices; baffle free edges (1 baffle triangle, not on the boundary) are feature edges
usage: check_baffle_input.py <dir>
"""
import sys, json, collections
import numpy as np
from scipy.spatial import cKDTree

d = sys.argv[1]
A = np.load(f"{d}/tet_arrays.npz"); info = json.load(open(f"{d}/input_info.json"))
P, Tt, T, L, E, FV = A["P"], A["tets"], A["tris"], A["tlab"], A["edges"], set(A["fv"].tolist())
baf = np.isin(L, info["baffle_labels"])
res = {}
fc = collections.Counter()
for t in Tt:
    for f in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)):
        fc[tuple(sorted(t[list(f)]))] += 1
bfaces = {k for k, c in fc.items() if c == 1}
tri_keys = [tuple(sorted(t)) for t in T]
res["tets"] = int(len(Tt)); res["nodes"] = int(len(P))
res["tet_faces_multiplicity"] = dict(collections.Counter(fc.values()))
res["boundary_faces"] = len(bfaces)
res["boundary_faces_without_triangle"] = len(bfaces - set(k for k, b in zip(tri_keys, baf) if not b))
res["nonbaffle_triangles_not_on_boundary"] = int(sum(1 for k, b in zip(tri_keys, baf) if not b and k not in bfaces))
res["baffle_triangles"] = int(baf.sum())
res["baffle_triangles_interior(2 tets)"] = int(sum(1 for k, b in zip(tri_keys, baf) if b and fc[k] == 2))
res["baffle_triangles_on_boundary"] = int(sum(1 for k, b in zip(tri_keys, baf) if b and fc[k] == 1))
dd, _ = cKDTree(P).query(P, k=2)
res["min_node_distance"] = float(dd[:, 1].min())
res["duplicate_nodes(<1e-9)"] = int((dd[:, 1] < 1e-9).sum())
v = P[Tt]; vol = np.einsum('ij,ij->i', np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0]), v[:, 3] - v[:, 0]) / 6
res["tet_volume_sign"] = {"pos": int((vol > 0).sum()), "neg": int((vol < 0).sum()), "zero": int((vol == 0).sum())}
res["domain_volume"] = float(np.abs(vol).sum())
ec = collections.Counter()
for t, b in zip(T, baf):
    if b: continue
    for x, y in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
        ec[(min(x, y), max(x, y))] += 1
res["boundary_edge_multiplicity"] = dict(collections.Counter(ec.values()))
# baffle free edges
bec = collections.Counter()
for t in T[baf]:
    for x, y in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
        bec[(min(x, y), max(x, y))] += 1
fe = {(min(x, y), max(x, y)) for x, y in E}
free = {k for k, c in bec.items() if c == 1 and k not in ec}
res["baffle_free_edges"] = len(free)
res["baffle_free_edges_not_feature"] = len(free - fe)
res["baffle_boundary_junction_edges"] = sum(1 for k in bec if k in ec)
res["baffle_junction_edges_not_feature"] = sum(1 for k in bec if k in ec and k not in fe)
val = collections.Counter(int(x) for e in E for x in e)
inc = collections.defaultdict(list)
for x, y in E:
    inc[int(x)].append(int(y)); inc[int(y)].append(int(x))
res["curve_ends_not_feature_vertex"] = sum(1 for v_, k in val.items() if k != 2 and v_ not in FV)
kink = 0
for v_, nb in inc.items():
    if len(nb) == 2 and v_ not in FV:
        d1, d2 = P[v_] - P[nb[0]], P[nb[1]] - P[v_]
        if np.degrees(np.arccos(np.clip(d1 @ d2 / np.linalg.norm(d1) / np.linalg.norm(d2), -1, 1))) > 30:
            kink += 1
res["kinks_gt30_without_feature_vertex"] = kink
res["feature_edges"] = int(len(E)); res["feature_vertices"] = len(FV)
res["curve_vertex_valence"] = dict(sorted(collections.Counter(val.values()).items()))
# per baffle: area and free-edge length
for lab in info["baffle_labels"]:
    S = T[L == lab]
    ar = 0.5 * np.linalg.norm(np.cross(P[S[:, 1]] - P[S[:, 0]], P[S[:, 2]] - P[S[:, 0]]), axis=1).sum()
    res[f"baffle_{lab}_area"] = float(ar)
ok = (res["boundary_faces_without_triangle"] == 0 and res["nonbaffle_triangles_not_on_boundary"] == 0
      and res["baffle_triangles_on_boundary"] == 0 and res["baffle_triangles"] == res["baffle_triangles_interior(2 tets)"]
      and res["duplicate_nodes(<1e-9)"] == 0 and set(res["boundary_edge_multiplicity"]) == {2}
      and res["baffle_free_edges_not_feature"] == 0 and res["baffle_junction_edges_not_feature"] == 0
      and res["curve_ends_not_feature_vertex"] == 0 and res["kinks_gt30_without_feature_vertex"] == 0
      and res["tet_volume_sign"]["zero"] == 0 and min(res["tet_volume_sign"]["pos"], res["tet_volume_sign"]["neg"]) == 0)
res["PASS"] = bool(ok)
print(json.dumps(res, indent=1))
json.dump(res, open(f"{d}/input_check.json", "w"), indent=1)
