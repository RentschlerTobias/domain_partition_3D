"""Hausdorff part 2: CAD samples -> mesh boundary (exact point-triangle distance, boundary quads split in 2)."""
import json, numpy as np
from scipy.spatial import cKDTree
D = np.load("boundary_pts.npz"); O = np.load("hausdorff_g_occ.npz")
P = D["pts"]; Q = D["quads"]
T = np.vstack([Q[:, [0, 1, 2]], Q[:, [0, 2, 3]]]); A, B, C = P[T[:, 0]], P[T[:, 1]], P[T[:, 2]]
tree = cKDTree((A + B + C) / 3)
def pt_tri(p, a, b, c):
    # Ericson, closest point on triangle (vectorised over triangles)
    ab, ac, ap = b - a, c - a, p - a
    d1 = (ab * ap).sum(-1); d2 = (ac * ap).sum(-1)
    bp = p - b; d3 = (ab * bp).sum(-1); d4 = (ac * bp).sum(-1)
    cp = p - c; d5 = (ab * cp).sum(-1); d6 = (ac * cp).sum(-1)
    va = d3 * d6 - d5 * d4; vb = d5 * d2 - d1 * d6; vc = d1 * d4 - d3 * d2
    den = va + vb + vc; den = np.where(np.abs(den) < 1e-300, 1e-300, den)
    v = vb / den; w = vc / den
    res = a + ab * v[..., None] + ac * w[..., None]
    # region tests
    m = (d1 <= 0) & (d2 <= 0); res[m] = a[m]
    m = (d3 >= 0) & (d4 <= d3); res[m] = b[m]
    m = (d6 >= 0) & (d5 <= d6); res[m] = c[m]
    m = (vc <= 0) & (d1 >= 0) & (d3 <= 0); t = d1 / np.where(d1 - d3 == 0, 1, d1 - d3); res[m] = (a + ab * t[..., None])[m]
    m = (vb <= 0) & (d2 >= 0) & (d6 <= 0); t = d2 / np.where(d2 - d6 == 0, 1, d2 - d6); res[m] = (a + ac * t[..., None])[m]
    m = (va <= 0) & ((d4 - d3) >= 0) & ((d5 - d6) >= 0); t = (d4 - d3) / np.where((d4 - d3) + (d5 - d6) == 0, 1, (d4 - d3) + (d5 - d6))
    res[m] = (b + (c - b) * t[..., None])[m]
    return np.linalg.norm(p - res, axis=-1)
S = O["cad"]
# symmetry-plane faces come from the mid-surface domain and still cover the inside of the real walls: drop plane
# samples that lie inside a real wall section (rings: closed polygon; lip: open chain closed by its end points)
import sys; sys.path.insert(0, "/work/canada/blocking_study/sweep_phi")
from layout2d import point_in_poly
Wsec = {ph: [np.array(w["pts"]) for w in json.load(open(f"/work/canada/blocking_study/sweep_phi/cgrid/walls/wall_{ph:05.1f}.json"))["walls"]] for ph in (0.0, 90.0)}
keep = np.ones(len(S), bool)
for i, q in enumerate(S):
    if O["cad_grp"][i] != "plane":
        continue
    ph = 0.0 if abs(q[1]) < 1e-7 else (90.0 if abs(q[0]) < 1e-7 else None)
    if ph is None:
        continue
    rz = (np.hypot(q[0], q[1]), q[2])
    if any(point_in_poly(rz, Wp) for Wp in Wsec[ph]):
        keep[i] = False
print("plane samples inside the real walls dropped:", int((~keep).sum()))
S = S[keep]; O = dict(O); O["cad_grp"] = O["cad_grp"][keep]; O["cad_face"] = O["cad_face"][keep]
_, nn = tree.query(S, k=24)
d = np.array([pt_tri(np.repeat(S[i][None], 24, 0), A[nn[i]], B[nn[i]], C[nn[i]]).min() for i in range(len(S))])
fwd = max(O["nodes"].max(), O["centres"].max())
i = d.argmax()
out = dict(mesh_to_cad_nodes_max=float(O["nodes"].max()), mesh_to_cad_centres_max=float(O["centres"].max()),
           mesh_to_cad_centres_p99=float(np.percentile(O["centres"], 99)),
           cad_to_mesh_max=float(d.max()), cad_to_mesh_p99=float(np.percentile(d, 99)), cad_to_mesh_mean=float(d.mean()),
           cad_to_mesh_argmax_point=S[i].tolist(), cad_to_mesh_argmax_face=int(O["cad_face"][i]),
           hausdorff=float(max(fwd, d.max())), n_cad_samples=int(len(S)))
g = O["cad_grp"]; out["cad_to_mesh_wall_max"] = float(d[g == "wall"].max()); out["cad_to_mesh_wall_p99"] = float(np.percentile(d[g == "wall"], 99))
iw = np.where(g == "wall")[0][np.argmax(d[g == "wall"])]; out["cad_to_mesh_wall_argmax_point"] = S[iw].tolist()
cg = O["centres_grp"]; out["mesh_to_cad_centres_wall_max"] = float(O["centres"][cg == "wall"].max())
out["mesh_to_cad_nodes_wall_max"] = float(O["nodes"][O["nodes_grp"] == "wall"].max())
jw = np.where(cg == "wall")[0][np.argmax(O["centres"][cg == "wall"])]; out["mesh_to_cad_centres_wall_argmax"] = O["centres_pts"][jw].tolist()
j = O["centres"].argmax(); out["mesh_to_cad_centres_argmax_point"] = O["centres_pts"][j].tolist()
k = O["nodes"].argmax(); out["mesh_to_cad_nodes_argmax_point"] = O["nodes_pts"][k].tolist()
print(json.dumps(out, indent=1)); json.dump(out, open("hausdorff_g.json", "w"), indent=1)
