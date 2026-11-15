"""Hausdorff part 2: CAD samples -> mesh boundary (exact point-triangle distance, boundary quads split in 2)."""
import json, numpy as np
from scipy.spatial import cKDTree
D = np.load("boundary_pts.npz"); O = np.load("hausdorff_occ.npz")
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
S = O["cad"]; _, nn = tree.query(S, k=24)
d = np.array([pt_tri(np.repeat(S[i][None], 24, 0), A[nn[i]], B[nn[i]], C[nn[i]]).min() for i in range(len(S))])
fwd = max(O["nodes"].max(), O["centres"].max())
i = d.argmax()
out = dict(mesh_to_cad_nodes_max=float(O["nodes"].max()), mesh_to_cad_centres_max=float(O["centres"].max()),
           mesh_to_cad_centres_p99=float(np.percentile(O["centres"], 99)),
           cad_to_mesh_max=float(d.max()), cad_to_mesh_p99=float(np.percentile(d, 99)), cad_to_mesh_mean=float(d.mean()),
           cad_to_mesh_argmax_point=S[i].tolist(), cad_to_mesh_argmax_face=int(O["cad_face"][i]),
           hausdorff=float(max(fwd, d.max())), n_cad_samples=int(len(S)))
j = O["centres"].argmax(); out["mesh_to_cad_centres_argmax_point"] = O["centres_pts"][j].tolist()
k = O["nodes"].argmax(); out["mesh_to_cad_nodes_argmax_point"] = O["nodes_pts"][k].tolist()
print(json.dumps(out, indent=1)); json.dump(out, open("hausdorff.json", "w"), indent=1)
