import sys, json, collections, numpy as np
sys.argv = ["x", "in_a/hex_n20000.ovm", "in_a"]
src = open("eval_dt.py").read().split("f2h, e2h = bc.build_topology(H)\nsing0")[0]
src = src.replace("print(json.dumps", "#")
exec(src)
lab = {int(l.split()[0]): l.split()[1] for l in open("in_a/labels.txt")}
# vertex-on-baffle distances for interior quads near baffles
for L_ in info["baffle_labels"]:
    idx = np.where(TL == L_)[0]
    lo, hi = Q[np.unique(T[idx])].min(0) - 0.1, Q[np.unique(T[idx])].max(0) + 0.1
    cand = np.where(((P[V] >= lo) & (P[V] <= hi)).all(axis=(1, 2)))[0]
    X = P[V[cand]]
    dv = pt_tri_dist(X.reshape(-1, 3), idx).reshape(-1, 4).max(1); dc = pt_tri_dist(X.mean(1), idx)
    for tv in (1e-4, 1e-3, tol, 0.02):
        on = cand[(dv < tv)]
        ar = quad_area(P[V[on]]).sum()
        tri_area = 0.5 * np.linalg.norm(np.cross(Q[T[idx, 1]] - Q[T[idx, 0]], Q[T[idx, 2]] - Q[T[idx, 0]]), axis=1).sum()
        print(lab[L_], f"vertex tol {tv:.4f}: quads {len(on)} cov {ar/tri_area:.3f}; centroid dist p50/p95 {np.median(dc[dv<tv]) if len(on) else 0:.4f}/{np.percentile(dc[dv<tv],95) if len(on) else 0:.4f}")
# baffle vertices: distance of every input baffle vertex area covered? -> fraction of baffle triangle centroids within 0.5h of a covered quad
# inverted hex locations
bad = np.where(sj <= 0)[0]
C = P[H[bad]].mean(1)
nb_tree = {}
for name, sel in (("baffle", baf), ("boundary", ~baf)):
    idx = np.where(sel)[0]; d = pt_tri_dist(C, idx); nb_tree[name] = d
print("inverted:", len(bad))
for i, c in enumerate(C):
    near_l = None
    print(f"  r {np.hypot(*c[:2]):.3f} z {c[2]:.3f} phi {np.degrees(np.arctan2(c[1], c[0])):.1f} sJ {sj[bad[i]]:.2f} d_baffle {nb_tree['baffle'][i]:.3f} d_bnd {nb_tree['boundary'][i]:.3f}")
