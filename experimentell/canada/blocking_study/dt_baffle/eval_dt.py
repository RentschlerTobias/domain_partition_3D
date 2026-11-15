"""Evaluate an AlgoHex hex mesh of the draft-tube baffle input and derive its block structure.
 - hexes, inverted (scaled Jacobian <= 0), min sJ, singular edges (base_complex, sharp 40 deg)
 - boundary fit: distance of boundary-quad vertices / centroids to the non-baffle input triangles
 - baffle coverage: interior hex faces whose 4 vertices + centroid lie on a baffle (point-triangle distance
   < tol), area / baffle area, per baffle label
 - split along the baffle faces (tip vertices on free baffle edges stay shared), then base complex with the
   singular edges + the 360-deg crack-tip/concave boundary edges -> blocks (split_blocks.py method, curved)
 - writes <ovm>_eval.json, <ovm>_split.vtk (cell data: block id, sJ), <ovm>_baffle_quads.vtk
usage: eval_dt.py <hex.ovm> <input dir with tet_arrays.npz, input_info.json>
"""
import sys, json, collections
import numpy as np
sys.path.insert(0, "/opt/stack/dp3d_keep_bl/experimentell/hex3d_algohex")
import ovm_io, base_complex as bc
from scipy.spatial import cKDTree

ovm, ind = sys.argv[1], sys.argv[2]
A = np.load(f"{ind}/tet_arrays.npz"); info = json.load(open(f"{ind}/input_info.json"))
Q, T, TL = A["P"], A["tris"], A["tlab"]
P, e, f, poly = ovm_io.read_ovm(ovm); H, skipped = ovm_io.ovm_to_cells(P, e, f, poly)
H = np.asarray(H)
sj = np.array([ovm_io.scaled_jacobian(P, c) for c in H])
h_avg = float(np.mean(np.linalg.norm(P[H[:, 0]] - P[H[:, 1]], axis=1)))
res = dict(hexes=int(len(H)), skipped=int(skipped), inverted=int((sj <= 0).sum()), min_sj=float(sj.min()),
           sj_p1=float(np.percentile(sj, 1)), mean_edge=h_avg)


def pt_tri_dist(X, tri_idx):
    """exact distance of points X to the triangle set (nearest 12 by centroid)"""
    C = Q[T[tri_idx]].mean(1); tree = cKDTree(C)
    _, nn = tree.query(X, k=min(12, len(C)))
    best = np.full(len(X), np.inf)
    for k in range(nn.shape[1]):
        tri = Q[T[tri_idx[nn[:, k]]]]
        best = np.minimum(best, _pt_tri(X, tri[:, 0], tri[:, 1], tri[:, 2]))
    return best


def _pt_tri(p, a, b, c):
    # Ericson, Real-Time Collision Detection 5.1.5 (vectorised)
    ab, ac, ap = b - a, c - a, p - a
    d1, d2 = (ab * ap).sum(1), (ac * ap).sum(1)
    bp = p - b; d3, d4 = (ab * bp).sum(1), (ac * bp).sum(1)
    cp = p - c; d5, d6 = (ab * cp).sum(1), (ac * cp).sum(1)
    va = d3 * d6 - d5 * d4; vb = d5 * d2 - d1 * d6; vc = d1 * d4 - d3 * d2
    den = va + vb + vc; den[den == 0] = 1e-300
    v = vb / den; w = vc / den
    q = a + ab * v[:, None] + ac * w[:, None]
    m = (d1 <= 0) & (d2 <= 0); q[m] = a[m]
    m = (d3 >= 0) & (d4 <= d3); q[m] = b[m]
    m = (d6 >= 0) & (d5 <= d6); q[m] = c[m]
    m = (vc <= 0) & (d1 >= 0) & (d3 <= 0); t = d1 / np.where(d1 - d3 == 0, 1, d1 - d3); q[m] = (a + ab * t[:, None])[m]
    m = (vb <= 0) & (d2 >= 0) & (d6 <= 0); t = d2 / np.where(d2 - d6 == 0, 1, d2 - d6); q[m] = (a + ac * t[:, None])[m]
    m = (va <= 0) & ((d4 - d3) >= 0) & ((d5 - d6) >= 0)
    t = (d4 - d3) / np.where((d4 - d3) + (d5 - d6) == 0, 1, (d4 - d3) + (d5 - d6)); q[m] = (b + (c - b) * t[:, None])[m]
    return np.linalg.norm(p - q, axis=1)


FACES = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
fc = collections.defaultdict(list)
for i, h in enumerate(H):
    for fa in FACES:
        fc[tuple(sorted(h[list(fa)]))].append((i, h[list(fa)]))
bq = [own[0][1] for k, own in fc.items() if len(own) == 1]
iq = {k: own[0][1] for k, own in fc.items() if len(own) == 2}
bq = np.array(bq)
baf = np.isin(TL, info["baffle_labels"])
nb_idx = np.where(~baf)[0]
dv = pt_tri_dist(P[np.unique(bq)], nb_idx); dc = pt_tri_dist(P[bq].mean(1), nb_idx)
res.update(boundary_quads=int(len(bq)), boundary_vertex_dist_max=float(dv.max()),
           boundary_centroid_dist_max=float(dc.max()), boundary_centroid_dist_p99=float(np.percentile(dc, 99)),
           boundary_quads_off_surface=int((dc > 0.25 * h_avg).sum()))
# baffle quads
tol = 0.02 * h_avg
keys = list(iq); V = np.array([iq[k] for k in keys])
baffle_q = {}
quad_area = lambda X: 0.5 * (np.linalg.norm(np.cross(X[:, 1] - X[:, 0], X[:, 2] - X[:, 0]), axis=1)
                             + np.linalg.norm(np.cross(X[:, 2] - X[:, 0], X[:, 3] - X[:, 0]), axis=1))
for lab in info["baffle_labels"]:
    idx = np.where(TL == lab)[0]
    lo, hi = Q[np.unique(T[idx])].min(0) - 0.1, Q[np.unique(T[idx])].max(0) + 0.1
    cand = np.where(((P[V] >= lo) & (P[V] <= hi)).all(axis=(1, 2)))[0]
    X = P[V[cand]]
    dvv = pt_tri_dist(X.reshape(-1, 3), idx).reshape(-1, 4).max(1)
    dcc = pt_tri_dist(X.mean(1), idx)
    on = cand[(dvv < tol) & (dcc < 0.1 * h_avg)]   # vertices on the baffle; centroid may sag (planar quad on a curved sheet)
    for j in on:
        baffle_q[keys[j]] = lab
    tri_area = 0.5 * np.linalg.norm(np.cross(Q[T[idx, 1]] - Q[T[idx, 0]], Q[T[idx, 2]] - Q[T[idx, 0]]), axis=1).sum()
    res[f"baffle{lab}_quads"] = int(len(on))
    res[f"baffle{lab}_coverage"] = float(quad_area(P[V[on]]).sum() / tri_area) if len(on) else 0.0
f2h, e2h = bc.build_topology(H)
sing0 = bc.singular_edges(H, P, f2h, e2h, 40.0)
res["singular_edges_raw"] = int(len(sing0))
cut0 = bc.sheet_faces(H, f2h, set(sing0)); blocks0 = bc.blocks_from_cut(H, f2h, cut0)
res["raw_blocks_unsplit"] = int(len(blocks0))
# ---- split along baffle quads
bnd_edges = set()
qedges = lambda q: [tuple(sorted((int(q[i]), int(q[(i + 1) % 4])))) for i in range(4)]
for q in bq:
    bnd_edges |= set(qedges(q))
be = collections.Counter(x for k in baffle_q for x in qedges(iq[k]))
tip = {v for x, c in be.items() if c == 1 and x not in bnd_edges for v in x}
verts = {int(v) for k in baffle_q for v in k} - tip
v2h = collections.defaultdict(list)
for i, h in enumerate(H):
    for v in h:
        v2h[int(v)].append(i)
Pl = list(P); H2 = H.copy(); bad_split = 0
bset = set(baffle_q)
for v in sorted(verts):
    hs = v2h[v]; par = {i: i for i in hs}
    def fnd(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for i in hs:
        for fa in FACES:
            k = tuple(sorted(H[i][list(fa)]))
            if v in k and k not in bset:
                for j, _ in fc[k]:
                    if j != i and j in par:
                        par[fnd(i)] = fnd(j)
    comps = collections.defaultdict(list)
    for i in hs:
        comps[fnd(i)].append(i)
    if len(comps) < 2:
        bad_split += 1; continue
    for c in list(comps.values())[1:]:
        nid = len(Pl); Pl.append(P[v].copy())
        for i in c:
            H2[i] = np.where(H2[i] == v, nid, H2[i])
P2 = np.array(Pl)
res.update(baffle_quads=len(baffle_q), tip_vertices=len(tip), split_vertices=len(verts), split_failures=bad_split)
f2h, e2h = bc.build_topology(H2)
sing = set(bc.singular_edges(H2, P2, f2h, e2h, 40.0))
bf = {k for k, hs in f2h.items() if len(hs) == 1}
bE = set()
for c in H2:
    for fa in bc.HEX_FACES:
        if frozenset(int(c[i]) for i in fa) in bf:
            l = [int(c[i]) for i in fa]
            bE |= {(min(l[i], l[(i + 1) % 4]), max(l[i], l[(i + 1) % 4])) for i in range(4)}
extra = {x for x in bE if len(set(e2h[x])) >= 3}
res["singular_edges_split"] = int(len(sing)); res["concave_or_tip_boundary_edges"] = int(len(extra - sing))
sing |= extra
cut = bc.sheet_faces(H2, f2h, sing); blocks = bc.blocks_from_cut(H2, f2h, cut)
sizes = sorted((len(b) for b in blocks), reverse=True)
res.update(blocks=len(blocks), block_sizes_top=sizes[:30], blocks_le_4_cells=int(sum(s <= 4 for s in sizes)))
bid = np.zeros(len(H2), int)
for i, b in enumerate(blocks):
    for h in b:
        bid[h] = i
base = ovm.replace(".ovm", "")
np.save(base + "_split_blockid.npy", bid)
with open(base + "_split.vtk", "w") as fo:
    fo.write("# vtk DataFile Version 2.0\nsplit\nASCII\nDATASET UNSTRUCTURED_GRID\n")
    fo.write(f"POINTS {len(P2)} double\n" + "".join(f"{p[0]} {p[1]} {p[2]}\n" for p in P2))
    fo.write(f"CELLS {len(H2)} {9 * len(H2)}\n" + "".join("8 " + " ".join(map(str, h)) + "\n" for h in H2))
    fo.write(f"CELL_TYPES {len(H2)}\n" + "12\n" * len(H2))
    fo.write(f"CELL_DATA {len(H2)}\nSCALARS block int 1\nLOOKUP_TABLE default\n" + "".join(f"{b}\n" for b in bid))
    fo.write("SCALARS sJ double 1\nLOOKUP_TABLE default\n" + "".join(f"{s}\n" for s in sj))
with open(base + "_baffle_quads.vtk", "w") as fo:
    ks = list(baffle_q)
    fo.write("# vtk DataFile Version 2.0\nbaffle quads\nASCII\nDATASET UNSTRUCTURED_GRID\n")
    fo.write(f"POINTS {len(P)} double\n" + "".join(f"{p[0]} {p[1]} {p[2]}\n" for p in P))
    fo.write(f"CELLS {len(ks)} {5 * len(ks)}\n" + "".join("4 " + " ".join(map(str, iq[k])) + "\n" for k in ks))
    fo.write(f"CELL_TYPES {len(ks)}\n" + "9\n" * len(ks))
    fo.write(f"CELL_DATA {len(ks)}\nSCALARS label int 1\nLOOKUP_TABLE default\n" + "".join(f"{baffle_q[k]}\n" for k in ks))
# singular edges as lines
with open(base + "_singular.vtk", "w") as fo:
    S = sorted(sing0)
    fo.write("# vtk DataFile Version 2.0\nsingular edges\nASCII\nDATASET UNSTRUCTURED_GRID\n")
    fo.write(f"POINTS {len(P)} double\n" + "".join(f"{p[0]} {p[1]} {p[2]}\n" for p in P))
    fo.write(f"CELLS {len(S)} {3 * len(S)}\n" + "".join(f"2 {a} {b}\n" for a, b in S))
    fo.write(f"CELL_TYPES {len(S)}\n" + "3\n" * len(S))
print(json.dumps(res, indent=1))
json.dump(res, open(base + "_eval.json", "w"), indent=1)
