import sys, json, collections, numpy as np
sys.path.insert(0, "/opt/stack/dp3d_keep_bl/experimentell/hex3d_algohex")
import ovm_io, base_complex as bc
from scipy.spatial import cKDTree
n = sys.argv[1]
A = np.load("in_a/tet_arrays.npz"); info = json.load(open("in_a/input_info.json"))
Q, T, TL, E = A["P"], A["tris"], A["tlab"], A["edges"]
lab = {int(l.split()[0]): l.split()[1] for l in open("in_a/labels.txt")}
P, e, f, poly = ovm_io.read_ovm(f"in_a/hex_n{n}.ovm"); H, _ = ovm_io.ovm_to_cells(P, e, f, poly); H = np.asarray(H)
f2h, e2h = bc.build_topology(H)
S = bc.singular_edges(H, P, f2h, e2h, 40.0)
M = np.array([(P[a] + P[b]) / 2 for a, b in S]); val = np.array([len(set(e2h[s])) for s in S])
h = np.mean(np.linalg.norm(P[H[:, 0]] - P[H[:, 1]], axis=1))
# distances: to baffle free edges (feature lines touching baffle tris only), to feature lines, to boundary tris by label
bafl = set(info["baffle_labels"])
vt = collections.defaultdict(set)
for t, l in zip(T, TL):
    for v in t: vt[int(v)].add(int(l))
fe_baffle = np.array([(Q[a] + Q[b]) / 2 for a, b in E if vt[int(a)] & bafl and vt[int(b)] & bafl])
fe_all = np.array([(Q[a] + Q[b]) / 2 for a, b in E])
d_bfe = cKDTree(fe_baffle).query(M)[0]; d_fe = cKDTree(fe_all).query(M)[0]
tc = Q[T].mean(1); tri_tree = cKDTree(tc); d_s, it = tri_tree.query(M)
r = np.hypot(M[:, 0], M[:, 1])
cat = collections.Counter()
for i in range(len(S)):
    if r[i] < 0.6 and (M[i, 2] < -2.2 or M[i, 2] > 2.28): c = "axis region (r<0.6, outside hub)"
    elif d_bfe[i] < 1.5 * h: c = "near baffle edges/junctions (<1.5h)"
    elif d_fe[i] < 1.5 * h: c = f"near other feature curves (<1.5h)"
    elif d_s[i] < 1.5 * h: c = f"near surface {lab[int(TL[it[i]])]} (<1.5h)"
    else: c = "interior"
    cat[c] += 1
print(f"n={n}: {len(S)} singular edges, valence {dict(collections.Counter(val))}, h {h:.3f}")
for k, v in cat.most_common(): print(f"   {v:4d}  {k}")
# connected singular curves
g = collections.defaultdict(list)
for a, b in S: g[a].append(b); g[b].append(a)
seen = set(); comps = []
for s in g:
    if s in seen: continue
    st = [s]; seen.add(s); c = []
    while st:
        u = st.pop(); c.append(u)
        for w in g[u]:
            if w not in seen: seen.add(w); st.append(w)
    comps.append(c)
print("   singular curves:", len(comps), "lengths(edges):", sorted([len(c) - 1 for c in comps], reverse=True)[:15])
for c in sorted(comps, key=len, reverse=True):
    X = P[c]; rr = np.hypot(X[:, 0], X[:, 1]); ph = np.degrees(np.arctan2(X[:, 1], X[:, 0]))
    vv = collections.Counter(int(val[i]) for i, (a, b) in enumerate(S) if a in set(c))
    print(f"   curve {len(c)-1:3d} edges val {dict(vv)}: z {X[:,2].min():.2f}..{X[:,2].max():.2f}, r {rr.min():.2f}..{rr.max():.2f}, phi {ph.min():.0f}..{ph.max():.0f}")
