"""Block structure of a hex mesh that conforms to internal zero-thickness baffles: split the mesh along the
hex faces lying on the baffle (crack, tip vertices shared), then the base complex (sheets from singular
edges incl. the 360-deg crack-tip edges) -> blocks; also beam_collapse-ready VTK/OVM-free output."""
import sys, json, collections, numpy as np, meshio
sys.path.insert(0, "/opt/stack/dp3d_keep_bl/experimentell/hex3d_algohex")
import ovm_io, base_complex as bc
from scipy.spatial import cKDTree
ovm, vtk = sys.argv[1], sys.argv[2]
P, e, f, poly = ovm_io.read_ovm(ovm); H, _ = ovm_io.ovm_to_cells(P, e, f, poly); H0 = H.copy(); H = H.copy()
m = meshio.read(vtk); cd = {c.type: (c.data, np.asarray(d).ravel()) for c, d in zip(m.cells, m.cell_data["color"])}
T, TL = cd["triangle"]; Q = m.points
FACES = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
fc = collections.defaultdict(list)
for i, h in enumerate(H):
    for fa in FACES: fc[tuple(sorted(h[list(fa)]))].append((i, h[list(fa)]))
h_avg = np.mean(np.linalg.norm(P[H[:, 0]] - P[H[:, 1]], axis=1))
baffle = set()
for lab in (7, 8):
    sel = TL == lab; X = Q[np.unique(T[sel])]; c = X.mean(0); n = np.linalg.svd(X - c)[2][2]
    tc = cKDTree(Q[T[sel]].mean(1))
    for k, own in fc.items():
        if len(own) != 2: continue
        V = P[list(k)]
        if (np.abs((V - c) @ n) < 0.05 * h_avg).all() and tc.query(V.mean(0))[0] < 0.5 * h_avg: baffle.add(k)
bnd_edges = set()
def qedges(q): return [tuple(sorted((q[i], q[(i + 1) % 4]))) for i in range(4)]
for k, own in fc.items():
    if len(own) == 1: bnd_edges |= set(qedges(own[0][1]))
be = collections.Counter(x for k in baffle for x in qedges(fc[k][0][1]))
tip = {v for x, c in be.items() if c == 1 and x not in bnd_edges for v in x}
verts = {v for k in baffle for v in k} - tip
v2h = collections.defaultdict(list)
for i, h in enumerate(H):
    for v in h: v2h[int(v)].append(i)
Pl = list(P)
for v in sorted(verts):
    hs = v2h[v]; par = {i: i for i in hs}
    def fnd(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for k, own in fc.items():
        pass
    # join hexes around v through shared faces that contain v and are not baffle faces
    for i in hs:
        for fa in FACES:
            k = tuple(sorted(H0[i][list(fa)]))
            if v in k and k not in baffle:
                for j, _ in fc[k]:
                    if j != i and j in par: par[fnd(i)] = fnd(j)
    comps = collections.defaultdict(list)
    for i in hs: comps[fnd(i)].append(i)
    if len(comps) != 2:
        print("vertex", v, P[v], "comps", [len(c) for c in comps.values()], "hexes", len(hs), "tip?", v in tip); continue
    nid = len(Pl); Pl.append(P[v].copy())
    for i in list(comps.values())[1]: H[i] = np.where(H[i] == v, nid, H[i])
P2 = np.array(Pl)
f2h, e2h = bc.build_topology(H); sing = set(bc.singular_edges(H, P2, f2h, e2h, 40.0))
bf = {k for k, hs in f2h.items() if len(hs) == 1}
bE = set()
for k in bf:
    pass
for hi, c in enumerate(H):
    for fa in bc.HEX_FACES:
        if frozenset(int(c[i]) for i in fa) in bf:
            l = [int(c[i]) for i in fa]
            bE |= {(min(l[i], l[(i + 1) % 4]), max(l[i], l[(i + 1) % 4])) for i in range(4)}
extra = {e for e in bE if len(set(e2h[e])) >= 3}
print(f"[split] boundary edges with hex valence >= 3 (concave / crack tip): {len(extra)}")
sing |= extra
cut = bc.sheet_faces(H, f2h, sing); blocks = bc.blocks_from_cut(H, f2h, cut)
sizes = sorted((len(b) for b in blocks), reverse=True)
res = dict(baffle_quads=len(baffle), tip_vertices=len(tip), split_vertices=len(verts), singular_edges=len(sing),
           blocks=len(blocks), block_sizes=sizes)
print(json.dumps(res))
ovm_io.write_hex_vtk(ovm.replace(".ovm", "_split.vtk"), P2, H)
bid = np.zeros(len(H), int)
for i, b in enumerate(blocks):
    for h in b: bid[h] = i
np.save(ovm.replace(".ovm", "_split_blockid.npy"), bid)
json.dump(res, open(ovm.replace(".ovm", "_split.json"), "w"))
