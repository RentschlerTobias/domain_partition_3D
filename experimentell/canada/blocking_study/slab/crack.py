"""Turn internal feature surfaces (zero-thickness baffles, labels given) of an AlgoHex VTK into a crack:
duplicate every baffle vertex except those on the free (interior) edges of the baffle, so the baffle becomes
two coincident boundary sheets that close at the free edges (fluid angle 360 deg there).  Generic: the side
of each tet is found locally per vertex (components of the tets around the vertex separated by baffle faces)."""
import sys, collections, numpy as np, meshio
src, dst = sys.argv[1], sys.argv[2]; blabs = {int(x) for x in sys.argv[3].split(",")}
m = meshio.read(src)
P = m.points.copy()
cells = {c.type: (c.data, np.asarray(d).ravel()) for c, d in zip(m.cells, m.cell_data["color"])}
tets = cells["tetra"][0].copy(); T, TL = cells["triangle"]; E = cells["line"][0]; V = cells["vertex"][0].ravel()
isb = np.isin(TL, list(blabs)); BT = T[isb]; BL = TL[isb]; OT = T[~isb]; OL = TL[~isb]
def ekey(a, b): return (min(a, b), max(a, b))
bedge = collections.Counter(ekey(*e) for t in BT for e in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])))
oedge = {ekey(*e) for t in OT for e in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0]))}
tip = {v for e, k in bedge.items() if k == 1 and e not in oedge for v in e}
bverts = sorted(set(BT.ravel()) - tip)
print(f"[crack] baffle tris {len(BT)}, baffle verts {len(set(BT.ravel()))}, tip verts kept shared {len(tip)}, duplicated {len(bverts)}")
bface = {tuple(sorted(t)) for t in BT}
v2t = collections.defaultdict(list)
for i, t in enumerate(tets):
    for v in t: v2t[int(v)].append(i)
# normals of baffle triangles per vertex (for orienting the two sides)
nrm = np.cross(P[BT[:, 1]] - P[BT[:, 0]], P[BT[:, 2]] - P[BT[:, 0]])
vn = np.zeros_like(P)
for t, n in zip(BT, nrm):
    for v in t: vn[v] += n
newid = {}
P2 = list(P)
for v in bverts:
    ts = v2t[v]
    # union-find over the tets around v, joined through faces containing v that are not baffle faces
    par = {i: i for i in ts}
    def f(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    fm = collections.defaultdict(list)
    for i in ts:
        t = tets[i]
        for k in range(4):
            fc = tuple(sorted(int(x) for j, x in enumerate(t) if j != k))
            if v in fc: fm[fc].append(i)
    for fc, owners in fm.items():
        if len(owners) == 2 and fc not in bface: par[f(owners[0])] = f(owners[1])
    comps = collections.defaultdict(list)
    for i in ts: comps[f(i)].append(i)
    if len(comps) != 2:
        raise SystemExit(f"[crack] vertex {v}: {len(comps)} sides (expected 2)")
    n = vn[v] / np.linalg.norm(vn[v])
    PP = np.array(P2); up = max(comps.values(), key=lambda c: float(np.mean((PP[tets[c]].mean(1) - P[v]) @ n)))
    nid = len(P2); P2.append(P[v].copy()); newid[v] = nid
    for i in up:
        tets[i] = np.where(tets[i] == v, nid, tets[i])
P2 = np.array(P2)
# upper copy of each baffle triangle: the version whose vertices belong to upper tets
tri_up = []
tetface = collections.Counter()
for t in tets:
    for k in range(4):
        tetface[tuple(sorted(int(x) for j, x in enumerate(t) if j != k))] += 1
for t in BT:
    up = [newid.get(int(v), int(v)) for v in t]
    assert tetface[tuple(sorted(up))] == 1 and tetface[tuple(sorted(int(v) for v in t))] == 1
    tri_up.append(up)
# feature edges: copy edges that lie on duplicated vertices (junction curves) to the upper side too
E2 = [tuple(e) for e in E]
bset = set(bverts)
for a, b in E:
    if a in bset and b in bset: E2.append((newid[a], newid[b]))
    elif a in bset and b in tip or b in bset and a in tip:
        E2.append((newid.get(a, a), newid.get(b, b)))
V2 = list(V) + [newid[v] for v in V if v in newid]
# wall triangles touching the baffle junction: remap vertices of those on the upper side
OT2 = []
for t in OT:
    fc = tuple(sorted(int(x) for x in t))
    if tetface[fc] == 1: OT2.append(list(t)); continue
    alt = tuple(sorted(newid.get(int(x), int(x)) for x in t))
    # mixed: try remapping subsets
    import itertools
    done = False
    dv = [k for k, x in enumerate(t) if int(x) in newid]
    for r in range(1, len(dv) + 1):
        for sub in itertools.combinations(dv, r):
            tt = [newid[int(x)] if k in sub else int(x) for k, x in enumerate(t)]
            if tetface[tuple(sorted(tt))] == 1: OT2.append(tt); done = True; break
        if done: break
    assert done, t
cells_out = [(1, [int(v)], 1) for v in V2] + [(3, [int(a), int(b)], 1) for a, b in E2] + \
    [(5, [int(x) for x in t], int(l)) for t, l in zip(OT2, OL)] + [(5, [int(x) for x in t], int(l)) for t, l in zip(BT, BL)] + \
    [(5, [int(x) for x in t], int(l) + 100) for t, l in zip(tri_up, BL)] + [(10, [int(x) for x in t], 0) for t in tets]
# boundary sanity: every boundary tet face is a triangle
nb = sum(1 for k, c in tetface.items() if c == 1)
print(f"[crack] points {len(P)} -> {len(P2)}; boundary faces {nb}, triangles written {len(OT2) + 2 * len(BT)}")
with open(dst, "w") as f:
    f.write("# vtk DataFile Version 2.0\ncrack\nASCII\nDATASET UNSTRUCTURED_GRID\n")
    f.write(f"POINTS {len(P2)} double\n" + "".join(f"{p[0]:.12g} {p[1]:.12g} {p[2]:.12g}\n" for p in P2))
    f.write(f"CELLS {len(cells_out)} {sum(1 + len(c[1]) for c in cells_out)}\n")
    f.write("".join(f"{len(c[1])} " + " ".join(map(str, c[1])) + "\n" for c in cells_out))
    f.write(f"CELL_TYPES {len(cells_out)}\n" + "".join(f"{c[0]}\n" for c in cells_out))
    f.write(f"CELL_DATA {len(cells_out)}\nSCALARS color int 1\nLOOKUP_TABLE default\n" + "".join(f"{c[2]}\n" for c in cells_out))
