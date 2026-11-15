"""2D block layout of a planar section: gmsh quasi-structured quad mesh (cross field, algorithm 11), then
the 2D base complex: separatrices traced from irregular vertices and from label/corner points; patches =
blocks of the 2D layout (= blocks of the swept 3D structure).  Not minimal, an upper estimate."""
import sys, json, collections, numpy as np, gmsh
brep, h = sys.argv[1], float(sys.argv[2])
labels = json.load(open(brep.replace(".brep", "_edges.json")))
gmsh.initialize(); gmsh.option.setNumber("General.Verbosity", 2)
gmsh.model.occ.importShapes(brep); gmsh.model.occ.synchronize()
gmsh.option.setNumber("Mesh.MeshSizeMax", h); gmsh.option.setNumber("Mesh.MeshSizeMin", h / 20)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 12)
gmsh.option.setNumber("Mesh.Algorithm", 11)
gmsh.option.setNumber("Mesh.RecombineAll", 1)
gmsh.model.mesh.generate(2)
nt, xyz, _ = gmsh.model.mesh.getNodes(); P = xyz.reshape(-1, 3)[:, [0, 2]]; idx = {int(k): i for i, k in enumerate(nt)}
Q = []
for typ, _, nn in zip(*gmsh.model.mesh.getElements(2)):
    if typ == 3: Q += [[idx[int(v)] for v in nn[i:i + 4]] for i in range(0, len(nn), 4)]
    elif typ == 2: print("[quad2d] WARNING: triangles left", len(nn) // 3)
Q = np.array(Q)
# corner points of the CAD section (gmsh points) and curve membership of boundary nodes
corner = set()
for _, t in gmsh.model.getEntities(0):
    n_, _, _ = gmsh.model.mesh.getNodes(0, t); corner |= {idx[int(x)] for x in n_}
gmsh.write(brep.replace(".brep", "_quad.msh")); gmsh.finalize()
E = collections.defaultdict(list)
for qi, q in enumerate(Q):
    for k in range(4):
        a, b = q[k], q[(k + 1) % 4]; E[(min(a, b), max(a, b))].append(qi)
bnd = {e for e, l in E.items() if len(l) == 1}
bv = {v for e in bnd for v in e}
val = collections.Counter(v for e in E for v in e)      # edge valence
nq = collections.Counter(v for q in Q for v in q)
irreg = {v for v in val if (v not in bv and val[v] != 4)}
# boundary: at a CAD corner the block layout needs a separatrix if the number of quads != round(angle/90)
src = set(irreg) | {v for v in corner if v in bv}
adj = collections.defaultdict(set)
for a, b in E: adj[a].add(b); adj[b].add(a)
def straight_next(prev, cur):
    """continue a quad-edge line through a regular vertex: the edge 'opposite' to (cur,prev)"""
    if cur in bv or val[cur] != 4: return None
    quads_prev = set(E[(min(prev, cur), max(prev, cur))])
    for n in adj[cur]:
        if n == prev: continue
        if not (set(E[(min(n, cur), max(n, cur))]) & quads_prev): return n
    return None
sep = set()
for v in src:
    for n in adj[v]:
        e = (min(v, n), max(v, n))
        if e in bnd: continue
        prev, cur = v, n; sep.add(e)
        while True:
            nx = straight_next(prev, cur)
            if nx is None: break
            e = (min(cur, nx), max(cur, nx))
            if e in sep: break
            sep.add(e); prev, cur = cur, nx
cut = sep | bnd
par = list(range(len(Q)))
def f(x):
    while par[x] != x: par[x] = par[par[x]]; x = par[x]
    return x
for e, l in E.items():
    if len(l) == 2 and e not in cut: par[f(l[0])] = f(l[1])
blocks = len({f(i) for i in range(len(Q))})
print(json.dumps(dict(quads=len(Q), irregular_interior=len(irreg),
                      valences=dict(collections.Counter(val[v] for v in irreg)),
                      corners=len(corner), blocks_2d=blocks)))
