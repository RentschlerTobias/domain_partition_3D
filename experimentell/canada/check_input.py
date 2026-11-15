"""Acceptance check of an AlgoHex input tet.vtk (feature curves and surface).

PASS requires: closed 2-manifold boundary, no dangling curve ends (valence-1 vertex that is not a
feature vertex), no turn > 30 deg at a valence-2 curve vertex that is not a feature vertex, no curve
segment shorter than 1/3 of the median edge length of its incident boundary triangles.
usage: check_input.py <tet.vtk> [labels.txt]
"""
import collections, sys
import numpy as np
import meshio

t = meshio.read(sys.argv[1]); P = t.points; key = next(iter(t.cell_data))
T = np.vstack([c.data for c in t.cells if c.type == "triangle"])
L = np.concatenate([np.asarray(d).ravel() for c, d in zip(t.cells, t.cell_data[key]) if c.type == "triangle"]).astype(int)
lines = np.vstack([c.data for c in t.cells if c.type == "line"])
FV = set(np.concatenate([c.data.ravel() for c in t.cells if c.type == "vertex"]).tolist())
names = {}
if len(sys.argv) > 2:
    names = {int(l.split()[0]): l.split()[1] for l in open(sys.argv[2])}
ok = True
ec = collections.Counter()
e2t = collections.defaultdict(list)
for i, tri in enumerate(T):
    for a, b in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
        k = tuple(sorted((int(a), int(b)))); ec[k] += 1; e2t[k].append(i)
mani = collections.Counter(ec.values())
print(f"boundary: {len(T)} triangles, edge multiplicity {dict(mani)}")
ok &= set(mani) == {2}
adj = collections.defaultdict(list)
for a, b in lines:
    adj[int(a)].append(int(b)); adj[int(b)].append(int(a))
ends = [v for v, nb in adj.items() if len(nb) == 1 and v not in FV]
# stricter: a curve end must sit where >= 3 surface labels meet (a real junction); an end that is only
# DECLARED a feature vertex hides a broken curve (the TE ridge once came in 9-15 pieces this way)
vlab = collections.defaultdict(set)
for tri, l in zip(T, L):
    for v in tri:
        vlab[int(v)].add(int(l))
loose = [v for v, nb in adj.items() if len(nb) == 1 and len(vlab[v]) < 3]
kinks = []
for v, nb in adj.items():
    if len(nb) == 2 and v not in FV:
        d1, d2 = P[v] - P[nb[0]], P[nb[1]] - P[v]
        ang = np.degrees(np.arccos(np.clip(d1 @ d2 / np.linalg.norm(d1) / np.linalg.norm(d2), -1, 1)))
        if ang > 30:
            kinks.append((ang, v))
tel = np.linalg.norm(P[T[:, [1, 2, 0]]] - P[T], axis=2)
short = []
for a, b in lines:
    k = tuple(sorted((int(a), int(b))))
    loc = np.median(tel[e2t[k]]) if e2t[k] else 0
    seg = np.linalg.norm(P[a] - P[b])
    if loc > 0 and seg < loc / 3:
        short.append((seg, loc, P[a]))
print(f"feature curves: {len(lines)} edges, {len(FV)} feature vertices; dangling ends {len(ends)}, "
      f"curve ends not at a surface junction {len(loose)}, "
      f"kinks > 30 deg without feature vertex {len(kinks)}, segments < 1/3 local size {len(short)}")
for ang, v in sorted(kinks, reverse=True)[:5]:
    print(f"   kink {ang:.0f} deg at r={np.hypot(*P[v][:2]):.3f} z={P[v][2]:.3f}")
for seg, loc, x in sorted(short)[:5]:
    print(f"   short segment {seg:.4f} (local {loc:.4f}) at r={np.hypot(*x[:2]):.3f} z={x[2]:.3f}")
ok &= not ends and not kinks and not short and not loose
# per label pair: components and dangling ends
pairs = collections.defaultdict(list)
for a, b in lines:
    k = tuple(sorted((int(a), int(b))))
    ls = tuple(sorted(set(L[e2t[k]].tolist())))
    pairs[ls].append((int(a), int(b)))
print("curves by label pair (edges, components):")
for ls, E in sorted(pairs.items(), key=lambda x: -len(x[1])):
    g = collections.defaultdict(list)
    for a, b in E:
        g[a].append(b); g[b].append(a)
    seen, nc = set(), 0
    for s0 in g:
        if s0 in seen:
            continue
        nc += 1; st = [s0]; seen.add(s0)
        while st:
            u = st.pop()
            for w in g[u]:
                if w not in seen:
                    seen.add(w); st.append(w)
    nm = " | ".join(names.get(l, str(l)) for l in ls)
    print(f"   {nm:45s} {len(E):5d} {nc:3d}")
print("PASS" if ok else "FAIL")
