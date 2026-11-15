"""Are the narrow slots between non-adjacent walls filled by the hex mesh, and with how many layers?
Slot samples: for boundary vertices of the given labels whose nearest non-adjacent other-label point is
closer than --max-gap, the midpoint of the two. A sample is 'covered' if it lies inside a hex (5-tet
split); layers across = gap / local hex edge length.

usage: slot_check.py <hexraw.vtu> <tet.vtk> <labels.txt> [--labels ENV_LIP,ENV_RING2,ENV_RING3] [--max-gap 0.15]
"""
import argparse, collections
import numpy as np
import meshio
from scipy.spatial import cKDTree

ap = argparse.ArgumentParser()
ap.add_argument("hex"); ap.add_argument("tet"); ap.add_argument("labels")
ap.add_argument("--labels-use", default="ENV_LIP,ENV_RING2,ENV_RING3")
ap.add_argument("--max-gap", type=float, default=0.15)
a = ap.parse_args()
names = {int(l.split()[0]): l.split()[1] for l in open(a.labels)}
lid = {v: k for k, v in names.items()}
use = [lid[n] for n in a.labels_use.split(",")]
t = meshio.read(a.tet); key = next(iter(t.cell_data))
T = np.vstack([c.data for c in t.cells if c.type == "triangle"])
L = np.concatenate([np.asarray(d).ravel() for c, d in zip(t.cells, t.cell_data[key]) if c.type == "triangle"]).astype(int)
X = t.points
vl = collections.defaultdict(set)
for tr, l in zip(T, L):
    for v in tr:
        vl[int(v)].add(int(l))
samples, gaps, pair = [], [], []
for A in use:
    va = np.array([v for v, s in vl.items() if A in s and len(s) == 1])
    for B in use:
        if B <= A:
            continue
        vb = np.array([v for v, s in vl.items() if B in s and len(s) == 1])
        d, j = cKDTree(X[vb]).query(X[va])
        m = d < a.max_gap
        samples.append(0.5 * (X[va[m]] + X[vb[j[m]]])); gaps.append(d[m]); pair += [f"{names[A]}|{names[B]}"] * m.sum()
S = np.vstack(samples); G = np.concatenate(gaps); pair = np.array(pair)
h = meshio.read(a.hex); H = h.cells_dict["hexahedron"]; P = h.points
C = P[H].mean(1)
edge = np.linalg.norm(P[H[:, 0]] - P[H[:, 1]], axis=1)
TETS = [(0, 1, 3, 4), (1, 2, 3, 6), (1, 4, 5, 6), (3, 4, 6, 7), (1, 3, 4, 6)]
def inside(p, hx):
    for tt in TETS:
        q = P[hx[list(tt)]]
        M = np.column_stack([q[1] - q[0], q[2] - q[0], q[3] - q[0]])
        try:
            l = np.linalg.solve(M, p - q[0])
        except np.linalg.LinAlgError:
            continue
        if l.min() >= -1e-9 and l.sum() <= 1 + 1e-9:
            return True
    return False
tree = cKDTree(C)
cov = np.zeros(len(S), bool); hloc = np.zeros(len(S))
for i, p in enumerate(S):
    _, cand = tree.query(p, k=12)
    hloc[i] = np.median(edge[cand])
    cov[i] = any(inside(p, H[c]) for c in cand)
print(f"{a.hex}: {len(H)} hexes")
for pr in sorted(set(pair)):
    m = pair == pr
    print(f"   {pr:22s} slot samples {m.sum():5d}  gap {G[m].min():.3f}..{G[m].max():.3f}  covered {100 * cov[m].mean():5.1f} %  "
          f"layers across (gap/h) median {np.median(G[m] / hloc[m]):.2f}, min {np.min(G[m] / hloc[m]):.2f}")
