"""Where is a raw AlgoHex hex mesh broken? For one .ovm: hex count, inverted hexes (corner Jacobian < 0),
and the hex-mesh boundary faces that are NOT on the input surface (holes / gaps from flipped
parametrisation), clustered by position. Writes <out>.vtu with the hexes (field min_jac) and
<out>_holes.vtu with the off-surface boundary quads.

usage: hex_defects.py <hex.ovm> <tet.vtk> <labels.txt> <out_prefix>
"""
import collections, sys
import numpy as np
import meshio
from scipy.spatial import cKDTree
sys.path.insert(0, "/opt/stack/domain_partition_3D/experimentell/hex3d_algohex")
from ovm_io import read_ovm, ovm_to_cells

ovm, tvtk, labf, out = sys.argv[1:5]
P, E, F, C = read_ovm(ovm)
H, skipped = ovm_to_cells(P, E, F, C)
P = np.asarray(P, float); H = np.asarray(H, int)
# VTK hex: 0-1-2-3 bottom, 4-7 top; corner Jacobian = det of the three edge vectors at each corner
nb = [(1, 3, 4), (2, 0, 5), (3, 1, 6), (0, 2, 7), (7, 5, 0), (4, 6, 1), (5, 7, 2), (6, 4, 3)]
J = np.empty((len(H), 8))
for c, (a, b, d) in enumerate(nb):
    v = P[H[:, [a, b, d]]] - P[H[:, c]][:, None]
    n = np.linalg.norm(v, axis=2) + 1e-30
    J[:, c] = np.linalg.det(v / n[..., None])
mj = J.min(1)
# boundary quads of the hex mesh
QF = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
cnt = collections.Counter(); rep = {}
for h in H:
    for q in QF:
        k = tuple(sorted(h[list(q)])); cnt[k] += 1; rep[k] = h[list(q)]
bq = np.array([rep[k] for k, c in cnt.items() if c == 1])
nonman = sum(1 for c in cnt.values() if c > 2)
# distance of boundary quad centres to the input surface
t = meshio.read(tvtk); key = next(iter(t.cell_data))
T = np.vstack([c.data for c in t.cells if c.type == "triangle"])
L = np.concatenate([np.asarray(d).ravel() for c, d in zip(t.cells, t.cell_data[key]) if c.type == "triangle"]).astype(int)
names = {int(l.split()[0]): l.split()[1] for l in open(labf)}
tc = t.points[T].mean(1); tree = cKDTree(np.vstack([tc, t.points[np.unique(T)]]))
lab_pt = np.concatenate([L, np.full(len(np.unique(T)), -1)])
qc = P[bq].mean(1)
d, j = tree.query(qc)
h_loc = np.linalg.norm(P[bq[:, 0]] - P[bq[:, 2]], axis=1)
off = d > 0.5 * h_loc
print(f"{ovm}: {len(H)} hexes ({skipped} skipped cells), inverted (min corner Jacobian < 0): {(mj < 0).sum()}, "
      f"non-manifold faces {nonman}")
print(f"   boundary quads {len(bq)}, off the input surface (> h/2): {off.sum()}")
def where(X):
    r = np.hypot(X[:, 0], X[:, 1]); ph = np.degrees(np.arctan2(X[:, 1], X[:, 0]))
    return r, ph, X[:, 2]
if off.any():
    X = qc[off]
    # cluster off-surface quads (single linkage, 3 h)
    tr = cKDTree(X); lab = -np.ones(len(X), int); k = 0
    for i in range(len(X)):
        if lab[i] >= 0:
            continue
        st = [i]; lab[i] = k
        while st:
            u = st.pop()
            for w in tr.query_ball_point(X[u], 3 * np.median(h_loc)):
                if lab[w] < 0:
                    lab[w] = k; st.append(w)
        k += 1
    print(f"   {k} hole clusters (largest first):")
    for c in sorted(range(k), key=lambda c: -(lab == c).sum())[:12]:
        Y = X[lab == c]; r, ph, z = where(Y)
        dd, jj = tree.query(Y.mean(0))
        print(f"     {len(Y):5d} quads  r {r.min():.2f}..{r.max():.2f}  phi {ph.min():.0f}..{ph.max():.0f}  "
              f"z {z.min():.2f}..{z.max():.2f}  nearest input label {names.get(int(lab_pt[jj]), '?')}")
if (mj < 0).any():
    X = P[H[mj < 0]].mean(1); r, ph, z = where(X)
    print(f"   inverted hexes: r {np.percentile(r, [5, 50, 95]).round(2)}, z {np.percentile(z, [5, 50, 95]).round(2)}")
meshio.write(out + ".vtu", meshio.Mesh(P, [("hexahedron", H)], cell_data={"min_jac": [mj]}))
if off.any():
    meshio.write(out + "_holes.vtu", meshio.Mesh(P, [("quad", bq[off])], cell_data={"dist": [d[off]]}))
