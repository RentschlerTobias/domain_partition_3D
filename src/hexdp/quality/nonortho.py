import sys, numpy as np, meshio
from collections import defaultdict
FACES = np.array([(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)])
m = meshio.read(sys.argv[1]); P = m.points; H = m.cells_dict['hexahedron']; B = np.asarray(m.cell_data['block_id'][0]).ravel()
C = P[H].mean(1)
F = H[:, FACES].reshape(-1, 4); own = np.repeat(np.arange(len(H)), 6)
key = np.sort(F, axis=1); order = np.lexsort(key.T[::-1]); ks = key[order]
same = np.all(ks[1:] == ks[:-1], axis=1); i = np.where(same)[0]
a, b = order[i], order[i + 1]
fa = F[a]; X = P[fa]; n = np.cross(X[:, 2] - X[:, 0], X[:, 3] - X[:, 1]); n /= np.linalg.norm(n, axis=1)[:, None]
d = C[own[b]] - C[own[a]]; d /= np.linalg.norm(d, axis=1)[:, None]
ang = np.degrees(np.arccos(np.clip(np.abs((n * d).sum(1)), 0, 1)))
inter = B[own[a]] != B[own[b]]
r = np.hypot(X.mean(1)[:, 0], X.mean(1)[:, 1])
print('internal faces', len(ang), 'max %.1f' % ang.max(), ' >65:', (ang > 65).sum(), ' >50:', (ang > 50).sum())
for t in (65, 50):
    s = ang > t
    print(f'  >{t}: on block interfaces {int((s & inter).sum())}, inside blocks {int((s & ~inter).sum())}; r of them: hub-near(<0.55) {int((s & (r<0.55)).sum())}, shroud-near(>1.85) {int((s & (r>1.85)).sum())}')
w = np.argsort(-ang)[:5]
print('worst:', [(round(float(ang[k]), 1), int(B[own[a[k]]]), int(B[own[b[k]]]), round(float(r[k]), 3), round(float(X[k].mean(0)[2]), 3)) for k in w])
