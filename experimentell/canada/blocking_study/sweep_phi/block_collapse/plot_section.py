"""phi=0 plane block edges of two block structures (z, r): baseline vs collapsed (reporting picture).
Reads the ASCII POLYDATA *_block_edges.vtk written by block_collapse.py."""
import sys, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
def read(fn):
    L = open(fn).read().split("\n"); i = next(k for k, l in enumerate(L) if l.startswith("POINTS")); n = int(L[i].split()[1])
    P = np.array([list(map(float, l.split())) for l in L[i + 1:i + 1 + n]])
    j = next(k for k, l in enumerate(L) if l.startswith("LINES")); m = int(L[j].split()[1])
    return P, [list(map(int, l.split()[1:])) for l in L[j + 1:j + 1 + m]]
fig, ax = plt.subplots(2, 1, figsize=(14, 9))
for a, (fn, title) in zip(ax, ((sys.argv[1], "sweep baseline (72 blocks)"), (sys.argv[2], "after block-level collapse (70 blocks)"))):
    P, lines = read(fn)
    for ln in lines:
        Q = P[ln]
        if np.all(np.abs(Q[:, 1]) < 1e-6):
            a.plot(Q[:, 2], Q[:, 0], "k-", lw=0.6)
    a.set_title(f"{title}: block edges in the phi=0 plane"); a.set_xlabel("z"); a.set_ylabel("r"); a.set_aspect("equal")
plt.tight_layout(); plt.savefig(sys.argv[3], dpi=110)
