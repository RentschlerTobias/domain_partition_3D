import sys, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, collections
from matplotlib.collections import PolyCollection, LineCollection
sys.path.insert(0, "/work/canada/blocking_study/sweep_phi/cgrid")
from cgrid2d import quad_sj
def bedges(P, Q):
    ec = collections.defaultdict(list)
    for q in Q:
        for k in range(4):
            a, b = q[k], q[(k+1)%4]; ec[(min(a,b), max(a,b))].append(q[4])
    return np.array([(P[a], P[b]) for (a, b), l in ec.items() if len(l) == 1 or l[0] != l[1]])
def show(ax, P, Q, xl, yl, title):
    sj = quad_sj(P, Q)
    ax.add_collection(PolyCollection(P[Q[:, :4]], facecolors="none", edgecolors="0.7", linewidths=0.2))
    if (sj <= 0).any(): ax.add_collection(PolyCollection(P[Q[sj <= 0, :4]], facecolors="red", edgecolors="red", linewidths=0.3))
    ax.add_collection(LineCollection(bedges(P, Q), colors="k", linewidths=0.8))
    ax.set_xlim(*xl); ax.set_ylim(*yl); ax.set_aspect("equal"); ax.set_title(title, fontsize=9)
if __name__ == "__main__":
    views = [((-0.2, 11.2), (-4.2, 10.2)), ((1.6, 3.2), (-0.2, 2.8)), ((1.7, 2.5), (-1.9, -1.4)), ((2.5, 4.2), (1.8, 4.7)), ((3.9, 6.5), (3.9, 6.6)), ((2.2, 2.7), (2.1, 2.6))]
    for f in sys.argv[1:]:
        D = np.load(f); P, Q = D["P"], D["Q"]
        fig, axs = plt.subplots(2, 3, figsize=(24, 16))
        for ax, (xl, yl) in zip(axs.ravel(), views):
            show(ax, P, Q, xl, yl, f)
        plt.savefig(f.replace(".npz", ".png"), dpi=70, bbox_inches="tight")
