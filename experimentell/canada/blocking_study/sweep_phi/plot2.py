import sys; sys.path.insert(0, ".")
from layout2d import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
secs = load_sections(); secs5 = [x for x in secs if abs(x[0] / 5 - round(x[0] / 5)) < 1e-9]; ml = mean_loops(secs5); ab = abstract(ml)
mg = json.load(open('layout_final.json'))['merges'] if 'final' in (sys.argv[1] if len(sys.argv)>1 else '') else []
lay = Layout(ml, ab, [tuple(x) for x in mg]); nb = lay.arrangement()
m = Mesh2D(lay, 0.1)
fig, axs = plt.subplots(1, 4, figsize=(28, 14))
plot_mesh(m, m.P, axs[0], title=f"parameter layout, {nb} blocks")
K = len(secs) - 1
for ax, k in zip(axs[1:], (0, K // 2, K)):
    phi, lp = secs[k]; X = m.solve(lp); sj = m.quality(X)
    plot_mesh(m, X, ax, sj, title=f"phi {phi}: inverted {(sj<=0).sum()}, min sJ {sj.min():.2f}")
plt.savefig(sys.argv[1] if len(sys.argv) > 1 else "layout_initial.png", dpi=80, bbox_inches="tight")
