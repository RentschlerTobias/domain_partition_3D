import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, sys
sys.path.insert(0, "/work/canada/blocking_study/sweep_phi/cgrid")
from plot_c import show
fig, axs = plt.subplots(3, 4, figsize=(28, 24))
for row, phi in enumerate((0, 45, 90)):
    D = np.load(f"sec2d_phi{phi:02d}.npz"); P, Q = D["P"], D["Q"]
    W = [np.array(w["pts"]) for w in json.load(open(f"walls/wall_{float(phi):05.1f}.json"))["walls"]]
    lip, r2, r3 = W[0], W[1], W[2]
    tipL = lip[np.argmax(lip[:, 1])]; le2 = r2[np.argmin(r2[:, 1])]; te3 = r3[np.argmax(r3[:, 1])]; bot = lip[-1]
    wins = [((-0.2, P[:, 0].max() + 0.2), (-4.2, 10.2)), ((tipL[0] - 0.35, tipL[0] + 0.45), (tipL[1] - 0.45, tipL[1] + 0.25)),
            ((te3[0] - 0.4, te3[0] + 0.25), (te3[1] - 0.3, te3[1] + 0.2)), ((bot[0] - 0.25, bot[0] + 0.25), (bot[1] - 0.3, bot[1] + 0.3))]
    for ax, (xl, yl), name in zip(axs[row], wins, ("full", "lip TE / ring-2 LE", "ring-3 TE", "mantle end + strip")):
        show(ax, P, Q, xl, yl, f"phi {phi}: {name}")
        for Wp in W:
            ax.plot(Wp[:, 0], Wp[:, 1], "b-", lw=1.2)
plt.savefig("cgrid_sections.png", dpi=60, bbox_inches="tight")
