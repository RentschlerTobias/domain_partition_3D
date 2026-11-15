import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, axs = plt.subplots(2, 4, figsize=(24, 12))
for row, phi in enumerate((0.0, 90.0)):
    W = json.load(open(f"walls/wall_{phi:05.1f}.json")); S = json.load(open(f"../sections/sec_{phi:05.1f}.json"))
    for ax, (xl, yl) in zip(axs[row], (((1.7, 3.0), (-1.9, 2.6)), ((1.7, 2.7), (1.9, 2.5)), ((1.8, 3.2), (3.9, 4.6)), ((1.9, 6.5), (5.9, 6.5)))):
        for w in W["walls"]:
            P = np.array(w["pts"]); ax.plot(P[:, 0], P[:, 1], "b-", lw=1)
        for d in S["slits"]:
            P = np.array(d["pts"]); ax.plot(P[:, 0], P[:, 1], "r-", lw=1)
        for d in S["loops"][0]:
            P = np.array(d["pts"]); ax.plot(P[:, 0], P[:, 1], "k-", lw=1)
        ax.set_xlim(*xl); ax.set_ylim(*yl); ax.set_aspect("equal"); ax.set_title(f"phi {phi}")
plt.savefig("walls_zoom.png", dpi=70, bbox_inches="tight")
