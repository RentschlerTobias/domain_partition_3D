import json, sys, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, axs = plt.subplots(1, len(sys.argv)-1, figsize=(7*(len(sys.argv)-1), 9))
for ax, f in zip(np.atleast_1d(axs), sys.argv[1:]):
    S = json.load(open(f))
    for L in S["loops"]:
        for d in L:
            P = np.array(d["pts"]); ax.plot(P[:,0], P[:,1], lw=1); ax.text(*P[len(P)//2], str(d["face"]), fontsize=8)
            ax.plot(*P[0], 'k.', ms=4)
    for d in S["slits"]:
        P = np.array(d["pts"]); ax.plot(P[:,0], P[:,1], 'r-', lw=1.5); ax.text(*P[len(P)//2], str(d["face"]), color='r')
    ax.set_aspect('equal'); ax.set_title(f"phi {S['phi']}")
plt.savefig("sections.png", dpi=110, bbox_inches='tight')
