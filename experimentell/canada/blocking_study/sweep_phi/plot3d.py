import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
L = open("sweep_block_edges.vtk").read().split("\n")
i = next(k for k,l in enumerate(L) if l.startswith("POINTS")); n = int(L[i].split()[1])
P = np.array([[float(x) for x in l.split()] for l in L[i+1:i+1+n]])
j = next(k for k,l in enumerate(L) if l.startswith("LINES")); m = int(L[j].split()[1])
lines = [[int(x) for x in l.split()[1:]] for l in L[j+1:j+1+m]]
fig = plt.figure(figsize=(18, 9))
ax = fig.add_subplot(1, 2, 1, projection="3d")
for ln in lines:
    Q = P[ln]; ax.plot(Q[:,0], Q[:,1], Q[:,2], "k-", lw=0.5)
ax.set_box_aspect((11, 11, 14)); ax.view_init(20, 35); ax.set_title("sweep: 3D block edges (72 blocks)")
ax2 = fig.add_subplot(1, 2, 2)
z0 = 8.0
# cross-section of block edges with z = z0 (edges along phi / in the cross-section near z0)
for ln in lines:
    Q = P[ln]
    if np.ptp(Q[:,2]) < 1e-6 and abs(Q[0,2] - z0) < 1.5 and False: pass
# draw edges whose z range contains z0 approximately: intersection points
pts = []
for ln in lines:
    Q = P[ln]
    for a, b in zip(Q[:-1], Q[1:]):
        if (a[2]-z0)*(b[2]-z0) < 0:
            t = (z0-a[2])/(b[2]-a[2]); pts.append(a + t*(b-a))
pts = np.array(pts)
# also draw phi-direction edges lying in planes near z0 (block edges at constant section line)
for ln in lines:
    Q = P[ln]
    if np.ptp(Q[:,2]) < 0.02 * (1 + np.ptp(Q[:,0]) + np.ptp(Q[:,1])) and abs(Q[:,2].mean() - 6.3) < 0.05:
        ax2.plot(Q[:,0], Q[:,1], "b-", lw=0.6)
ax2.plot(pts[:,0], pts[:,1], "r.", ms=3); ax2.set_aspect("equal"); ax2.set_title("blue: block edges in the plane z~6.3 (ring-3 TE level); red: block edges crossing z=8")
plt.savefig("blocks3d.png", dpi=80, bbox_inches="tight")
