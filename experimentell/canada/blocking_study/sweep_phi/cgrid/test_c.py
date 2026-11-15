import sys, time; sys.path.insert(0, "/work/canada/blocking_study/sweep_phi/cgrid")
from cgrid2d import *
t = time.time(); C = CGrid(0.1); print("setup", round(time.time() - t, 1), "s; h1", C.h1, "q", C.q)
for c in C.chains:
    print("chain slit", c["slit"], "wall", c["wall"], "closed", c["closed"], "edges", c["edges"], "nodes", len(c["ids"]),
          "start", None if c["start"] is None else c["start"]["kind"], "end", None if c["end"] is None else c["end"]["kind"])
for k in (0, 7, 8, 18, 30, 36):
    X, LAY, STRIP, Xs = C.geometry(k)
    X, it = winslow(C.m, X)
    P, Q, B = C.assemble(X, LAY, STRIP)
    sj = quad_sj(P, Q)
    print(f"phi {C.phis[k]}: its {it}, blocks {C.nb2}+{len(B)}, quads {len(Q)}, inverted {(sj<=0).sum()}, min sJ {sj.min():.3f}, p1 {np.percentile(sj,1):.3f}")
    np.savez(f"test_{int(C.phis[k])}.npz", P=P, Q=Q)
