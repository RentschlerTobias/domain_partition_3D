import sys, time; sys.path.insert(0, "/work/canada/blocking_study/sweep_phi/cgrid_grading")
from cgrid2d_g import *
t = time.time(); C = CGrid(H_CORE)
print("setup", round(time.time() - t, 1), "s; nodes", C.NT, "nr", C.m.nr, "nz", C.m.nz, flush=True)
C.assemble(*C.geometry(18)[:3])
for k in [int(a) for a in sys.argv[1:]]:
    t = time.time(); X, LAY, STRIP, Xs = C.geometry(k); X, it = winslow(C.m, X)
    P, Q, B = C.assemble(X, LAY, STRIP); sj = quad_sj(P, Q)
    print(f"phi {C.phis[k]}: its {it} {time.time()-t:.0f}s, quads {len(Q)}, inverted {(sj<=0).sum()}, min sJ {sj.min():.3f}, p1 {np.percentile(sj,1):.3f}", flush=True)
    np.savez(f"test_{int(C.phis[k])}.npz", P=P, Q=Q)
