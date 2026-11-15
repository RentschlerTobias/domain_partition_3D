import sys, time; sys.path.insert(0, "/work/canada/blocking_study/sweep_phi/cgrid_grading")
from cgrid2d_g import *
C = CGrid(H_CORE); C.UNT = False
free = np.zeros(C.NT, bool); free[C.m.inner] = True
for c in C.chains:
    free[c["lay_ids"][:, 1:-1].ravel()] = True
for f in sys.argv[1:]:
    D = np.load(f); P, Q = D["P"].copy(), D["Q"]
    bad = np.abs(P).max(1) > 1e3
    print(f, "diverged nodes", bad.sum())
    t = time.time(); sj0 = quad_sj(P, Q); P = untangle(P, Q, free); sj = quad_sj(P, Q)
    print(f"  before: inv {(sj0<=0).sum()} min {sj0.min():.3f} -> after: inv {(sj<=0).sum()} min {sj.min():.3f} p1 {np.percentile(sj,1):.3f} ({time.time()-t:.0f}s)", flush=True)
