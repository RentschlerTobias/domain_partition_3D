import sys, os, time, multiprocessing as mp; sys.path.insert(0, "/work/canada/blocking_study/sweep_phi/cgrid_grading")
from cgrid2d_g import *
C = CGrid(H_CORE)
C.assemble(*C.geometry(18, Xs=np.zeros((C.m0.N, 2)) + 1)[:3]) if False else None
def run(k):
    t = time.time(); X, LAY, STRIP, Xs = C.geometry(k); X, it = winslow(C.m, X)
    P, Q, B = C.assemble(X, LAY, STRIP); sj = quad_sj(P, Q)
    np.savez(f"test_{os.environ.get('TAG','')}{int(C.phis[k])}.npz", P=P, Q=Q)
    return f"phi {C.phis[k]}: its {it} {time.time()-t:.0f}s, quads {len(Q)}, inverted {(sj<=0).sum()}, min sJ {sj.min():.3f}, p1 {np.percentile(sj,1):.3f}"
if __name__ == "__main__":
    ks = [int(a) for a in sys.argv[1:]]
    # orientation cache from phi 45 first (cheap: geometry without solve is enough for the summed areas)
    X, LAY, STRIP, Xs = C.geometry(18); C.assemble(X, LAY, STRIP)
    with mp.Pool(len(ks)) as p:
        for r in p.map(run, ks): print(r, flush=True)
