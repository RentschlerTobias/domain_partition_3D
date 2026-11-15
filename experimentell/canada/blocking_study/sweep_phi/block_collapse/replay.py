"""Replay a collapse path (records.json / best_*.json 'path') on the input structure at another resolution,
evaluate, check conformity, write ParaView files and the boundary points for the exact-CAD Hausdorff check.
usage: python replay.py --vtk in.vtk --cad cad.npz --path best.json --h 0.1 --out DIR --tag best
       (--path none = baseline re-fit)"""
import sys, os, json, argparse, collections, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import block_collapse as bc
from bstruct import *


def conformity(st):
    """shared block faces: max node mismatch; fine mesh: quads shared by 1/2 hexes, odd boundary edges"""
    ft = {}; worst = 0.0
    for X, c8 in zip(st["X"], st["c8s"]):
        C = corner_cube(c8)
        for a, s in FACES:
            T = np.take(C, s, axis=a); k = frozenset(T.ravel().tolist()); G = face_slice(X, a, s)
            if k in ft:
                G0, C0 = ft[k]; r = orient(G0, C0, T)
                worst = max(worst, float(np.abs(r[0] - G).max()) if r[0].shape == G.shape else np.inf)
            else:
                ft[k] = (G, T)
    # global node merge by coordinates
    pts = np.vstack([X.reshape(-1, 3) for X in st["X"]])
    key = np.round(pts / 1e-7).astype(np.int64)
    _, gid, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    inv = inv.ravel(); P = pts[gid]
    H = []; off = 0
    for X in st["X"]:
        n1, n2, n3 = (s - 1 for s in X.shape[:3])
        idx = inv[off:off + X.shape[0] * X.shape[1] * X.shape[2]].reshape(X.shape[:3]); off += idx.size
        c = [idx[:-1, :-1, :-1], idx[1:, :-1, :-1], idx[1:, 1:, :-1], idx[:-1, 1:, :-1],
             idx[:-1, :-1, 1:], idx[1:, :-1, 1:], idx[1:, 1:, 1:], idx[:-1, 1:, 1:]]
        H.append(np.stack([x.ravel() for x in c], 1))
    H = np.vstack(H)
    HF = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
    Fq = np.vstack([H[:, f] for f in HF]); Ks = np.sort(Fq, 1)
    u, first, cnt = np.unique(Ks, axis=0, return_index=True, return_counts=True)
    bq = Fq[first[cnt == 1]]
    E = np.sort(np.vstack([bq[:, [0, 1]], bq[:, [1, 2]], bq[:, [2, 3]], bq[:, [3, 0]]]), 1)
    _, ec = np.unique(E, axis=0, return_counts=True)
    return dict(max_shared_face_mismatch=worst, quad_count_hist={int(k): int(v) for k, v in zip(*np.unique(cnt, return_counts=True))},
                boundary_quads=int(len(bq)), odd_boundary_edges=int((ec % 2).sum()), nodes=int(len(P)), hexes=int(len(H))), P, bq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vtk", required=True); ap.add_argument("--cad", required=True); ap.add_argument("--path", required=True)
    ap.add_argument("--h", type=float, default=0.1); ap.add_argument("--out", required=True); ap.add_argument("--tag", default="best")
    ap.add_argument("--tol", type=float, default=0.01); ap.add_argument("--nrev", type=int, default=60000)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); t0 = time.time()
    cad = CAD(a.cad); P, grids = read_blocks_vtk(a.vtk)
    st0, surf, _ = bc.initial_state(P, grids, cad, a.tol); ref_cov = bc.coverage(st0, cad)
    allQ = [(S.A, S.B, S.C) for S in surf.values()]
    surf_all = TriSurface(*(np.vstack([x[i] for x in allQ]) for i in range(3)))
    rev = cad.fp[np.random.default_rng(0).choice(len(cad.fp), min(a.nrev, len(cad.fp)), replace=False)]
    path = [] if a.path == "none" else json.load(open(a.path))["path"]
    Xb, ib = bc.realise(st0, dict(c8s=st0["c8s"], ncls=st0["ncls"], ecls=st0["ecls"], fcls=st0["fcls"]), None, cad, surf, a.h)
    st = dict(st0); st["X"] = Xb; ev = bc.evaluate(Xb, st, ib, cad, surf_all, rev)
    print(f"[{a.tag}] step 0: {ev}", flush=True)
    for k, step in enumerate(path):
        want = set(step["sheet_origin"])
        sheets = sheets_of(st["c8s"])
        s = [s for s in sheets if {st["origin"][b] for b, _ in s["members"]} == want]
        assert len(s) == 1, f"step {k}: sheet not found"
        child, geo = bc.collapse_topo(st, s[0], cad, ref_cov, step["variant"])
        Xn, inf = bc.realise(st, child, geo, cad, surf, a.h)
        ev = bc.evaluate(Xn, child, inf, cad, surf_all, rev)
        child["X"] = Xn; child["path"] = st.get("path", []) + [step]; st = child
        print(f"[{a.tag}] step {k + 1}: {ev}", flush=True)
    conf, Pm, bq = conformity(st)
    print(f"[{a.tag}] conformity {conf}", flush=True)
    st.setdefault("path", [])
    bc.write_state(f"{a.out}/{a.tag}_h{a.h}", st, dict(ev, conformity=conf), cad)
    np.savez(f"{a.out}/{a.tag}_h{a.h}_boundary_pts.npz", nodes=Pm[np.unique(bq)], centres=Pm[bq].mean(1), quads=bq, pts=Pm)
    print(f"[{a.tag}] done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
