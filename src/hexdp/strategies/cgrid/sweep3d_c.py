"""2.5D phi sweep of the 2D section layout into 3D hex blocks + TFI fill + evaluation.

- layout: from the 5-deg sections (mean section, abstraction, merges of layout_final.json)
- sections: all phi levels in sections/ (2.5 deg) re-fitted with the SAME fine 2D grid (Winslow, warm start
  from the neighbouring phi)
- 3D blocks: every 2D block x [0,45] and x [45,90] deg (conforming quarter butterfly needs the 45-deg split, see
  report); 2D blocks on the axis: quarter butterfly = core square block + two half blocks (the halves meet
  the swept neighbours on the outer arc).
- every block is then re-filled by 3D TFI from its 6 boundary faces ("tfi"); the swept/butterfly node
  positions before that are kept as "sweep" reference.
- evaluation: inverted hexes, scaled Jacobian, conformity (block faces shared 1:1, quad faces 2/1), boundary
  points written for the Hausdorff check (hausdorff.py, pythonOCC).
usage: python3 sweep3d_c.py [h2d] [core_a] [--cgrid-dir DIR] [--sweep-dir DIR] [--out DIR] [--workers N]
       (cgrid dir holds walls/ and walls.json, default: the working directory; sweep dir holds sections/ and
        layout_final.json, default: its parent; out default: the cgrid dir)
"""
import argparse
import json
import multiprocessing as mp
import os
import sys
import time
import collections

from hexdp.strategies.sweep.layout2d import *  # noqa: F401,F403
from hexdp.strategies.cgrid import cgrid2d as cg

# per-phi fit runs in worker processes (fork): the C-grid model is handed over through this module global
_STATE = {}


def fit(k):
    C, m, phis = _STATE["C"], _STATE["m"], _STATE["phis"]
    Xo, LAY, STRIP, Xs = C.geometry(k)
    Xo, it = cg.winslow(m, Xo)
    P, Q, B = C.assemble(Xo, LAY, STRIP)
    return k, P, it, {kk: v for kk, v in C.blunt_info.items() if kk[0] == phis[k]}


def run(cgrid_dir=".", sweep_dir=None, out=None, h2d=0.1, a_core=0.45, workers=8):
    """C-grid variant of the phi sweep. Returns the metrics dict (also written to ``out``/sweep_metrics.json)."""
    H = h2d
    A_CORE = a_core
    OUT = out or cgrid_dir
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    C = cg.CGrid(H, sweep_dir=sweep_dir, cgrid_dir=cgrid_dir)
    secs = C.secs; lay = C.lay; m = C.m; NB2 = C.nb2
    phis = C.phis; K = len(phis) - 1
    assert K % 2 == 0 and abs(phis[K // 2] - 45) < 1e-9, phis
    K2 = K // 2

    _STATE.update(C=C, m=m, phis=phis)


    k45, P45, _, _ = fit(K2)                              # fixes the block orientation cache before forking
    _, Q2D, XBLK = C.assemble(*C.geometry(K2)[:3])
    with mp.Pool(workers) as pool:
        res_fit = pool.map(fit, range(K + 1))
    X = np.array([r[1] for r in res_fit])                  # (K+1, NT, 2)
    blunt_info = {}
    for r in res_fit:
        blunt_info.update(r[3])
    q2 = {}
    for k in range(K + 1):
        sj = cg.quad_sj(X[k], Q2D); q2[float(phis[k])] = (int((sj <= 0).sum()), float(sj.min()), float(np.percentile(sj, 1)))
    NB2C = NB2 + len(XBLK)
    print(f"[2d] {NB2} outside + {len(XBLK)} C-layer/strip = {NB2C} blocks; nodes {C.NT}, quads {len(Q2D)}; "
          f"layer delta {cg.DELTA}, N_L {cg.N_L}, h1 {C.h1:.3e}, ratio {C.q:.4f}", flush=True)
    print("[2d] per-phi (inverted, min sJ, p1):", {p: (a, round(b, 3), round(c, 3)) for p, (a, b, c) in q2.items()}, flush=True)
    print("[2d] knife-edge blunting (phi, chain): removed wall length", {f"{k[0]}/{k[1]}": round(v, 4) for k, v in blunt_info.items() if v > 0}, flush=True)


    def xyz(rz, phi):
        p = np.radians(phi)
        return np.stack([rz[..., 0] * np.cos(p), rz[..., 0] * np.sin(p), rz[..., 1]], -1)


    R = lay.R
    blocks2d = []
    for b, (i0, i1, j0, j1) in enumerate(lay.blocks):
        I0, I1, J0, J1 = m.oi[i0], m.oi[i1], m.oj[j0], m.oj[j1]
        ids = np.zeros((I1 - I0 + 1, J1 - J0 + 1), int)
        for I in range(I0, I1 + 1):
            for J in range(J0, J1 + 1):
                ci = I + 0.5 if I < I1 else I - 0.5; cj = J + 0.5 if J < J1 else J - 0.5
                ids[I - I0, J - J0] = m.keys[m.key(I, J, ci, cj)]
        blocks2d.append(dict(b=b, ids=ids, axis=bool(R[i0] == 0.0), kind2="outside"))
    for xb in XBLK:
        blocks2d.append(dict(b=xb["bid"], ids=xb["ids"], axis=False, kind2=xb["kind"], chain=xb["chain"]))
    n_axis = sum(b["axis"] for b in blocks2d)
    print(f"[2d] axis-adjacent blocks: {n_axis}", flush=True)

    # ---------------------------------------------------------------- 3D blocks
    blocks3d = []          # dict(name, X (n1+1, n2+1, n3+1, 3) right-handed, nid (same shape) global node keys)


    def swept(bd, k0, k1, tag):
        ids = bd["ids"]
        rz = X[k0:k1 + 1][:, ids]                            # (nk, ni, nj, 2)
        P = xyz(rz, phis[k0:k1 + 1][:, None, None])          # (nk, ni, nj, 3)
        P = P.transpose(1, 0, 2, 3)                          # (i, k, j): r, phi, z -> right-handed
        nid = np.empty(P.shape[:3], object)
        for a in range(P.shape[0]):
            for c in range(P.shape[1]):
                for d in range(P.shape[2]):
                    nid[a, c, d] = ("s", int(ids[a, d]), k0 + c)
        return dict(name=f"b{bd['b']}{tag}", b2=bd["b"], X=P, nid=nid, kind="swept", kind2=bd["kind2"])


    def interp_section(ids_line, s, phi):
        """(r,z) on the radial grid line ids_line (axis -> outer) at fraction s, phi (deg), linear in both"""
        n = len(ids_line) - 1
        f = np.clip(s, 0, 1) * n; i = np.minimum(np.floor(f).astype(int), n - 1); w = f - i
        kf = np.interp(phi, phis, np.arange(K + 1)); k = np.minimum(np.floor(kf).astype(int), K - 1); v = kf - k
        def at(kk):
            return (1 - w)[..., None] * X[kk, ids_line[i]] + w[..., None] * X[kk, ids_line[i + 1]]
        return (1 - v)[..., None] * at(k) + v[..., None] * at(k + 1)


    def butterfly(bd):
        ids = bd["ids"]; nI = ids.shape[0] - 1; nJ = ids.shape[1] - 1
        nR = max(2, int(round(nI * (1 - A_CORE))))
        a = A_CORE
        out = []
        t = np.linspace(0, 1, K2 + 1); l = np.linspace(0, 1, nR + 1)
        # core: (p, q) in [0,a]^2
        xi, eta = np.meshgrid(a * t, a * t, indexing="ij")
        core = (xi, eta)
        # half A (phi 0..45): l radial (inner = square side xi=a), k along phi
        ph = phis[:K2 + 1]
        inA = np.stack([np.full(K2 + 1, a), a * t], -1); outA = np.stack([np.cos(np.radians(ph)), np.sin(np.radians(ph))], -1)
        hA = (1 - l)[:, None, None] * inA[None] + l[:, None, None] * outA[None]          # (l, k, 2)
        phB = phis[K2:]
        inB = np.stack([a * (1 - t), np.full(K2 + 1, a)], -1); outB = np.stack([np.cos(np.radians(phB)), np.sin(np.radians(phB))], -1)
        hB = (1 - l)[:, None, None] * inB[None] + l[:, None, None] * outB[None]
        for name, (U, V), kind in (("core", core, "core"), ("hA", (hA[..., 0], hA[..., 1]), "half"), ("hB", (hB[..., 0], hB[..., 1]), "half")):
            s = np.hypot(U, V); ph_ = np.degrees(np.arctan2(V, U))
            P = np.zeros(U.shape + (nJ + 1, 3)); nid = np.empty(U.shape + (nJ + 1,), object)
            for J in range(nJ + 1):
                rz = interp_section(ids[:, J], s, ph_)
                P[:, :, J] = xyz(rz, ph_)
                for a1 in range(U.shape[0]):
                    for a2 in range(U.shape[1]):
                        if kind == "half" and a1 == nR:                 # outer arc = swept neighbour nodes
                            kk = a2 if name == "hA" else K2 + a2
                            nid[a1, a2, J] = ("s", int(ids[nI, J]), kk)
                        else:
                            nid[a1, a2, J] = ("c", tuple(np.round(P[a1, a2, J], 9)))
            out.append(dict(name=f"b{bd['b']}{name}", b2=bd["b"], X=P, nid=nid, kind=kind, kind2=kind))
        return out


    for bd in blocks2d:
        if bd["axis"]:
            blocks3d += butterfly(bd)
        else:
            blocks3d += [swept(bd, 0, K2, "lo"), swept(bd, K2, K, "hi")]
    print(f"[3d] {len(blocks3d)} blocks ({sum(b['kind'] == 'swept' for b in blocks3d)} swept, "
          f"{sum(b['kind'] == 'core' for b in blocks3d)} butterfly cores, {sum(b['kind'] == 'half' for b in blocks3d)} butterfly halves)", flush=True)


    # ---------------------------------------------------------------- TFI fill
    def _frac(E):
        """normalised arc length along the first axis of an edge array (n, 3)"""
        d = np.r_[0, np.cumsum(np.linalg.norm(np.diff(E, axis=0), axis=-1))]
        return d / d[-1] if d[-1] > 0 else np.linspace(0, 1, len(E))


    def tfi3d(P):
        """3D TFI (Boolean sum of the three linear projectors) with arc-length based blending parameters: u(i,j,k)
        is the bilinear blend (in the index fractions of j,k) of the arc-length fractions along the 4 block edges in i"""
        n1, n2, n3 = (s_ - 1 for s_ in P.shape[:3])
        I = np.linspace(0, 1, n1 + 1)[:, None, None]; J = np.linspace(0, 1, n2 + 1)[None, :, None]; Kk = np.linspace(0, 1, n3 + 1)[None, None, :]
        e = lambda E: _frac(E)
        u = ((1 - J) * (1 - Kk) * e(P[:, 0, 0])[:, None, None] + J * (1 - Kk) * e(P[:, -1, 0])[:, None, None]
             + (1 - J) * Kk * e(P[:, 0, -1])[:, None, None] + J * Kk * e(P[:, -1, -1])[:, None, None])
        v = ((1 - I) * (1 - Kk) * e(P[0, :, 0])[None, :, None] + I * (1 - Kk) * e(P[-1, :, 0])[None, :, None]
             + (1 - I) * Kk * e(P[0, :, -1])[None, :, None] + I * Kk * e(P[-1, :, -1])[None, :, None])
        w = ((1 - I) * (1 - J) * e(P[0, 0, :])[None, None, :] + I * (1 - J) * e(P[-1, 0, :])[None, None, :]
             + (1 - I) * J * e(P[0, -1, :])[None, None, :] + I * J * e(P[-1, -1, :])[None, None, :])
        u, v, w = u[..., None], v[..., None], w[..., None]
        Pu = (1 - u) * P[:1] + u * P[-1:]
        Pv = (1 - v) * P[:, :1] + v * P[:, -1:]
        Pw = (1 - w) * P[:, :, :1] + w * P[:, :, -1:]
        PuPv = (1 - u) * ((1 - v) * P[:1, :1] + v * P[:1, -1:]) + u * ((1 - v) * P[-1:, :1] + v * P[-1:, -1:])
        PvPw = (1 - v) * ((1 - w) * P[:, :1, :1] + w * P[:, :1, -1:]) + v * ((1 - w) * P[:, -1:, :1] + w * P[:, -1:, -1:])
        PuPw = (1 - u) * ((1 - w) * P[:1, :, :1] + w * P[:1, :, -1:]) + u * ((1 - w) * P[-1:, :, :1] + w * P[-1:, :, -1:])
        PuPvPw = 0
        for a_, ia in ((1 - u, 0), (u, -1)):
            for b_, ib in ((1 - v, 0), (v, -1)):
                for c_, ic in ((1 - w, 0), (w, -1)):
                    PuPvPw = PuPvPw + a_ * b_ * c_ * P[ia, ib, ic][None, None, None]
        T = Pu + Pv + Pw - PuPv - PvPw - PuPw + PuPvPw
        Q = P.copy(); Q[1:-1, 1:-1, 1:-1] = T[1:-1, 1:-1, 1:-1]
        return Q


    def hex_sj(P):
        """min corner scaled Jacobian per hex of a structured block (right-handed index order)"""
        c = [P[:-1, :-1, :-1], P[1:, :-1, :-1], P[:-1, 1:, :-1], P[1:, 1:, :-1],
             P[:-1, :-1, 1:], P[1:, :-1, 1:], P[:-1, 1:, 1:], P[1:, 1:, 1:]]
        # corner (dx,dy,dz) -> neighbours along +-each axis
        idx = lambda x, y, z: c[x + 2 * y + 4 * z]
        sj = None
        for x in (0, 1):
            for y in (0, 1):
                for z in (0, 1):
                    p = idx(x, y, z)
                    e1 = (idx(1 - x, y, z) - p) * (1 if x == 0 else -1)
                    e2 = (idx(x, 1 - y, z) - p) * (1 if y == 0 else -1)
                    e3 = (idx(x, y, 1 - z) - p) * (1 if z == 0 else -1)
                    d = np.einsum("...i,...i->...", np.cross(e1, e2), e3)
                    s = d / (np.linalg.norm(e1, axis=-1) * np.linalg.norm(e2, axis=-1) * np.linalg.norm(e3, axis=-1) + 1e-300)
                    sj = s if sj is None else np.minimum(sj, s)
        return sj


    def untangle_block(P, thr=0.05, rounds=8):
        """move block-interior nodes of hexes with sJ < thr to maximise the min sJ of the 8 surrounding hexes"""
        from scipy.optimize import minimize
        n1, n2, n3 = (x - 1 for x in P.shape[:3]); moved = 0
        for _ in range(rounds):
            sj = hex_sj(P); bad = np.argwhere(sj < thr)
            if not len(bad):
                break
            nodes = {(i + a, j + b, k + c) for i, j, k in bad for a in (0, 1) for b in (0, 1) for c in (0, 1)}
            nodes = [x for x in nodes if 0 < x[0] < n1 and 0 < x[1] < n2 and 0 < x[2] < n3]
            if not nodes:
                break
            for (i, j, k) in nodes:
                sl = (slice(max(i - 1, 0), i + 2), slice(max(j - 1, 0), j + 2), slice(max(k - 1, 0), k + 2))
                x0 = P[i, j, k].copy(); h = np.linalg.norm(P[min(i + 1, n1), j, k] - P[i, j, k]) + 1e-12
                def f(x):
                    P[i, j, k] = x; return -hex_sj(P[sl]).min()
                f0 = f(x0)
                r = minimize(f, x0, method="Nelder-Mead", options=dict(xatol=1e-5 * h, fatol=1e-7, maxiter=400,
                             initial_simplex=[x0, x0 + [0.3 * h, 0, 0], x0 + [0, 0.3 * h, 0], x0 + [0, 0, 0.3 * h]]))
                P[i, j, k] = r.x if r.fun <= f0 else x0; moved += 1
        return P, moved


    res = {}
    for mode in ("sweep", "tfi", "tfi_untangled"):
        allsj = []
        for B in blocks3d:
            if mode == "sweep":
                P = B["X"]
            elif mode == "tfi":
                P = tfi3d(B["X"])
            else:
                P = B["tfi"].copy()
                if (B["sj_tfi"] < 0.05).any():
                    P, nm = untangle_block(P); B["untangled_nodes"] = nm
            B[mode] = P
            sj = hex_sj(P); B[f"sj_{mode}"] = sj; allsj.append(sj.ravel())
        allsj = np.concatenate(allsj)
        for B in blocks3d:
            bad = np.argwhere(B[f"sj_{mode}"] <= 0)
            if len(bad):
                Pm = B[mode]
                print(f"   [{mode}] inverted in {B['name']} shape {B[f'sj_{mode}'].shape}: idx {bad[:8].tolist()} at",
                      [np.round(Pm[tuple(i)], 3).tolist() for i in bad[:3]], flush=True)
        res[mode] = dict(hexes=int(allsj.size), inverted=int((allsj <= 0).sum()), min_sJ=float(allsj.min()),
                         p1_sJ=float(np.percentile(allsj, 1)), median_sJ=float(np.median(allsj)))
        worst = sorted(((float(B[f"sj_{mode}"].min()), B["name"]) for B in blocks3d))[:6]
        res[mode]["worst_blocks"] = worst
        print(f"[3d {mode}] {res[mode]}", flush=True)
    # TFI fill of the 2D sections alone (block-wise 2D TFI of the phi=0/45/90 faces is part of the 3D faces)

    # ---------------------------------------------------------------- global mesh (node merge by keys)
    def gkey(k):
        if k[0] == "s":
            return k
        return k


    keymap = {}; pts = []; hexes = []; hblock = []; hsj = []
    for B in blocks3d:
        B["tfi_raw"] = B["tfi"]; B["tfi"] = B["tfi_untangled"]; B["sj_tfi_raw"] = B["sj_tfi"]; B["sj_tfi"] = B["sj_tfi_untangled"]
    res["untangled_blocks"] = {B["name"]: B.get("untangled_nodes", 0) for B in blocks3d if B.get("untangled_nodes")}
    print("[3d] untangled interior nodes per block:", res["untangled_blocks"], flush=True)
    for bi, B in enumerate(blocks3d):
        P = B["tfi"]; nid = B["nid"]
        gid = np.zeros(nid.shape, int)
        for ix in np.ndindex(nid.shape):
            kk = nid[ix]
            if kk not in keymap:
                keymap[kk] = len(pts); pts.append(P[ix])
            gid[ix] = keymap[kk]
        B["gid"] = gid
        c = [gid[:-1, :-1, :-1], gid[1:, :-1, :-1], gid[1:, 1:, :-1], gid[:-1, 1:, :-1],
             gid[:-1, :-1, 1:], gid[1:, :-1, 1:], gid[1:, 1:, 1:], gid[:-1, 1:, 1:]]
        H8 = np.stack([x.ravel() for x in c], 1); hexes.append(H8); hblock += [bi] * len(H8); hsj.append(B["sj_tfi"].ravel())
    pts = np.array(pts); hexes = np.vstack(hexes); hblock = np.array(hblock); hsj = np.concatenate(hsj)
    # position consistency of shared nodes (TFI keeps boundary faces, so shared nodes must coincide)
    maxdev = 0.0
    for B in blocks3d:
        maxdev = max(maxdev, float(np.abs(pts[B["gid"]] - B["tfi"]).max()))
    print(f"[mesh] {len(pts)} nodes, {len(hexes)} hexes, max shared-node mismatch {maxdev:.2e}", flush=True)

    # quad faces: interior (2 hexes) / boundary (1)
    FACES = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
    fc = collections.Counter()
    for f in FACES:
        for row in np.sort(hexes[:, f], axis=1):
            fc[tuple(row)] += 1
    cnt = collections.Counter(fc.values())
    bfaces = [f for f, n in fc.items() if n == 1]
    # boundary edges: each must be shared by an even number (>=2) of boundary quads (watertight incl. baffles)
    ec = collections.Counter()
    quad_of = {}
    for f in FACES:
        for row in hexes[:, f]:
            key = tuple(sorted(row))
            if fc[key] == 1:
                quad_of[key] = row
    for key, row in quad_of.items():
        for a, b in zip(row, np.roll(row, -1)):
            ec[(min(a, b), max(a, b))] += 1
    odd = sum(1 for v in ec.values() if v % 2)
    print(f"[mesh] quad faces shared by n hexes: {dict(cnt)}; boundary quads {len(bfaces)}; boundary edges with odd "
          f"face count: {odd}", flush=True)

    # block-level conformity: every block face is either all-boundary or matches exactly one face of one other block
    def bfaces_of(B):
        g = B["gid"]
        return {"i0": g[0], "i1": g[-1], "j0": g[:, 0], "j1": g[:, -1], "k0": g[:, :, 0], "k1": g[:, :, -1]}
    face_index = collections.defaultdict(list)
    for bi, B in enumerate(blocks3d):
        for fn, F in bfaces_of(B).items():
            face_index[frozenset(F.ravel().tolist())].append((bi, fn))
    nonconf = 0; shared = 0; bnd_faces = 0
    for bi, B in enumerate(blocks3d):
        for fn, F in bfaces_of(B).items():
            S = frozenset(F.ravel().tolist())
            if len(face_index[S]) == 2:
                shared += 1; continue
            # boundary face: all its quads must be boundary quads
            G = F
            quads = [tuple(sorted((G[a, b], G[a + 1, b], G[a + 1, b + 1], G[a, b + 1]))) for a in range(G.shape[0] - 1) for b in range(G.shape[1] - 1)]
            if all(fc[q] == 1 for q in quads):
                bnd_faces += 1
            else:
                nonconf += 1
    print(f"[blocks] block faces: {shared // 2} shared pairs, {bnd_faces} on the boundary, {nonconf} non-conforming", flush=True)
    res["conformity"] = dict(quad_faces_by_count={int(k): v for k, v in cnt.items()}, boundary_quads=len(bfaces),
                             odd_boundary_edges=odd, shared_block_face_pairs=shared // 2, boundary_block_faces=bnd_faces,
                             nonconforming_block_faces=nonconf, max_shared_node_mismatch=maxdev)
    res["counts"] = dict(blocks2d_outside=NB2, blocks2d_extra=len(XBLK), blocks2d=NB2C, layer_blocks2d=sum(b["kind"] == "layer" for b in XBLK),
                         strip_blocks2d=sum(b["kind"] == "strip" for b in XBLK), delta=cg.DELTA, n_layer=cg.N_L, h1=C.h1, ratio=C.q, axis_blocks2d=n_axis, blocks3d=len(blocks3d), phi_levels=K + 1, h2d=H,
                         a_core=A_CORE, nodes=int(len(pts)), hexes=int(len(hexes)), dtoo_volumes=98)
    # ---------------------------------------------------------------- outputs
    def write_vtk_hex(fn, P, Hx, cdata):
        with open(fn, "wb") as f:
            f.write(b"# vtk DataFile Version 3.0\nsweep_phi\nBINARY\nDATASET UNSTRUCTURED_GRID\n")
            f.write(f"POINTS {len(P)} double\n".encode()); f.write(P.astype(">f8").tobytes()); f.write(b"\n")
            f.write(f"CELLS {len(Hx)} {9 * len(Hx)}\n".encode())
            f.write(np.c_[np.full(len(Hx), 8), Hx].astype(">i4").tobytes()); f.write(b"\n")
            f.write(f"CELL_TYPES {len(Hx)}\n".encode()); f.write(np.full(len(Hx), 12, ">i4").tobytes()); f.write(b"\n")
            f.write(f"CELL_DATA {len(Hx)}\n".encode())
            for name, arr, typ in cdata:
                f.write(f"SCALARS {name} {typ} 1\nLOOKUP_TABLE default\n".encode())
                f.write(arr.astype(">i4" if typ == "int" else ">f8").tobytes()); f.write(b"\n")


    def write_vtk_lines(fn, polylines, cdata):
        P = np.vstack(polylines); off = np.cumsum([0] + [len(p) for p in polylines])
        with open(fn, "w") as f:
            f.write("# vtk DataFile Version 3.0\nblock edges\nASCII\nDATASET POLYDATA\n")
            f.write(f"POINTS {len(P)} double\n" + "".join(f"{x:.9g} {y:.9g} {z:.9g}\n" for x, y, z in P))
            f.write(f"LINES {len(polylines)} {len(P) + len(polylines)}\n")
            for a, b in zip(off[:-1], off[1:]):
                f.write(f"{b - a} " + " ".join(map(str, range(a, b))) + "\n")
            f.write(f"CELL_DATA {len(polylines)}\n")
            for name, arr in cdata:
                f.write(f"SCALARS {name} int 1\nLOOKUP_TABLE default\n" + "".join(f"{int(v)}\n" for v in arr))


    kindid = {"outside": 0, "core": 1, "half": 2, "layer": 3, "strip": 4}
    write_vtk_hex(f"{OUT}/sweep_tfi_blocks.vtk", pts, hexes,
                  [("block", hblock, "int"), ("block2d", np.array([blocks3d[b]["b2"] for b in hblock]), "int"),
                   ("kind", np.array([kindid[blocks3d[b]["kind2"]] for b in hblock]), "int"), ("sJ", hsj, "double")])
    # pre-TFI (sweep) positions for comparison
    pts_s = pts.copy()
    for B in blocks3d:
        pts_s[B["gid"]] = B["sweep"]
    write_vtk_hex(f"{OUT}/sweep_noTFI_blocks.vtk", pts_s, hexes,
                  [("block", hblock, "int"), ("sJ", np.concatenate([B["sj_sweep"].ravel() for B in blocks3d]), "double")])
    # block edges
    lines, lid = [], []
    for bi, B in enumerate(blocks3d):
        P = B["tfi"]
        for ax in range(3):
            sl = [slice(None)] * 3
            for e1 in (0, -1):
                for e2 in (0, -1):
                    o = [d for d in range(3) if d != ax]
                    idx = [None] * 3; idx[ax] = slice(None); idx[o[0]] = e1; idx[o[1]] = e2
                    lines.append(P[tuple(idx)]); lid.append(bi)
    write_vtk_lines(f"{OUT}/sweep_block_edges.vtk", lines, [("block", lid)])
    # boundary points for the Hausdorff check: nodes + quad centres of boundary quads
    bq = np.array([quad_of[k] for k in quad_of])
    np.savez(f"{OUT}/boundary_pts.npz", nodes=pts[np.unique(bq)], centres=pts[bq].mean(1), quads=bq, pts=pts)
    # 2D layout VTK (phi = 0, 45, 90) with block id
    for k in (0, K2, K):
        P2 = X[k]; NT = len(P2)
        with open(f"{OUT}/layout2d_c_phi{int(phis[k]):02d}.vtk", "w") as f:
            f.write("# vtk DataFile Version 3.0\nlayout2d\nASCII\nDATASET UNSTRUCTURED_GRID\n")
            f.write(f"POINTS {NT} double\n" + "".join(f"{a_:.9g} {b_:.9g} 0\n" for a_, b_ in P2))
            f.write(f"CELLS {len(Q2D)} {5 * len(Q2D)}\n" + "".join("4 " + " ".join(map(str, q[:4])) + "\n" for q in Q2D))
            f.write(f"CELL_TYPES {len(Q2D)}\n" + "9\n" * len(Q2D))
            sj = cg.quad_sj(P2, Q2D)
            f.write(f"CELL_DATA {len(Q2D)}\nSCALARS block int 1\nLOOKUP_TABLE default\n" + "".join(f"{q[4]}\n" for q in Q2D))
            f.write("SCALARS sJ double 1\nLOOKUP_TABLE default\n" + "".join(f"{v:.6g}\n" for v in sj))
        np.savez(f"{OUT}/sec2d_phi{int(phis[k]):02d}.npz", P=P2, Q=Q2D)
    # first cell height on the walls (after TFI): layer blocks, distance wall node -> first layer node
    fc = []
    for B in blocks3d:
        if B["kind2"] != "layer":
            continue
        P = B["tfi"]                                      # (a, phi, b) with b = layer direction (0 wall) or a
        bd = next(x for x in blocks2d if x["b"] == B["b2"])
        ids = bd["ids"]
        if np.isin(ids[:, 0], C.chains[bd["chain"]]["Lid"][:, 0]).all():
            d = np.linalg.norm(P[:, :, 1] - P[:, :, 0], axis=-1)
        elif np.isin(ids[:, -1], C.chains[bd["chain"]]["Lid"][:, 0]).all():
            d = np.linalg.norm(P[:, :, -2] - P[:, :, -1], axis=-1)
        else:
            continue
        fc.append(d.ravel())
    fc = np.concatenate(fc)
    res["first_cell"] = dict(min=float(fc.min()), median=float(np.median(fc)), max=float(fc.max()), target_h1=C.h1)
    print("[wall] first cell height after TFI:", res["first_cell"], flush=True)
    res["per_phi_2d"] = {str(k): v for k, v in q2.items()}
    res["blunting"] = {f"{k[0]}/{k[1]}": v for k, v in blunt_info.items() if v > 0}
    res["time_s"] = time.time() - t0
    json.dump(res, open(f"{OUT}/sweep_metrics.json", "w"), indent=1)
    print(f"[done] {time.time() - t0:.0f} s", flush=True)
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("h2d", type=float, nargs="?", default=0.1)
    ap.add_argument("core_a", type=float, nargs="?", default=0.45)
    ap.add_argument("--cgrid-dir", default=os.getcwd())
    ap.add_argument("--sweep-dir", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args(argv)
    run(a.cgrid_dir, a.sweep_dir, a.out, a.h2d, a.core_a, a.workers)


if __name__ == "__main__":
    main()
