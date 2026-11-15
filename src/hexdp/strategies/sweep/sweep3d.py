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
usage: python3 sweep3d.py [h2d] [core_a]
"""
import sys, json, time, collections
sys.path.insert(0, "/work/canada/blocking_study/sweep_phi")
from layout2d import *

H = float(sys.argv[1]) if len(sys.argv) > 1 else 0.1
A_CORE = float(sys.argv[2]) if len(sys.argv) > 2 else 0.45
OUT = HERE
t0 = time.time()

secs = load_sections()
phis = np.array([p for p, _ in secs]); K = len(phis) - 1
assert K % 2 == 0 and abs(phis[K // 2] - 45) < 1e-9, phis
K2 = K // 2
secs5 = [s for s in secs if abs(s[0] / 5 - round(s[0] / 5)) < 1e-9]
ml = mean_loops(secs5); ab = abstract(ml)
merges = [tuple(x) for x in json.load(open(f"{OUT}/layout_final.json"))["merges"]]
lay = Layout(ml, ab, merges); NB2 = lay.arrangement()
assert lay.tjunctions == 0
m = Mesh2D(lay, H)
print(f"[2d] {NB2} blocks, {m.N} nodes, {len(m.cells)} quads, phi levels {K + 1}", flush=True)

# ---------------------------------------------------------------- per-phi re-fit (warm start from neighbour)
X = [None] * (K + 1); q2 = {}
order = list(range(K2, -1, -1)) + list(range(K2 + 1, K + 1))
prev = None
for k in order:
    if k == K2 + 1:
        prev = X[K2]
    X[k] = m.solve(secs[k][1], iters=8000, X0=prev); prev = X[k]
    sj = m.quality(X[k]); q2[float(phis[k])] = (int((sj <= 0).sum()), float(sj.min()))
print("[2d] per-phi (inverted, min sJ):", {p: (a, round(b, 3)) for p, (a, b) in q2.items()}, flush=True)
X = np.array(X)                                         # (K+1, N, 2) = (r, z)


def xyz(rz, phi):
    p = np.radians(phi)
    return np.stack([rz[..., 0] * np.cos(p), rz[..., 0] * np.sin(p), rz[..., 1]], -1)


# ---------------------------------------------------------------- 2D blocks -> fine index boxes + node ids
R = lay.R
blocks2d = []
for b, (i0, i1, j0, j1) in enumerate(lay.blocks):
    I0, I1, J0, J1 = m.oi[i0], m.oi[i1], m.oj[j0], m.oj[j1]
    ids = np.zeros((I1 - I0 + 1, J1 - J0 + 1), int)
    for I in range(I0, I1 + 1):
        for J in range(J0, J1 + 1):
            ci = I + 0.5 if I < I1 else I - 0.5; cj = J + 0.5 if J < J1 else J - 0.5
            ids[I - I0, J - J0] = m.keys[m.key(I, J, ci, cj)]
    blocks2d.append(dict(b=b, box=(I0, I1, J0, J1), ids=ids, axis=bool(R[i0] == 0.0)))
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
    return dict(name=f"b{bd['b']}{tag}", b2=bd["b"], X=P, nid=nid, kind="swept")


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
        out.append(dict(name=f"b{bd['b']}{name}", b2=bd["b"], X=P, nid=nid, kind=kind))
    return out


for bd in blocks2d:
    if bd["axis"]:
        blocks3d += butterfly(bd)
    else:
        blocks3d += [swept(bd, 0, K2, "lo"), swept(bd, K2, K, "hi")]
print(f"[3d] {len(blocks3d)} blocks ({sum(b['kind'] == 'swept' for b in blocks3d)} swept, "
      f"{sum(b['kind'] == 'core' for b in blocks3d)} butterfly cores, {sum(b['kind'] == 'half' for b in blocks3d)} butterfly halves)", flush=True)


# ---------------------------------------------------------------- TFI fill
def tfi3d(P):
    n1, n2, n3 = (s - 1 for s in P.shape[:3])
    u = np.linspace(0, 1, n1 + 1)[:, None, None, None]; v = np.linspace(0, 1, n2 + 1)[None, :, None, None]
    w = np.linspace(0, 1, n3 + 1)[None, None, :, None]
    F = lambda a: a
    Pu = (1 - u) * P[:1] + u * P[-1:]
    Pv = (1 - v) * P[:, :1] + v * P[:, -1:]
    Pw = (1 - w) * P[:, :, :1] + w * P[:, :, -1:]
    PuPv = (1 - u) * ((1 - v) * P[:1, :1] + v * P[:1, -1:]) + u * ((1 - v) * P[-1:, :1] + v * P[-1:, -1:])
    PvPw = (1 - v) * ((1 - w) * P[:, :1, :1] + w * P[:, :1, -1:]) + v * ((1 - w) * P[:, -1:, :1] + w * P[:, -1:, -1:])
    PuPw = (1 - u) * ((1 - w) * P[:1, :, :1] + w * P[:1, :, -1:]) + u * ((1 - w) * P[-1:, :, :1] + w * P[-1:, :, -1:])
    PuPvPw = 0
    for a, ia in ((1 - u, 0), (u, -1)):
        for b, ib in ((1 - v, 0), (v, -1)):
            for c, ic in ((1 - w, 0), (w, -1)):
                PuPvPw = PuPvPw + a * b * c * P[ia, ib, ic][None, None, None]
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


res = {}
for mode in ("sweep", "tfi"):
    allsj = []
    for B in blocks3d:
        P = B["X"] if mode == "sweep" else tfi3d(B["X"])
        B[mode] = P
        sj = hex_sj(P); B[f"sj_{mode}"] = sj; allsj.append(sj.ravel())
    allsj = np.concatenate(allsj)
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
res["counts"] = dict(blocks2d=NB2, axis_blocks2d=n_axis, blocks3d=len(blocks3d), phi_levels=K + 1, h2d=H,
                     a_core=A_CORE, nodes=int(len(pts)), hexes=int(len(hexes)), dtoo_volumes=98)
res["per_phi_2d"] = {str(k): v for k, v in q2.items()}

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


kindid = {"swept": 0, "core": 1, "half": 2}
write_vtk_hex(f"{OUT}/sweep_tfi_blocks.vtk", pts, hexes,
              [("block", hblock, "int"), ("block2d", np.array([blocks3d[b]["b2"] for b in hblock]), "int"),
               ("kind", np.array([kindid[blocks3d[b]["kind"]] for b in hblock]), "int"), ("sJ", hsj, "double")])
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
    P3 = np.c_[X[k], np.zeros(m.N)]
    with open(f"{OUT}/layout2d_phi{int(phis[k]):02d}.vtk", "w") as f:
        f.write("# vtk DataFile Version 3.0\nlayout2d\nASCII\nDATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {m.N} double\n" + "".join(f"{a:.9g} {b:.9g} 0\n" for a, b, _ in P3))
        f.write(f"CELLS {len(m.cells)} {5 * len(m.cells)}\n" + "".join("4 " + " ".join(map(str, q[:4])) + "\n" for q in m.cells))
        f.write(f"CELL_TYPES {len(m.cells)}\n" + "9\n" * len(m.cells))
        sj = m.quality(X[k])
        f.write(f"CELL_DATA {len(m.cells)}\nSCALARS block int 1\nLOOKUP_TABLE default\n" + "".join(f"{q[4]}\n" for q in m.cells))
        f.write("SCALARS sJ double 1\nLOOKUP_TABLE default\n" + "".join(f"{v:.6g}\n" for v in sj))
res["time_s"] = time.time() - t0
json.dump(res, open(f"{OUT}/sweep_metrics.json", "w"), indent=1)
print(f"[done] {time.time() - t0:.0f} s", flush=True)
