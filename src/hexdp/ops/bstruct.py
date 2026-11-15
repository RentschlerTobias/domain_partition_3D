"""Generic multi-block hex structure: I/O, structured-grid helpers, CAD classification, sheets.

A block structure here = conforming coarse hex mesh (one coarse hex per block, block corners = coarse nodes)
plus, per block, a structured fine grid X (n1+1, n2+1, n3+1, 3) whose corners are the coarse nodes in VTK
hexahedron order:  c0=(0,0,0) c1=(1,0,0) c2=(1,1,0) c3=(0,1,0) c4=(0,0,1) c5=(1,0,1) c6=(1,1,1) c7=(0,1,1).
Nothing in here knows the canada geometry: input = a hex VTK with a per-cell block id + a CAD export
(cad_export.py).
"""
import itertools, collections, json
import numpy as np
from scipy.spatial import cKDTree

# corner (i,j,k) offsets in VTK order and the parallel classes (axis 0 = i, 1 = j, 2 = k)
OFF = np.array([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)])
PARALLEL = ([(0, 1), (3, 2), (7, 6), (4, 5)], [(0, 3), (1, 2), (5, 6), (4, 7)], [(0, 4), (1, 5), (2, 6), (3, 7)])
HEX_EDGES = ((0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7))
# block faces as (axis, side) and the slice of the 2x2x2 corner array
FACES = [(a, s) for a in range(3) for s in (0, 1)]


def corner_cube(c8):
    """8 VTK-ordered ids -> 2x2x2 array C[i,j,k]"""
    C = np.empty((2, 2, 2), dtype=np.int64)
    for n, (i, j, k) in enumerate(OFF):
        C[i, j, k] = c8[n]
    return C


def cube_to_c8(C):
    return tuple(int(C[i, j, k]) for i, j, k in OFF)


def face_slice(X, axis, side):
    return np.take(X, 0 if side == 0 else -1, axis=axis)


# ------------------------------------------------------------------------------------------------ I/O
def read_blocks_vtk(path, field="block"):
    """hex VTK with a per-cell block id -> points, list of structured node-id grids (one per block id).
    The (i,j,k) indexing is recovered by walking face neighbours (no assumption on the writer)."""
    import meshio
    m = meshio.read(path)
    H = np.vstack([c.data for c in m.cells if c.type == "hexahedron"]).astype(np.int64)
    bid = np.concatenate([d for d, c in zip(m.cell_data[field], m.cells) if c.type == "hexahedron"]).astype(int)
    P = np.asarray(m.points, float)
    LOC = OFF
    nbr_local = {}          # local node -> its 3 edge neighbours (local)
    for a, b in HEX_EDGES:
        nbr_local.setdefault(a, []).append(b); nbr_local.setdefault(b, []).append(a)
    grids = []
    for b in np.unique(bid):
        cells = np.where(bid == b)[0]
        f2c = collections.defaultdict(list)
        for c in cells:
            for ax in range(3):
                for s in (0, 1):
                    loc = [n for n in range(8) if LOC[n][ax] == s]
                    f2c[frozenset(H[c, loc].tolist())].append(c)
        ijk = {}
        c0 = cells[0]
        for n in range(8):
            ijk[int(H[c0, n])] = LOC[n].copy()
        seen = {c0}; dq = collections.deque([c0])
        while dq:
            c = dq.popleft(); hc = H[c]
            for ax in range(3):
                for s in (0, 1):
                    loc = [n for n in range(8) if LOC[n][ax] == s]
                    key = frozenset(hc[loc].tolist())
                    for g in f2c[key]:
                        if g in seen:
                            continue
                        hg = H[g]; S = set(key)
                        lg = {int(v): i for i, v in enumerate(hg)}; lc = {int(v): i for i, v in enumerate(hc)}
                        for n in key:
                            m_ = [int(hc[x]) for x in nbr_local[lc[n]] if int(hc[x]) not in S][0]
                            p_ = [int(hg[x]) for x in nbr_local[lg[n]] if int(hg[x]) not in S][0]
                            ijk[p_] = 2 * ijk[n] - ijk[m_]
                        seen.add(g); dq.append(g)
        assert len(seen) == len(cells), f"block {b}: not face-connected"
        K = np.array(list(ijk.values())); lo = K.min(0); dims = K.max(0) - lo + 1
        G = -np.ones(dims, np.int64)
        for n, v in ijk.items():
            G[tuple(v - lo)] = n
        assert (G >= 0).all() and G.size == len(ijk) and np.prod(dims - 1) == len(cells), f"block {b} not structured"
        X = P[G]
        if np.median(hex_sj(X)) < 0:            # make right-handed
            G = G[:, :, ::-1]
        grids.append(G)
    return P, grids


def write_hex_vtk(fn, P, Hx, cdata, pdata=()):
    with open(fn, "wb") as f:
        f.write(b"# vtk DataFile Version 3.0\nblock_collapse\nBINARY\nDATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {len(P)} double\n".encode()); f.write(np.asarray(P, ">f8").tobytes()); f.write(b"\n")
        f.write(f"CELLS {len(Hx)} {9 * len(Hx)}\n".encode())
        f.write(np.c_[np.full(len(Hx), 8), Hx].astype(">i4").tobytes()); f.write(b"\n")
        f.write(f"CELL_TYPES {len(Hx)}\n".encode()); f.write(np.full(len(Hx), 12, ">i4").tobytes()); f.write(b"\n")
        f.write(f"CELL_DATA {len(Hx)}\n".encode())
        for name, arr in cdata:
            arr = np.asarray(arr); typ = "int" if arr.dtype.kind in "iub" else "double"
            f.write(f"SCALARS {name} {typ} 1\nLOOKUP_TABLE default\n".encode())
            f.write(arr.astype(">i4" if typ == "int" else ">f8").tobytes()); f.write(b"\n")
        if pdata:
            f.write(f"POINT_DATA {len(P)}\n".encode())
            for name, arr in pdata:
                arr = np.asarray(arr); typ = "int" if arr.dtype.kind in "iub" else "double"
                f.write(f"SCALARS {name} {typ} 1\nLOOKUP_TABLE default\n".encode())
                f.write(arr.astype(">i4" if typ == "int" else ">f8").tobytes()); f.write(b"\n")


def write_lines_vtk(fn, polylines, cdata):
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


# ------------------------------------------------------------------------------------------------ grids
def hex_sj(P):
    """min corner scaled Jacobian per hex of a structured grid (right-handed index order)"""
    c = {}
    for x in (0, 1):
        for y in (0, 1):
            for z in (0, 1):
                c[x, y, z] = P[x:P.shape[0] - 1 + x, y:P.shape[1] - 1 + y, z:P.shape[2] - 1 + z]
    sj = None
    for x in (0, 1):
        for y in (0, 1):
            for z in (0, 1):
                p = c[x, y, z]
                e1 = (c[1 - x, y, z] - p) * (1 if x == 0 else -1)
                e2 = (c[x, 1 - y, z] - p) * (1 if y == 0 else -1)
                e3 = (c[x, y, 1 - z] - p) * (1 if z == 0 else -1)
                d = np.einsum("...i,...i->...", np.cross(e1, e2), e3)
                s = d / (np.linalg.norm(e1, axis=-1) * np.linalg.norm(e2, axis=-1) * np.linalg.norm(e3, axis=-1) + 1e-300)
                sj = s if sj is None else np.minimum(sj, s)
    return sj


def arc_param(G, axis):
    """normalised parameter along `axis`: cumulative segment length averaged over all other lines"""
    D = np.linalg.norm(np.diff(G, axis=axis), axis=-1)
    oth = tuple(a for a in range(D.ndim) if a != axis)
    d = D.mean(axis=oth) if oth else D
    if d.sum() <= 1e-300:
        d = np.ones_like(d)
    d = np.maximum(d, 1e-6 * d.max())
    s = np.r_[0, np.cumsum(d)]
    return s / s[-1]


def resample_axis(G, axis, n):
    """resample a structured grid to n intervals along axis, uniform in the averaged arc-length parameter"""
    s = arc_param(G, axis)
    f = np.interp(np.linspace(0, 1, n + 1), s, np.arange(len(s)))
    return take_frac(G, axis, f)


def take_frac(G, axis, f):
    """G sampled at fractional indices f (1D array) along axis (linear)"""
    m = G.shape[axis] - 1
    i0 = np.clip(np.floor(f).astype(int), 0, max(m - 1, 0)); w = f - i0
    A = np.take(G, i0, axis=axis); B = np.take(G, np.minimum(i0 + 1, m), axis=axis)
    sh = [1] * G.ndim; sh[axis] = len(f)
    w = w.reshape(sh)
    return (1 - w) * A + w * B


def resample(G, counts):
    for ax, n in enumerate(counts):
        if G.shape[ax] - 1 != n:
            G = resample_axis(G, ax, n)
    return G


def sample_surface_frac(H, T):
    """H (nu+1, nv+1, nw+1, 3); T (nu+1, nv+1) fractional parameter in [0,1] along w -> (nu+1, nv+1, 3)"""
    nw = H.shape[2] - 1
    f = np.clip(T, 0, 1) * nw; i0 = np.clip(np.floor(f).astype(int), 0, max(nw - 1, 0)); w = (f - i0)[..., None]
    A = np.take_along_axis(H, i0[..., None, None], axis=2)[:, :, 0]
    B = np.take_along_axis(H, np.minimum(i0 + 1, nw)[..., None, None], axis=2)[:, :, 0]
    return (1 - w) * A + w * B


def coons2(e_u0, e_u1, e_v0, e_v1):
    """2D TFI: e_v0 = G[:,0], e_v1 = G[:,-1] (len nu+1); e_u0 = G[0,:], e_u1 = G[-1,:] (len nv+1)"""
    nu, nv = len(e_v0) - 1, len(e_u0) - 1
    u = np.linspace(0, 1, nu + 1)[:, None, None]; v = np.linspace(0, 1, nv + 1)[None, :, None]
    G = (1 - v) * e_v0[:, None] + v * e_v1[:, None] + (1 - u) * e_u0[None] + u * e_u1[None]
    G -= ((1 - u) * (1 - v) * e_v0[0] + u * (1 - v) * e_v0[-1] + (1 - u) * v * e_v1[0] + u * v * e_v1[-1])
    return G


def tfi_correct2(G, E):
    """E = dict with keys (axis, side) -> boundary line of the target; returns G corrected by the Coons patch
    of the boundary differences (G already at the target counts), boundary set exactly"""
    d = {k: E[k] - face_slice(G, *k) for k in E}
    G = G + coons2(d[(0, 0)], d[(0, 1)], d[(1, 0)], d[(1, 1)])
    for (a, s), L in E.items():
        if a == 0:
            G[0 if s == 0 else -1] = L
        else:
            G[:, 0 if s == 0 else -1] = L
    return G


def tfi3_full(P):
    n1, n2, n3 = (s - 1 for s in P.shape[:3])
    u = np.linspace(0, 1, n1 + 1)[:, None, None, None]; v = np.linspace(0, 1, n2 + 1)[None, :, None, None]
    w = np.linspace(0, 1, n3 + 1)[None, None, :, None]
    Pu = (1 - u) * P[:1] + u * P[-1:]
    Pv = (1 - v) * P[:, :1] + v * P[:, -1:]
    Pw = (1 - w) * P[:, :, :1] + w * P[:, :, -1:]
    PuPv = (1 - u) * ((1 - v) * P[:1, :1] + v * P[:1, -1:]) + u * ((1 - v) * P[-1:, :1] + v * P[-1:, -1:])
    PvPw = (1 - v) * ((1 - w) * P[:, :1, :1] + w * P[:, :1, -1:]) + v * ((1 - w) * P[:, -1:, :1] + w * P[:, -1:, -1:])
    PuPw = (1 - u) * ((1 - w) * P[:1, :, :1] + w * P[:1, :, -1:]) + u * ((1 - w) * P[-1:, :, :1] + w * P[-1:, :, -1:])
    T = Pu + Pv + Pw - PuPv - PvPw - PuPw
    for a, ia in ((1 - u, 0), (u, -1)):
        for b, ib in ((1 - v, 0), (v, -1)):
            for c, ic in ((1 - w, 0), (w, -1)):
                T = T + a * b * c * P[ia, ib, ic][None, None, None]
    return T


def tfi_correct3(G, F):
    """F = dict (axis, side) -> target face grid; G at target counts"""
    D = np.zeros_like(G)
    for (a, s), S in F.items():
        idx = [slice(None)] * 3; idx[a] = 0 if s == 0 else -1
        D[tuple(idx)] = S - G[tuple(idx)]
    G = G + tfi3_full(D)
    for (a, s), S in F.items():
        idx = [slice(None)] * 3; idx[a] = 0 if s == 0 else -1
        G[tuple(idx)] = S
    return G


def orient(G, C, target):
    """re-index grid G (nd structured, + coord axis) and its corner-id array C (shape (2,)*nd) by an axis
    permutation + flips so that C matches `target` (same shape, -1 = don't care). Returns (G', C') or None."""
    nd = C.ndim
    for perm in itertools.permutations(range(nd)):
        Cp = np.transpose(C, perm); Gp = np.transpose(G, perm + (nd,))
        for fl in itertools.product((False, True), repeat=nd):
            Cf, Gf = Cp, Gp
            for a, f in enumerate(fl):
                if f:
                    Cf = np.flip(Cf, a); Gf = np.flip(Gf, a)
            if np.all((target < 0) | (Cf == target)):
                return np.ascontiguousarray(Gf), Cf.copy()
    return None


def canon2(G, C):
    """canonical orientation of a face grid: C[0,0] = min id, C[1,0] < C[0,1]"""
    i, j = np.unravel_index(np.argmin(C), C.shape)
    T = -np.ones((2, 2), np.int64); T[0, 0] = C[i, j]
    a, b = C[1 - i, j], C[i, 1 - j]
    T[1, 0], T[0, 1] = min(a, b), max(a, b)
    return orient(G, C, T)


# ------------------------------------------------------------------------------------------------ geometry
def closest_on_tris(p, a, b, c):
    """Ericson closest point; p (..., 3), a/b/c broadcastable"""
    ab, ac, ap = b - a, c - a, p - a
    d1 = (ab * ap).sum(-1); d2 = (ac * ap).sum(-1)
    bp = p - b; d3 = (ab * bp).sum(-1); d4 = (ac * bp).sum(-1)
    cp = p - c; d5 = (ab * cp).sum(-1); d6 = (ac * cp).sum(-1)
    va = d3 * d6 - d5 * d4; vb = d5 * d2 - d1 * d6; vc = d1 * d4 - d3 * d2
    den = va + vb + vc; den = np.where(np.abs(den) < 1e-300, 1e-300, den)
    v = vb / den; w = vc / den
    res = a + ab * v[..., None] + ac * w[..., None]
    res = np.broadcast_to(res, np.broadcast_shapes(res.shape, p.shape)).copy()
    A_, B_, C_ = (np.broadcast_to(x, res.shape) for x in (a, b, c))
    def put(m, val):
        res[m] = np.broadcast_to(val, res.shape)[m]
    put((d1 <= 0) & (d2 <= 0), A_)
    put((d3 >= 0) & (d4 <= d3), B_)
    put((d6 >= 0) & (d5 <= d6), C_)
    m = (vc <= 0) & (d1 >= 0) & (d3 <= 0); t = d1 / np.where(d1 - d3 == 0, 1, d1 - d3); put(m, a + ab * t[..., None])
    m = (vb <= 0) & (d2 >= 0) & (d6 <= 0); t = d2 / np.where(d2 - d6 == 0, 1, d2 - d6); put(m, a + ac * t[..., None])
    m = (va <= 0) & ((d4 - d3) >= 0) & ((d5 - d6) >= 0)
    t = (d4 - d3) / np.where((d4 - d3) + (d5 - d6) == 0, 1, (d4 - d3) + (d5 - d6)); put(m, b + (c - b) * t[..., None])
    return res


class TriSurface:
    """triangle set with closest-point queries (KD tree on centroids, exact point-triangle on k nearest)"""
    def __init__(self, A, B, C, k=12):
        self.A, self.B, self.C = A, B, C; self.k = min(k, len(A))
        self.tree = cKDTree((A + B + C) / 3)
        # candidate triangles also via their vertices (centroids alone miss long, thin boundary-layer triangles)
        self.vtree = cKDTree(np.vstack([A, B, C])); self.ntri = len(A)

    def project(self, Q):
        Q = np.asarray(Q, float); sh = Q.shape; Q = Q.reshape(-1, 3)
        if not len(Q):
            return Q.reshape(sh), np.zeros(sh[:-1])
        _, nn = self.tree.query(Q, k=self.k)
        _, nv = self.vtree.query(Q, k=min(self.k, 3 * self.ntri))
        nn = np.concatenate([nn.reshape(len(Q), -1), nv.reshape(len(Q), -1) % self.ntri], 1)
        R = closest_on_tris(Q[:, None], self.A[nn], self.B[nn], self.C[nn])
        d = np.linalg.norm(R - Q[:, None], axis=-1); j = d.argmin(1)
        out = R[np.arange(len(Q)), j]
        return out.reshape(sh), d[np.arange(len(Q)), j].reshape(sh[:-1])


class Polylines:
    """closest point on a set of polylines (segments between consecutive samples of the same id)"""
    def __init__(self, P, ids):
        ok = ids[1:] == ids[:-1]
        self.a = P[:-1][ok]; self.b = P[1:][ok]; self.id = ids[:-1][ok]
        self.tree = cKDTree((self.a + self.b) / 2)
        self.seglen = np.linalg.norm(self.b - self.a, axis=1).max()

    def project(self, Q, only=None, k=8):
        Q = np.asarray(Q, float).reshape(-1, 3)
        if only is not None:
            if not hasattr(self, "_cache"):
                self._cache = {}
            if only not in self._cache:
                sel = np.where(self.id == only)[0]
                self._cache[only] = (sel, cKDTree((self.a[sel] + self.b[sel]) / 2))
            sel, tr = self._cache[only]
            _, nn = tr.query(Q, k=min(k, len(sel)))
            nn = sel[nn.reshape(len(Q), -1)]
        else:
            _, nn = self.tree.query(Q, k=k); nn = nn.reshape(len(Q), -1)
        a, b = self.a[nn], self.b[nn]; ab = b - a
        t = np.clip(((Q[:, None] - a) * ab).sum(-1) / np.maximum((ab * ab).sum(-1), 1e-300), 0, 1)
        R = a + ab * t[..., None]; d = np.linalg.norm(R - Q[:, None], axis=-1); j = d.argmin(1)
        r = np.arange(len(Q))
        return R[r, j], d[r, j], self.id[nn[r, j]]


class CAD:
    """CAD export (cad_export.py): samples, topology, closure relation"""
    def __init__(self, path):
        D = np.load(path, allow_pickle=False)
        self.fp, self.fid = D["face_pts"], D["face_id"]
        self.ep, self.eid = D["edge_pts"], D["edge_id"]
        self.V = D["verts"]
        t = json.loads(str(D["topo"]))
        self.edge_faces = [set(x) for x in t["edge_faces"]]
        self.vert_edges = [set(x) for x in t["vert_edges"]]
        self.nf = int(self.fid.max()) + 1
        self.ftree = cKDTree(self.fp); self.vtree = cKDTree(self.V)
        self.edges = Polylines(self.ep, self.eid)
        self.edge_verts = collections.defaultdict(set)
        for v, es in enumerate(self.vert_edges):
            for e in es:
                self.edge_verts[e].add(v)

    def label_points(self, Q, k=8, cand=4):
        """CAD face of each point: distance to the local tangent plane (PCA of the k nearest samples) of each of
        the `cand` nearest faces; robust where a thin strip of one face lies next to another face"""
        if not hasattr(self, "_ftrees"):
            self._ftrees = {f: cKDTree(self.fp[self.fid == f]) for f in np.unique(self.fid)}
            self._fpts = {f: self.fp[self.fid == f] for f in self._ftrees}
        d0, j0 = self.ftree.query(Q, k=64)
        best_d = np.full(len(Q), np.inf); best_f = self.fid[j0[:, 0]].copy()
        for f in self._ftrees:
            m = (self.fid[j0] == f).any(1)
            if not m.any():
                continue
            kk = min(k, len(self._fpts[f]))
            dn, jn = self._ftrees[f].query(Q[m], k=kk)
            S = self._fpts[f][jn.reshape(m.sum(), -1)]
            c = S.mean(1); U = S - c[:, None]
            _, _, Vt = np.linalg.svd(U, full_matrices=False)
            n = Vt[:, -1]
            dist = np.abs(((Q[m] - c) * n).sum(-1))
            spacing = dn.reshape(m.sum(), -1)[:, -1]
            dist = np.where(dn.reshape(m.sum(), -1)[:, 0] <= 1.5 * spacing + 1e-9, dist, np.inf)
            dist = np.maximum(dist, 0) + 1e-3 * dn.reshape(m.sum(), -1)[:, 0]
            sel = np.where(m)[0]; better = dist < best_d[sel]
            best_d[sel[better]] = dist[better]; best_f[sel[better]] = f
        return best_f, best_d

    def contains(self, lo, hi):
        """is CAD entity lo=(dim,id) in the closure of hi=(dim,id)? (dim 3 = interior contains all)"""
        (dl, il), (dh, ih) = lo, hi
        if dh == 3:
            return True
        if dl > dh:
            return False
        if dl == dh:
            return il == ih
        if dl == 0 and dh == 1:
            return il in self.edge_verts[ih] if ih >= 0 else False
        if dl == 0 and dh == 2:
            return any(ih in self.edge_faces[e] for e in self.vert_edges[il])
        if dl == 1 and dh == 2:
            return ih in self.edge_faces[il]
        return False


def merge_class(cad, a, b):
    """class of an entity obtained by merging two entities of classes a, b; None if incompatible"""
    lo, hi = (a, b) if a[0] <= b[0] else (b, a)
    return lo if cad.contains(lo, hi) else None


# ------------------------------------------------------------------------------------------------ topology
def coarse_hexes(blocks):
    return [tuple(b["c8"]) for b in blocks]


def block_entities(c8):
    C = corner_cube(c8)
    faces = {}
    for a, s in FACES:
        faces[(a, s)] = frozenset(np.take(C, s, axis=a).ravel().tolist())
    edges = []
    for ax, grp in enumerate(PARALLEL):
        for i, j in grp:
            edges.append((ax, frozenset((c8[i], c8[j]))))
    return faces, edges


def sheets_of(c8s):
    """dual sheets of a coarse hex mesh: equivalence classes of parallel block edges.
    Returns list of dict(edges=set of frozenset, members=[(block, axis)])"""
    par = {}
    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for c8 in c8s:
        for grp in PARALLEL:
            es = [frozenset((c8[i], c8[j])) for i, j in grp]
            for e in es:
                par.setdefault(e, e)
            for e in es[1:]:
                ra, rb = find(es[0]), find(e)
                if ra != rb:
                    par[rb] = ra
    cls = collections.defaultdict(set)
    for e in par:
        cls[find(e)].add(e)
    out = []
    for r, es in cls.items():
        mem = []
        for bi, c8 in enumerate(c8s):
            for ax, grp in enumerate(PARALLEL):
                i, j = grp[0]
                if frozenset((c8[i], c8[j])) in es:
                    mem.append((bi, ax))
        out.append(dict(edges=es, members=mem))
    out.sort(key=lambda s: (-len(s["members"]), min(min(e) for e in s["edges"])))
    return out


# ------------------------------------------------------------------------------------------------ parametric resampling
def line_param(Lp, Lt):
    """fractional index on prior line Lp of the points of target curve Lt (both describe the same curve):
    matched by normalised arc length"""
    def s_(L):
        d = np.linalg.norm(np.diff(L, axis=0), axis=1); s = np.r_[0, np.cumsum(d)]
        return s / s[-1] if s[-1] > 0 else np.linspace(0, 1, len(L))
    return np.interp(s_(Lt), s_(Lp), np.arange(len(Lp), dtype=float))


def _gather(G, idx_f):
    """multilinear interpolation of grid G (nd + coord) at fractional index arrays idx_f (list of nd arrays)"""
    nd = len(idx_f); lo, w = [], []
    for a, f in enumerate(idx_f):
        m = G.shape[a] - 1
        i0 = np.clip(np.floor(f).astype(int), 0, max(m - 1, 0)); lo.append(i0); w.append(np.clip(f - i0, 0, 1) if m > 0 else 0 * f)
    out = 0
    for corner in itertools.product((0, 1), repeat=nd):
        ix = tuple(np.minimum(lo[a] + c, G.shape[a] - 1) for a, c in enumerate(corner))
        wt = 1
        for a, c in enumerate(corner):
            wt = wt * (w[a] if c else 1 - w[a])
        out = out + wt[..., None] * G[ix]
    return out


def resample_face_param(G, E):
    """face prior G -> grid matching target boundary curves E[(axis, side)] (counts from E):
    boundary parameters by arc-length matching, interior parameters by linear blending (monotone)"""
    eb, et = E[(1, 0)], E[(1, 1)]          # along axis 0 at v = 0 / 1
    el, er = E[(0, 0)], E[(0, 1)]          # along axis 1 at u = 0 / 1
    ab = line_param(G[:, 0], eb); at = line_param(G[:, -1], et)
    bl = line_param(G[0, :], el); br = line_param(G[-1, :], er)
    n0, n1 = len(eb) - 1, len(el) - 1
    u = np.linspace(0, 1, n0 + 1)[:, None]; v = np.linspace(0, 1, n1 + 1)[None]
    U = (1 - v) * ab[:, None] + v * at[:, None]
    V = (1 - u) * bl[None] + u * br[None]
    return tfi_correct2(_gather(G, [U, V]), E)


def resample_block_param(X, edge_lines, F):
    """block prior X -> grid matching the target faces F; edge_lines[ax][(p, q)] = target curve of the edge
    along axis ax at corner (p, q) of the remaining axes (in index order)"""
    cnt = [len(edge_lines[ax][(0, 0)]) - 1 for ax in range(3)]
    grids = np.meshgrid(*[np.linspace(0, 1, n + 1) for n in cnt], indexing="ij")
    idx = []
    for ax in range(3):
        o = [a for a in range(3) if a != ax]
        par = {}
        for (p, q), Lt in edge_lines[ax].items():
            sl = [None] * 3; sl[ax] = slice(None); sl[o[0]] = 0 if p == 0 else -1; sl[o[1]] = 0 if q == 0 else -1
            par[(p, q)] = line_param(X[tuple(sl)], Lt)
        sh = [1, 1, 1]; sh[ax] = -1
        pv = {k: v.reshape(sh) for k, v in par.items()}
        s, t = grids[o[0]], grids[o[1]]
        idx.append((1 - s) * (1 - t) * pv[(0, 0)] + s * (1 - t) * pv[(1, 0)] + (1 - s) * t * pv[(0, 1)] + s * t * pv[(1, 1)])
    return tfi_correct3(_gather(X, idx), F)
