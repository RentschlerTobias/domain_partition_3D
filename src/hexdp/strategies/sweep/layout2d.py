"""2D block layout of the meridional section of the draft-tube quarter (baffle domain), shared by all phi.

Method (generic, no labels / face numbers used as decisions):
  1. Sections (section_rz.py) give per phi one outer loop (CAD edges, tagged by the 3D face they lie on) and
     open curves of the embedded baffles (zero-thickness slits). A slit whose end point lies on the loop is
     an attached slit and is inserted into the loop (out and back); the others are free slits = own loops.
  2. Mean section: every loop segment resampled by normalised arc length and averaged over phi.
  3. Rectilinear (2D polycube) abstraction on the mean section: every segment gets a direction in
     {+r,+z,-r,-z} = rounded tangent angle (with 10 deg hysteresis). Where the rounded direction changes
     inside a smooth CAD edge, a virtual corner is inserted (e.g. the bulb nose turns 90 deg). Corner
     valence v = 2 - (turn of the rounded direction)/90, slit tips v = 4.
  4. Parameter coordinates: vertices joined by an r-segment share z, by a z-segment share r (union-find);
     each group coordinate = mean physical coordinate of its members.
  5. Layout = full-extension separatrices (rays) from every vertex with v >= 2 (v-1 rays into the domain)
     until they hit the boundary; T-junction free by construction. Blocks = faces of the arrangement.
  6. 2D chord collapse: merge two parallel coordinate groups (= remove the chord of blocks between them)
     greedily, smallest physical offset first, if the block count drops, the parameter polygon stays valid
     and the re-fitted section meshes (Winslow, all checked phi) have no inverted quad and min sJ >= guard.
  7. Re-fit per phi: boundary nodes by normalised arc length per segment, interior nodes by Winslow
     smoothing in the index space of the (singularity-free) parameter grid.
"""
import json, glob, math, itertools, collections, sys, os, time
import numpy as np
import scipy.sparse as sp, scipy.sparse.linalg as spl

# Default data directory (sections/, layout_final.json): the working directory, as the scripts are run
# from their case directory. New code passes paths explicitly.
HERE = os.getcwd()
M_RES = 400            # samples per segment for the mean section


# ------------------------------------------------------------------ geometry helpers
def resample(P, n):
    P = np.asarray(P, float)
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
    t = np.linspace(0, s[-1], n)
    return np.c_[np.interp(t, s, P[:, 0]), np.interp(t, s, P[:, 1])]


def at_fraction(P, t):
    """points at normalised arc-length fractions t on polyline P"""
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]; s /= s[-1]
    return np.c_[np.interp(t, s, P[:, 0]), np.interp(t, s, P[:, 1])]


def area(P):
    return 0.5 * np.sum(P[:, 0] * np.roll(P[:, 1], -1) - np.roll(P[:, 0], -1) * P[:, 1])


# ------------------------------------------------------------------ 1. sections -> loops of segments
def build_loops(S, tol=1e-4):
    """loops: list of lists of segments {'pts', 'face', 'kind_start', 'slit'}; outer loop CCW (domain left)"""
    L = S["loops"][0]
    segs = [dict(pts=np.array(d["pts"]), face=d["face"], start="cad", slit=None) for d in L]
    allp = np.vstack([s["pts"][:-1] for s in segs])
    if area(allp) < 0:                                          # make CCW
        segs = [dict(s, pts=s["pts"][::-1]) for s in segs[::-1]]
    loops = [segs]
    for k, d in enumerate(S["slits"]):
        P = np.array(d["pts"])
        att = None
        for end in (0, -1):
            for i, s in enumerate(segs):
                if np.linalg.norm(s["pts"][0] - P[end]) < tol:
                    att = (i, end)
        if att is not None:
            i, end = att
            Q = P if end == 0 else P[::-1]
            segs[i:i] = [dict(pts=Q, face=d["face"], start="junction", slit=(k, 0)),
                         dict(pts=Q[::-1], face=d["face"], start="tip", slit=(k, 1))]
        else:
            loops.append([dict(pts=P, face=d["face"], start="tip", slit=(k, 0)),
                          dict(pts=P[::-1], face=d["face"], start="tip", slit=(k, 1))])
    return loops


def signature(loops):
    return [[(s["face"], s["start"], s["slit"]) for s in L] for L in loops]


def load_sections(pattern=f"{HERE}/sections/sec_*.json"):
    secs = []
    for f in sorted(glob.glob(pattern)):
        S = json.load(open(f)); secs.append((S["phi"], build_loops(S)))
    secs.sort(key=lambda x: x[0])
    # free slits: same orientation as in the neighbouring (previous) section (OCC returns arbitrary orientation)
    for k in range(1, len(secs)):
        for li in range(1, len(secs[k][1])):
            ref = secs[k - 1][1][li][0]["pts"]; cur = secs[k][1][li][0]["pts"]
            if np.linalg.norm(cur[0] - ref[0]) > np.linalg.norm(cur[0] - ref[-1]):
                L = secs[k][1][li]
                L[0]["pts"], L[1]["pts"] = L[1]["pts"], L[0]["pts"]
    sig0 = signature(secs[0][1])
    for phi, lp in secs:
        assert signature(lp) == sig0, f"section topology differs at phi {phi}"
    return secs


def mean_loops(secs):
    out = []
    for li, L in enumerate(secs[0][1]):
        ML = []
        for si, s in enumerate(L):
            P = np.mean([resample(lp[li][si]["pts"], M_RES) for _, lp in secs], axis=0)
            ML.append(dict(s, pts=P))
        out.append(ML)
    return out


# ------------------------------------------------------------------ 2. rectilinear abstraction
def tangent_angles(P):
    d = np.diff(P, axis=0); a = np.degrees(np.arctan2(d[:, 1], d[:, 0]))
    return np.degrees(np.unwrap(np.radians(a)))


def abstract(mloops, hyst=10.0):
    """split segments at virtual corners; returns per loop a list of param edges:
       dict(seg=index of the source segment, t0, t1 (fractions on it), D=direction (deg, multiple of 90),
            vkind=kind of the start vertex, v=valence of the start vertex)"""
    res = []
    for L in mloops:
        edges = []
        A_prev = None; D = None
        for si, s in enumerate(L):
            A = tangent_angles(s["pts"])
            if A_prev is None:
                D = 90.0 * round(A[:5].mean() / 90.0); A = A - 360 * round((A[0] - D) / 360)
                turn0 = None
            else:
                if s["start"] == "tip":
                    turn = -180.0
                else:
                    turn = ((A[0] - A_prev + 180) % 360) - 180
                A = A - 360 * round((A[0] - (A_prev + turn)) / 360)
            if A_prev is not None:
                Dn = 90.0 * round(A[0] / 90.0) if s["start"] != "tip" else D - 180
                v = int(round(2 - (Dn - D) / 90))
                edges[-1]["t1"] = 1.0
                edges.append(dict(seg=si, t0=0.0, t1=None, D=Dn, vkind=s["start"], v=v)); D = Dn
            else:
                edges.append(dict(seg=si, t0=0.0, t1=None, D=D, vkind=s["start"], v=None))
            # virtual corners inside the segment
            frac = np.r_[0, np.cumsum(np.linalg.norm(np.diff(s["pts"], axis=0), axis=1))]; frac /= frac[-1]
            for k in range(len(A)):
                dev = A[k] - D
                if abs(dev) > 45 + hyst:
                    sgn = np.sign(dev)
                    kk = k
                    while kk > 0 and (A[kk - 1] - D) * sgn > 45:
                        kk -= 1
                    t = 0.5 * (frac[kk] + frac[kk + 1])
                    edges[-1]["t1"] = t
                    Dn = D + 90 * sgn
                    edges.append(dict(seg=si, t0=t, t1=None, D=Dn, vkind="virtual", v=int(round(2 - (Dn - D) / 90))))
                    D = Dn
            A_prev = A[-1]
        # closing vertex (start of the loop)
        s0 = L[0]; A0 = tangent_angles(s0["pts"])
        if s0["start"] == "tip":
            Dn = D - 180
        else:
            turn = ((A0[0] - A_prev + 180) % 360) - 180
            Dn = 90.0 * round((A_prev + turn) / 90.0)
        edges[-1]["t1"] = 1.0
        edges[0]["v"] = int(round(2 - (Dn - D) / 90))
        total = Dn - edges[0]["D"]
        res.append(dict(edges=edges, total_turn=total))
    return res


class UF:
    def __init__(s): s.p = {}
    def f(s, x):
        s.p.setdefault(x, x)
        while s.p[x] != x:
            s.p[x] = s.p[s.p[x]]; x = s.p[x]
        return x
    def u(s, a, b): s.p[s.f(a)] = s.f(b)


def dir_vec(D):
    k = int(round(D / 90)) % 4
    return [(1, 0), (0, 1), (-1, 0), (0, -1)][k]


class Layout:
    """parameter-space layout; vertices = (loop, edge index) start vertex of each param edge"""

    def __init__(self, mloops, absd, merges=()):
        self.mloops, self.absd = mloops, absd
        self.verts = []                       # (li, ei)
        self.vpos = []                        # mean physical position
        for li, a in enumerate(absd):
            for ei, e in enumerate(a["edges"]):
                self.verts.append((li, ei))
                P = mloops[li][e["seg"]]["pts"]
                self.vpos.append(at_fraction(P, np.array([e["t0"]]))[0])
        self.vpos = np.array(self.vpos)
        self.vid = {v: i for i, v in enumerate(self.verts)}
        self.merges = list(merges)
        self.build_groups()

    def edge_list(self):
        out = []
        for li, a in enumerate(self.absd):
            n = len(a["edges"])
            for ei, e in enumerate(a["edges"]):
                out.append((self.vid[(li, ei)], self.vid[(li, (ei + 1) % n)], e["D"], li, ei))
        return out

    def build_groups(self):
        ur, uz = UF(), UF()
        for i in range(len(self.verts)):
            ur.f(i); uz.f(i)
        for a, b, D, li, ei in self.edge_list():
            dx, dy = dir_vec(D)
            if dy == 0: uz.u(a, b)
            else: ur.u(a, b)
        n = len(self.verts); self.canon = list(range(n))
        for i in range(n):                    # vertices at the same place (attached-slit junction) = one point
            for j in range(i):
                if np.linalg.norm(self.vpos[i] - self.vpos[j]) < 1e-6:
                    ur.u(i, j); uz.u(i, j); self.canon[i] = self.canon[j]
        for ax, a, b in self.merges:
            (ur if ax == 0 else uz).u(a, b)
        self.ur, self.uz = ur, uz
        U = np.zeros((len(self.verts), 2))
        for ax, uf in ((0, ur), (1, uz)):
            g = collections.defaultdict(list)
            for i in range(len(self.verts)):
                g[uf.f(i)].append(i)
            for root, mem in g.items():
                U[mem, ax] = self.vpos[mem, ax].mean()
            setattr(self, f"groups{ax}", g)
        self.U = np.round(U, 9) + 0.0          # snap (axis groups: mean of ~1e-16 values -> exactly 0)

    def valid(self, eps=1e-3):
        for a, b, D, li, ei in self.edge_list():
            dx, dy = dir_vec(D)
            L = (self.U[b] - self.U[a]) @ np.array([dx, dy])
            o = (self.U[b] - self.U[a]) @ np.array([-dy, dx])
            if L < eps or abs(o) > 1e-9:
                return False, f"edge {a}->{b} length {L:.4f} offset {o:.2e}"
        # simple polygon: no two non-adjacent param edges intersect / overlap (except slit pairs)
        E = self.edge_list()
        segs = [(self.U[a], self.U[b], a, b) for a, b, *_ in E]
        for i in range(len(segs)):
            for j in range(i + 1, len(segs)):
                p0, p1, a0, b0 = segs[i]; q0, q1, a1, b1 = segs[j]
                C = self.canon
                if {C[a0], C[b0]} & {C[a1], C[b1]}:
                    continue
                if seg_touch(p0, p1, q0, q1):
                    return False, f"edges {a0}-{b0} and {a1}-{b1} touch"
        return True, ""

    # -------------------------------------------------------------- arrangement in coarse index space
    def arrangement(self):
        U = self.U
        R = np.unique(np.round(U[:, 0], 12)); Z = np.unique(np.round(U[:, 1], 12))
        ir = {x: k for k, x in enumerate(R)}; iz = {x: k for k, x in enumerate(Z)}
        vij = [(ir[round(u[0], 12)], iz[round(u[1], 12)]) for u in U]
        bedges = set(); bnodes = set(); seg_nodes = []
        slit_edges = {}                     # coarse edge -> slit id (for node duplication)
        E = self.edge_list()
        for a, b, D, li, ei in E:
            (i0, j0), (i1, j1) = vij[a], vij[b]
            dx, dy = dir_vec(D)
            nodes = [(i0 + k * dx, j0 + k * dy) for k in range(max(abs(i1 - i0), abs(j1 - j0)) + 1)]
            seg_nodes.append(nodes)
            for p, q in zip(nodes[:-1], nodes[1:]):
                bedges.add(ekey(p, q))
            bnodes |= set(nodes)
        self.R, self.Z, self.vij, self.seg_nodes = R, Z, vij, seg_nodes
        # domain cells: centre inside the outer loop polygon (loop 0, slits have zero area)
        poly = np.array([self.U[a] for a, b, D, li, ei in E if li == 0])
        cells = set()
        for i in range(len(R) - 1):
            for j in range(len(Z) - 1):
                c = ((R[i] + R[i + 1]) / 2, (Z[j] + Z[j + 1]) / 2)
                if point_in_poly(c, poly):
                    cells.add((i, j))
        self.cells = cells
        # rays
        rays = []; redges = set()
        n_per_loop = [len(a["edges"]) for a in self.absd]
        for vi, (li, ei) in enumerate(self.verts):
            e = self.absd[li]["edges"][ei]; v = e["v"]
            if v < 2:
                continue
            Dout = e["D"]
            for m in range(1, v):
                dx, dy = dir_vec(Dout + 90 * m)
                p = vij[vi]; path = [p]
                while True:
                    q = (p[0] + dx, p[1] + dy)
                    if not self.edge_in_domain(p, q) or ekey(p, q) in bedges:
                        raise RuntimeError(f"ray from vertex {vi} leaves the domain at {p}->{q}")
                    redges.add(ekey(p, q)); path.append(q); p = q
                    if q in bnodes:
                        break
                rays.append((vi, Dout + 90 * m, path))
        self.rays = rays
        blockedges = bedges | redges
        # flood fill blocks
        comp = {}
        nb = 0
        for c in sorted(cells):
            if c in comp:
                continue
            stack = [c]; comp[c] = nb
            while stack:
                i, j = stack.pop()
                for (di, dj), ed in (((1, 0), ekey((i + 1, j), (i + 1, j + 1))), ((-1, 0), ekey((i, j), (i, j + 1))),
                                     ((0, 1), ekey((i, j + 1), (i + 1, j + 1))), ((0, -1), ekey((i, j), (i + 1, j)))):
                    n = (i + di, j + dj)
                    if n in cells and n not in comp and ed not in blockedges:
                        comp[n] = nb; stack.append(n)
            nb += 1
        blocks = []
        for b in range(nb):
            cs = [c for c, k in comp.items() if k == b]
            i0 = min(c[0] for c in cs); i1 = max(c[0] for c in cs) + 1
            j0 = min(c[1] for c in cs); j1 = max(c[1] for c in cs) + 1
            if (i1 - i0) * (j1 - j0) != len(cs):
                raise RuntimeError(f"block {b} not a rectangle")
            blocks.append((i0, i1, j0, j1))
        self.blocks, self.comp, self.blockedges, self.bedges = blocks, comp, blockedges, bedges
        self.check_conforming()
        return len(blocks)

    def edge_in_domain(self, p, q):
        """coarse edge p-q (unit) has a domain cell on both sides and is not a boundary edge"""
        (i0, j0), (i1, j1) = p, q
        if j0 == j1:
            i = min(i0, i1); a, b = (i, j0 - 1), (i, j0)
        else:
            j = min(j0, j1); a, b = (i0 - 1, j), (i0, j)
        return a in self.cells and b in self.cells

    def check_conforming(self):
        """T-junction = a non-corner node on a block side that is the end of a perpendicular block edge on the
        other side of that side, where the side is interior (not boundary / slit)"""
        bad = 0
        for i0, i1, j0, j1 in self.blocks:
            sides = [((i, j0) for i in range(i0 + 1, i1)), ((i, j1) for i in range(i0 + 1, i1)),
                     ((i0, j) for j in range(j0 + 1, j1)), ((i1, j) for j in range(j0 + 1, j1))]
            outs = [(0, -1), (0, 1), (-1, 0), (1, 0)]
            for nodes, (di, dj) in zip(sides, outs):
                for c in nodes:
                    q = (c[0] + di, c[1] + dj)
                    along = ((c[0] - 1, c[1]), c) if dj != 0 else ((c[0], c[1] - 1), c)
                    if ekey(c, q) in self.blockedges and ekey(*along) not in self.bedges and self.edge_in_domain(c, q):
                        bad += 1
        self.tjunctions = bad
        return bad


def ekey(p, q):
    return (min(p, q), max(p, q))


def seg_touch(p0, p1, q0, q1, eps=1e-9):
    """axis-aligned segments touch (intersect or overlap)"""
    ax0 = min(p0[0], p1[0]) - eps; bx0 = max(p0[0], p1[0]) + eps
    ay0 = min(p0[1], p1[1]) - eps; by0 = max(p0[1], p1[1]) + eps
    ax1 = min(q0[0], q1[0]); bx1 = max(q0[0], q1[0]); ay1 = min(q0[1], q1[1]); by1 = max(q0[1], q1[1])
    return not (bx1 < ax0 or ax1 > bx0 or by1 < ay0 or ay1 > by0)


def point_in_poly(c, poly):
    x, y = c; inside = False; n = len(poly)
    for k in range(n):
        x0, y0 = poly[k]; x1, y1 = poly[(k + 1) % n]
        if (y0 > y) != (y1 > y):
            xi = x0 + (y - y0) * (x1 - x0) / (y1 - y0)
            if xi > x:
                inside = not inside
    return inside


# ------------------------------------------------------------------ fine parameter grid + Winslow re-fit
class Mesh2D:
    def __init__(self, lay, h, nmin=1, extra_cuts=None, extra_b=None, spacing=None):
        """extra_cuts(mesh) -> {(I,J): ('z'|'r', index)}: additional fine nodes split by side (used by the C-grid)"""
        self.lay = lay
        R, Z = lay.R, lay.Z
        self.nr = [max(nmin, int(round((R[k + 1] - R[k]) / h))) for k in range(len(R) - 1)]
        self.nz = [max(nmin, int(round((Z[k + 1] - Z[k]) / h))) for k in range(len(Z) - 1)]
        self.oi = np.r_[0, np.cumsum(self.nr)]; self.oj = np.r_[0, np.cumsum(self.nz)]
        # parameter coordinate of fine index
        if spacing is None:
            self.ur = np.concatenate([np.linspace(R[k], R[k + 1], self.nr[k] + 1)[:-1] for k in range(len(R) - 1)] + [[R[-1]]])
            self.uz = np.concatenate([np.linspace(Z[k], Z[k + 1], self.nz[k] + 1)[:-1] for k in range(len(Z) - 1)] + [[Z[-1]]])
        else:                 # spacing(axis, k, lo, hi) -> node coordinates of coarse interval k incl. both ends
            dr = [np.asarray(spacing(0, k, R[k], R[k + 1])) for k in range(len(R) - 1)]
            dz = [np.asarray(spacing(1, k, Z[k], Z[k + 1])) for k in range(len(Z) - 1)]
            self.nr = [len(d) - 1 for d in dr]; self.nz = [len(d) - 1 for d in dz]
            self.oi = np.r_[0, np.cumsum(self.nr)]; self.oj = np.r_[0, np.cumsum(self.nz)]
            self.ur = np.concatenate([d[:-1] for d in dr] + [[R[-1]]]); self.uz = np.concatenate([d[:-1] for d in dz] + [[Z[-1]]])
        E = lay.edge_list()
        # slit-interior fine nodes
        self.slit_nodes = {}
        for (a, b, D, li, ei), nodes in zip(E, lay.seg_nodes):
            e = lay.absd[li]["edges"][ei]
            if lay.mloops[li][e["seg"]]["slit"] is None:
                continue
            fn = self.fine_path(nodes)
            dx, dy = dir_vec(D)
            for (I, J) in fn[1:-1]:
                self.slit_nodes[(I, J)] = ("z", I) if dx == 0 else ("r", J)
        if extra_cuts is not None:
            self.slit_nodes.update(extra_cuts(self))
        # cells and node keys
        keys = {}; cells = []
        for (ci, cj) in lay.cells:
            for I in range(self.oi[ci], self.oi[ci + 1]):
                for J in range(self.oj[cj], self.oj[cj + 1]):
                    q = []
                    for (A, B) in ((I, J), (I + 1, J), (I + 1, J + 1), (I, J + 1)):
                        k = self.key(A, B, I + 0.5, J + 0.5)
                        if k not in keys:
                            keys[k] = len(keys)
                        q.append(keys[k])
                    cells.append(q + [lay.comp[(ci, cj)]])
        self.keys = keys; self.cells = np.array(cells)
        self.N = len(keys)
        self.P = np.array([[self.ur[k[0]], self.uz[k[1]]] for k in keys])        # parameter position
        # boundary nodes: (node, edge index, tau) per param edge
        self.bnd = []
        for ei_g, ((a, b, D, li, ei), nodes) in enumerate(zip(E, lay.seg_nodes)):
            fn = self.fine_path(nodes)
            dx, dy = dir_vec(D)
            cx, cy = -dy * 0.5, dx * 0.5                      # left side (domain side) offset
            u = np.array([[self.ur[I], self.uz[J]] for I, J in fn])
            d = np.r_[0, np.cumsum(np.linalg.norm(np.diff(u, axis=0), axis=1))]; tau = d / d[-1]
            ids = [keys[self.key(I, J, I + cx, J + cy)] for I, J in fn]
            self.bnd.append((li, ei, np.array(ids), tau))
        self.is_b = np.zeros(self.N, bool)
        for _, _, ids, _ in self.bnd:
            self.is_b[ids] = True
        if extra_b is not None:                    # additional Dirichlet fine nodes (all side copies)
            for k, i in self.keys.items():
                if k[:2] in extra_b(self):
                    self.is_b[i] = True
        self.build_stencil()

    def fine_path(self, nodes):
        out = []
        for (i0, j0), (i1, j1) in zip(nodes[:-1], nodes[1:]):
            I0, J0, I1, J1 = self.oi[i0], self.oj[j0], self.oi[i1], self.oj[j1]
            n = max(abs(I1 - I0), abs(J1 - J0))
            seg = [(I0 + (I1 - I0) * k // n, J0 + (J1 - J0) * k // n) for k in range(n + 1)]
            out += seg if not out else seg[1:]
        return out

    def key(self, I, J, ci, cj):
        s = self.slit_nodes.get((I, J))
        if s is None:
            return (I, J)
        if s[0] == "z":
            return (I, J, int(ci > s[1]))
        return (I, J, int(cj > s[1]))

    def build_stencil(self):
        inner = [k for k, i in self.keys.items() if not self.is_b[i]]
        self.inner = np.array([self.keys[k] for k in inner])
        nb = np.zeros((len(inner), 8), int)
        offs = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, 1), (1, -1), (-1, -1)]
        for r, k in enumerate(inner):
            I, J = k[:2]
            for c, (di, dj) in enumerate(offs):
                nb[r, c] = self.keys[self.key(I + di, J + dj, I, J)]
        self.nb = nb

    def boundary_positions(self, loops):
        X = np.full((self.N, 2), np.nan)
        for li, ei, ids, tau in self.bnd:
            e = self.lay.absd[li]["edges"][ei]
            P = loops[li][e["seg"]]["pts"]
            X[ids] = at_fraction(np.asarray(P), e["t0"] + tau * (e["t1"] - e["t0"]))
        return X

    def solve(self, loops, iters=4000, tol=1e-9, X0=None):
        X = self.boundary_positions(loops)
        inn, nb = self.inner, self.nb
        n = len(inn)
        if X0 is None:
            # Laplace in index space (sparse) as start
            pos = {g: r for r, g in enumerate(inn)}
            rows, cols, vals = [], [], []
            rhs = np.zeros((n, 2))
            for r in range(n):
                rows.append(r); cols.append(r); vals.append(4.0)
                for c in range(4):
                    g = nb[r, c]
                    if g in pos:
                        rows.append(r); cols.append(pos[g]); vals.append(-1.0)
                    else:
                        rhs[r] += X[g]
            A = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
            lu = spl.splu(A.tocsc())
            X[inn] = np.c_[lu.solve(rhs[:, 0]), lu.solve(rhs[:, 1])]
        else:
            X[inn] = X0[inn]
        scale = np.ptp(X[self.is_b], axis=0).max()
        for it in range(iters):
            E, W, Nn, S, NE, NW, SE, SW = (X[nb[:, c]] for c in range(8))
            xi = 0.5 * (E - W); eta = 0.5 * (Nn - S)
            a = (eta ** 2).sum(1)[:, None]; g = (xi ** 2).sum(1)[:, None]; b = (xi * eta).sum(1)[:, None]
            new = (a * (E + W) + g * (Nn + S) - 0.5 * b * (NE - NW - SE + SW)) / (2 * (a + g) + 1e-300)
            d = np.abs(new - X[inn]).max()
            X[inn] = X[inn] + 1.0 * (new - X[inn])
            if d < tol * scale:
                break
        self.iters = it + 1
        return X

    def quality(self, X):
        """per cell min scaled Jacobian (corner-based), orientation from the parameter cell (CCW)"""
        q = self.cells[:, :4]; P = X[q]
        sj = np.full(len(q), np.inf)
        for c in range(4):
            p0 = P[:, c]; p1 = P[:, (c + 1) % 4]; p3 = P[:, (c + 3) % 4]
            e1 = p1 - p0; e2 = p3 - p0
            cr = e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0]
            sj = np.minimum(sj, cr / (np.linalg.norm(e1, axis=1) * np.linalg.norm(e2, axis=1) + 1e-300))
        return sj


# ------------------------------------------------------------------ plotting
def block_edge_polylines(m, X):
    """fine-mesh polylines of all block edges (cell sides on block boundaries) as segments"""
    lay = m.lay; segs = []
    blk = m.cells[:, 4]
    edge_cells = collections.defaultdict(list)
    for c, q in enumerate(m.cells[:, :4]):
        for k in range(4):
            a, b = q[k], q[(k + 1) % 4]
            edge_cells[(min(a, b), max(a, b))].append(blk[c])
    for (a, b), bl in edge_cells.items():
        if len(bl) == 1 or bl[0] != bl[1]:
            segs.append((X[a], X[b]))
    return np.array(segs)


def plot_mesh(m, X, ax, sj=None, fine=True, title=""):
    from matplotlib.collections import LineCollection, PolyCollection
    q = m.cells[:, :4]
    if fine:
        ax.add_collection(PolyCollection(X[q], facecolors="none", edgecolors="0.75", linewidths=0.2))
    if sj is not None and (sj <= 0).any():
        ax.add_collection(PolyCollection(X[q[sj <= 0]], facecolors="red", edgecolors="red", linewidths=0.5))
    ax.add_collection(LineCollection(block_edge_polylines(m, X), colors="k", linewidths=0.9))
    ax.autoscale(); ax.set_aspect("equal"); ax.set_title(title, fontsize=9)
