"""C/O-grid layer around the real thin walls, inserted at block level into the swept section layout (§9).

Generic construction (no labels / face numbers):
  * chains: for every slit (baffle mid-curve) of the layout, the chain of layout edges along it: a free slit gives a
    closed loop (both sides, O-layer); an attached slit gives an open chain (C-layer) that is extended along
    collinear boundary edges (v = 2 vertices) which lie inside the real wall (mid-curve of the obstacle wall, here
    the lip/mantle). Chain ends:
      - end vertex v = 1 (convex corner): the layer ends on the adjacent boundary edge (end edge on the boundary);
      - end vertex v = 3 (reflex corner): the layer cannot end on the boundary -> it continues as a strip ("H")
        along the straight-ahead separatrix of that vertex until that separatrix hits the boundary.
  * pillowing: the chain nodes of the existing 2D grid are moved from the mid-curve to the outer layer boundary
    W + delta*n (W = wall point by normalised arc length between the wall noses, n = smoothed wall normal); the
    layer (n_L cells, geometric grading from the wall) is inserted between the wall and these nodes; the outside
    grid is re-fitted by Winslow with the moved nodes as Dirichlet data. Every separatrix ending on the slit
    (incl. the 3 at each free tip) continues as a radial layer edge to the wall -> conforming, T-junction free.
    The old tips become interior vertices of valence 6 (3 rays + 2 layer sides + 1 radial layer edge).
  * layer / strip blocks: split wherever the adjacent outside block changes along the chain.
"""
import sys, os, json, glob, collections
import numpy as np
sys.path.insert(0, "/work/canada/blocking_study/sweep_phi")
from layout2d import *

CG = os.path.dirname(os.path.abspath(__file__))
DELTA = 0.03          # layer thickness
N_L = 12              # layer cells
TIP_TANGENT = True
CORNER_SNAP = False
TIP_BISECT = True
UNTANGLE = True
H1_REL = 1e-3         # first cell = H1_REL * mean free-ring chord


def polylen(P):
    return np.r_[0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]


def at_s(P, s):
    c = polylen(P)
    return np.c_[np.interp(s, c, P[:, 0]), np.interp(s, c, P[:, 1])]


def closest(P, q):
    """arc length of the point of polyline P closest to q, and the point"""
    a, b = P[:-1], P[1:]; d = b - a; L2 = (d ** 2).sum(1) + 1e-300
    t = np.clip(((q - a) * d).sum(1) / L2, 0, 1); X = a + t[:, None] * d
    k = np.argmin(np.linalg.norm(X - q, axis=1)); c = polylen(P)
    return c[k] + t[k] * np.sqrt(L2[k]), X[k]


def grading(h1, delta, n):
    lo, hi = 0.5, 2.0
    for _ in range(200):
        q = 0.5 * (lo + hi)
        s = h1 * n if abs(q - 1) < 1e-12 else h1 * (q ** n - 1) / (q - 1)
        lo, hi = (q, hi) if s < delta else (lo, q)
    d = np.r_[0, np.cumsum(h1 * q ** np.arange(n))]
    return d * delta / d[-1], q


class CGrid:
    def __init__(self, h=0.1):
        self.secs = load_sections()
        self.phis = np.array([p for p, _ in self.secs])
        secs5 = [s for s in self.secs if abs(s[0] / 5 - round(s[0] / 5)) < 1e-9]
        ml = mean_loops(secs5); ab = abstract(ml)
        merges = [tuple(x) for x in json.load(open(f"{HERE}/layout_final.json"))["merges"]]
        self.lay = lay = Layout(ml, ab, merges); self.nb2 = lay.arrangement()
        self.E = E = lay.edge_list()
        self.walls = {phi: json.load(open(f"{CG}/walls/wall_{phi:05.1f}.json"))["walls"] for phi in self.phis}
        self.wmeta = json.load(open(f"{CG}/walls.json"))
        self.m0 = Mesh2D(lay, h)                         # slit model (for the strip separatrix positions)
        self.find_chains()
        self.m = Mesh2D(lay, h, extra_cuts=self.extra_cuts, extra_b=self.extra_bnodes)
        self.setup()

    # ------------------------------------------------------------ chains (layout level)
    def slit_of(self, g):
        a, b, D, li, ei = self.E[g]
        return self.lay.mloops[li][self.lay.absd[li]["edges"][ei]["seg"]]["slit"]

    def vval(self, vid):
        li, ei = self.lay.verts[vid]; return self.lay.absd[li]["edges"][ei]["v"]

    def wall_for(self, pts):
        """index of the wall component nearest to the polyline pts (phi = 45 section)"""
        W = self.walls[45.0]
        return int(np.argmin([np.min(np.linalg.norm(np.array(w["pts"])[::5, None] - pts[None, ::10], axis=2)) for w in W]))

    def edge_pts(self, g, k_phi):
        a, b, D, li, ei = self.E[g]; e = self.lay.absd[li]["edges"][ei]
        P = np.asarray(self.secs[k_phi][1][li][e["seg"]]["pts"])
        c = polylen(P); c /= c[-1]
        s = np.linspace(e["t0"], e["t1"], 200)
        return np.c_[np.interp(s, c, P[:, 0]), np.interp(s, c, P[:, 1])]

    def find_chains(self):
        E = self.E; nE = len(E); k45 = int(np.argmin(abs(self.phis - 45)))
        slits = collections.defaultdict(list)
        for g in range(nE):
            s = self.slit_of(g)
            if s is not None:
                slits[s[0]].append(g)
        gap = {tuple(c["pair"]): c["max_gap"] for c in self.wmeta["components"]}
        self.chains = []
        for k, gs in sorted(slits.items()):
            li = E[gs[0]][3]
            mid = np.vstack([self.edge_pts(g, k45) for g in gs])
            wi = self.wall_for(mid); wpair = tuple(self.walls[45.0][wi]["pair"]); t = gap[wpair]
            if li != 0:                                   # free slit -> closed O-layer
                self.chains.append(dict(slit=k, wall=wi, closed=True, edges=sorted(gs), start=None, end=None, t=t))
                continue
            loop = [g for g in range(nE) if E[g][3] == 0]
            pos = {g: i for i, g in enumerate(loop)}
            first, last = min(gs, key=lambda g: pos[g]), max(gs, key=lambda g: pos[g])
            chain = [g for g in loop if pos[first] <= pos[g] <= pos[last]]
            W = np.array(self.walls[45.0][wi]["pts"])
            def walllike(g):                            # boundary edge inside the real wall (mid-curve of it)
                P = self.edge_pts(g, k45)
                return np.max([abs(closest(W, p)[0] * 0 + np.linalg.norm(closest(W, p)[1] - p)) for p in P[::10]]) < 1.5 * t
            # extend forward / backward over v = 2 vertices along wall-like boundary edges
            while True:
                nxt = loop[(pos[chain[-1]] + 1) % len(loop)]
                if self.vval(E[nxt][0]) == 2 and E[nxt][2] == E[chain[-1]][2] and walllike(nxt):
                    chain.append(nxt)
                else:
                    break
            while True:
                prv = loop[(pos[chain[0]] - 1) % len(loop)]
                if self.vval(E[chain[0]][0]) == 2 and E[prv][2] == E[chain[0]][2] and walllike(prv):
                    chain.insert(0, prv)
                else:
                    break
            ends = []
            for end_v, adj, Dlast in ((E[chain[0]][0], loop[(pos[chain[0]] - 1) % len(loop)], E[chain[0]][2] + 180),
                                      (E[chain[-1]][1], loop[(pos[chain[-1]] + 1) % len(loop)], E[chain[-1]][2])):
                v = self.vval(end_v)
                if v == 1:
                    ends.append(dict(kind="boundary", vertex=end_v, edge=adj))
                elif v == 3:
                    ray = [r for r in self.lay.rays if r[0] == end_v and int(round(r[1] - Dlast)) % 360 == 0]
                    assert len(ray) == 1, f"no straight-ahead separatrix at chain end vertex {end_v}"
                    ends.append(dict(kind="strip", vertex=end_v, edge=adj, ray=ray[0], D=Dlast))
                else:
                    raise RuntimeError(f"chain end vertex valence {v} not handled")
            self.chains.append(dict(slit=k, wall=wi, closed=False, edges=chain, start=ends[0], end=ends[1], t=t,
                                    junction=E[first][0]))

    def extra_bnodes(self, m):
        if not hasattr(self, "_eb"):
            self._eb = {tuple(x) for c in self.chains if not c["closed"] for e in (c["start"], c["end"])
                        if e["kind"] == "strip" for x in e["fine"]}
        return self._eb

    def extra_cuts(self, m):
        out = {}
        for c in self.chains:
            if c["closed"]:
                continue
            # junction of the attached slit: split like the slit interior
            a, b, D, li, ei = self.E[[g for g in c["edges"] if self.slit_of(g) is not None][0]]
            i, j = self.lay.vij[c["junction"]]; dx, dy = dir_vec(D)
            out[(m.oi[i], m.oj[j])] = ("z", m.oi[i]) if dx == 0 else ("r", m.oj[j])
            for e in (c["start"], c["end"]):
                if e["kind"] == "strip":
                    vi, Dr, path = e["ray"]; dx, dy = dir_vec(Dr)
                    fp = m.fine_path(path); e["fine"] = fp
                    for (I, J) in fp:
                        out[(I, J)] = ("z", I) if dx == 0 else ("r", J)
        return out

    # ------------------------------------------------------------ structure on the fine grid
    def setup(self):
        m = self.m; E = self.E
        self.node_of_edge = {}
        for g, (li, ei, ids, tau) in enumerate(m.bnd):
            self.node_of_edge[g] = ids
        e2c = collections.defaultdict(list)
        for c_, q in enumerate(m.cells[:, :4]):
            for k in range(4):
                a, b = q[k], q[(k + 1) % 4]; e2c[(min(a, b), max(a, b))].append(c_)
        self.e2c = e2c
        newid = m.N
        extra_b = set()
        for c in self.chains:
            ids = []
            for g in c["edges"]:
                x = list(self.node_of_edge[g]); ids += x if not ids else x[1:]
            if c["closed"]:
                assert ids[0] == ids[-1]; ids = ids[:-1]
            c["ids"] = np.array(ids)
            # chain fractions per original edge (for the wall mapping): edge index of every chain node
            # outside block adjacent to every chain edge -> layer block split points
            n = len(ids); nn = n if c["closed"] else n - 1
            blk = []
            for i in range(nn):
                a, b = ids[i], ids[(i + 1) % n]
                cs = e2c[(min(a, b), max(a, b))]; assert len(cs) == 1, cs
                blk.append(int(m.cells[cs[0], 4]))
            c["adj_block"] = blk
            # layer node ids: lay_ids[i, k], k = 0 wall ... N_L outer (= chain node)
            L = np.zeros((n, N_L + 1), int)
            for i in range(n):
                for k in range(N_L):
                    L[i, k] = newid; newid += 1
                L[i, N_L] = ids[i]
            c["lay_ids"] = L
            for e in (c["start"], c["end"]):
                if e is None or e["kind"] != "strip":
                    continue
                fp = e["fine"]; Dr = e["ray"][1]; dx, dy = dir_vec(Dr)
                # side of the cut where the wall continues: decided per phi geometrically; keys of both sides
                keys0 = [m.keys[(I, J, 0)] for I, J in fp]; keys1 = [m.keys[(I, J, 1)] for I, J in fp]
                e["keys"] = (np.array(keys0), np.array(keys1))
                extra_b |= set(keys0) | set(keys1)
                # strip interior nodes (rows 1..n, k = 1..N_L-1); row 0 = layer radial line at the chain end
                S = np.zeros((len(fp), N_L + 1), int)
                for j in range(1, len(fp)):
                    for k in range(1, N_L):
                        S[j, k] = newid; newid += 1
                e["strip_ids"] = S
                # outside blocks along the strip (left of the ray direction) for block splits
                blk = []
                for j in range(len(fp) - 1):
                    a, b = keys0[j], keys0[j + 1]
                    cs = e2c[(min(a, b), max(a, b))]; blk.append(int(m.cells[cs[0], 4]))
                e["adj_block"] = blk
        self.NT = newid
        # strip separatrix nodes become Dirichlet nodes of the outside grid
        m.is_b[list(extra_b)] = True
        m.build_stencil()
        # layer grading
        rings = [c for c in self.chains if c["closed"]]
        k45 = int(np.argmin(abs(self.phis - 45)))
        chord = np.mean([polylen(np.vstack([self.edge_pts(g, k) for g in c["edges"][:1]]))[-1]
                         for c in rings for k in range(len(self.phis))])
        self.h1 = H1_REL * chord
        self.dk, self.q = grading(self.h1, DELTA, N_L)
        self.mean_ring_chord = float(chord)

    # ------------------------------------------------------------ per-phi geometry
    def wall_normals(self, P, closed, sign_ref=None, win=None):
        win = win or 0.5 * DELTA
        c = polylen(P)
        if closed:
            Q = P[:-1] if np.linalg.norm(P[0] - P[-1]) < 1e-9 else P
            A = area(Q)
        T = np.zeros_like(P)
        for i in range(len(P)):
            lo = np.searchsorted(c, c[i] - win); hi = np.searchsorted(c, c[i] + win) - 1
            lo = max(0, min(lo, i - 1)); hi = min(len(P) - 1, max(hi, i + 1))
            if closed and (i == 0 or i == len(P) - 1):
                T[i] = P[1] - P[-2]
            else:
                T[i] = P[hi] - P[lo]
        T /= np.linalg.norm(T, axis=1)[:, None]
        N = np.c_[T[:, 1], -T[:, 0]]                          # right normal
        if closed and A < 0:
            N = -N
        return N, c

    def nose(self, W, tip, tdir, t):
        """wall point where the mid-curve leaves the wall: the farthest wall point along the tangent ray of the
        mid-curve tip among points within t/2 of that ray (flat caps: cap middle; knife edges: the apex)"""
        d = W - tip; proj = d @ tdir; cross = np.abs(d[:, 0] * tdir[1] - d[:, 1] * tdir[0])
        cand = np.where((cross < 0.5 * t) & (proj > -0.5 * t) & (proj < 10 * t))[0]
        if not len(cand):
            return int(np.argmin(np.linalg.norm(d, axis=1)))
        pm = proj[cand].max(); c2 = cand[proj[cand] > pm - 0.25 * t]
        return int(c2[np.argmin(cross[c2])])

    def remap_edge(self, X, g, loops, near_pt, target):
        """move the boundary nodes of layout edge g so that its end nearest to near_pt lands on target (on or just
        beyond the edge curve); the far end stays, the others are scaled along the curve"""
        a, b, D, li, ei = self.E[g]; ed = self.lay.absd[li]["edges"][ei]
        P = np.asarray(loops[li][ed["seg"]]["pts"])
        ids = self.node_of_edge[g]
        s_t, q = closest(P, target)
        if np.linalg.norm(q - target) > 1e-9:                     # target beyond the curve end: extend
            P = np.vstack([target, P]) if s_t < 0.5 * polylen(P)[-1] else np.vstack([P, target])
        sc = np.array([closest(P, p)[0] for p in X[ids]]); s_t, _ = closest(P, target)
        near = 0 if np.linalg.norm(X[ids[0]] - near_pt) < np.linalg.norm(X[ids[-1]] - near_pt) else -1
        s_near, s_far = sc[near], sc[-1 - near]
        X[ids] = at_s(P, s_far + (sc - s_far) * (s_t - s_far) / (s_near - s_far))

    def geometry(self, kphi, Xs=None):
        """all positions at section kphi: outside grid X (Dirichlet set, interior = slit-model start), layers, strips"""
        m, E = self.m, self.E
        phi, loops = self.secs[kphi]
        if Xs is None:
            Xs = self.m0.solve(loops, iters=6000)
        X = m.boundary_positions(loops)
        for key, i in m.keys.items():                 # warm start from the slit model
            if np.isnan(X[i, 0]):
                X[i] = Xs[self.m0.keys[key[:2]]] if key[:2] in self.m0.keys else Xs[self.m0.keys[key]]
        Wl = self.walls[phi]
        if not hasattr(self, "blunt_info"):
            self.blunt_info = {}
        LAY, STRIP = {}, {}
        for ci, c in enumerate(self.chains):
            W = np.array(Wl[c["wall"]]["pts"]); t = c["t"]
            W, cut = blunt(W, c["closed"], t)
            self.blunt_info[phi, ci] = cut
            ids = c["ids"]; Mid = X[ids].copy(); n = len(ids)
            if c["closed"]:
                if np.linalg.norm(W[0] - W[-1]) > 1e-9:
                    W = np.vstack([W, W[:1]])
                Nn, cw = self.wall_normals(W, True); Lw = cw[-1]
                tipA, tipB = 0, len(self.node_of_edge[c["edges"][0]]) - 1
                sT = []
                for ti, nb_ in ((tipA, Mid[1]), (tipB, Mid[tipB - 1])):
                    tdir = Mid[ti] - nb_; tdir /= np.linalg.norm(tdir)
                    sT.append(cw[self.nose(W[:-1], Mid[ti], tdir, t)])
                sA, sB = sT
                side0 = Mid[tipA:tipB + 1]
                tang = np.gradient(side0, axis=0); tang /= np.linalg.norm(tang, axis=1)[:, None] + 1e-300
                left = side0 + 0.5 * t * np.c_[-tang[:, 1], tang[:, 0]]
                f0 = polylen(side0); f0 /= f0[-1]
                cand = []
                for bb in (sB if sB > sA else sB + Lw, sB if sB < sA else sB - Lw):
                    P_ = at_s(W, (sA + f0 * (bb - sA)) % Lw)
                    cand.append((np.mean(np.linalg.norm(P_[1:-1] - left[1:-1], axis=1)), bb))
                b0 = min(cand)[1]
                side1 = np.vstack([Mid[tipB:], Mid[:1]]); f1 = polylen(side1); f1 /= f1[-1]
                b1 = sA + Lw if b0 > sA else sA - Lw
                s_nodes = np.r_[sA + f0 * (b0 - sA), (b0 + f1 * (b1 - b0))[1:-1]]
                snaps = {}
                if CORNER_SNAP:                      # wall end of the radial next to a tip -> cap / cut corner
                    for s_tip, nbi in ((sA, 1), (b1, n - 1), (b0, tipB - 1), (b0, tipB + 1)):
                        sc = corner_between(W, s_tip, s_nodes[nbi], t, True)
                        if sc is not None:
                            snaps[nbi] = at_s(W, np.array([sc % Lw]))[0]
                s_nodes = s_nodes % Lw
                Wp = at_s(W, s_nodes)
                Nw = np.c_[np.interp(s_nodes, cw, Nn[:, 0]), np.interp(s_nodes, cw, Nn[:, 1])]
                if TIP_TANGENT:                      # radial layer edge at a tip continues the mid-curve tangent
                    for ti, nb_ in ((tipA, Mid[1]), (tipB, Mid[tipB - 1])):
                        d = Mid[ti] - nb_; Nw[ti] = d / np.linalg.norm(d)
            else:
                st, en = c["start"], c["end"]
                gsl = [g for g in c["edges"] if self.slit_of(g) is not None]
                ti = 0
                for g in c["edges"]:
                    if g == gsl[1]:
                        break
                    ti += len(self.node_of_edge[g]) - 1
                tdir = Mid[ti] - Mid[ti - 1]; tdir /= np.linalg.norm(tdir)
                if np.linalg.norm(W[0] - Mid[0]) > np.linalg.norm(W[-1] - Mid[0]):
                    W = W[::-1]
                cw = polylen(W)
                s_st = 0.0 if st["kind"] == "boundary" else closest(W, Mid[0])[0]
                s_en = cw[-1] if en["kind"] == "boundary" else closest(W, Mid[-1])[0]
                s_tip = cw[self.nose(W, Mid[ti], tdir, t)]
                fA = polylen(Mid[:ti + 1]); fA /= fA[-1]; fB = polylen(Mid[ti:]); fB /= fB[-1]
                s_nodes = np.r_[s_st + fA * (s_tip - s_st), (s_tip + fB * (s_en - s_tip))[1:]]
                snaps = {}
                if CORNER_SNAP:
                    for nbi in (ti - 1, ti + 1):
                        sc = corner_between(W, s_nodes[ti], s_nodes[nbi], t, False)
                        if sc is not None:
                            snaps[nbi] = at_s(W, np.array([sc]))[0]
                Wp = at_s(W, s_nodes)
                Nn, _ = self.wall_normals(W, False)
                Nw = np.c_[np.interp(s_nodes, cw, Nn[:, 0]), np.interp(s_nodes, cw, Nn[:, 1])]
                sg = np.sign(np.median(np.sum((Wp - Mid) * Nw, axis=1)))
                Nw = Nw * (sg if sg != 0 else 1)
                if TIP_TANGENT:
                    Nw[ti] = tdir
            Nw /= np.linalg.norm(Nw, axis=1)[:, None]
            if TIP_BISECT:                           # tip radial = outward bisector of the two wall chords
                tips_ = (0, len(self.node_of_edge[c["edges"][0]]) - 1) if c["closed"] else (ti,)
                for t_ in tips_:
                    u1 = Wp[(t_ - 1) % n] - Wp[t_]; u2 = Wp[(t_ + 1) % n] - Wp[t_]
                    u1 /= np.linalg.norm(u1); u2 /= np.linalg.norm(u2)
                    bis = -(u1 + u2)
                    if np.linalg.norm(bis) > 1e-9:
                        bis /= np.linalg.norm(bis)
                        Nw[t_] = bis if bis @ Nw[t_] > 0 else -bis
            lay = Wp[:, None, :] + self.dk[None, :, None] * Nw[:, None, :]
            for i_, Wc in snaps.items():             # outer node kept, wall end on the corner
                O_ = lay[i_, -1].copy()
                lay[i_] = Wc[None, :] + (self.dk / DELTA)[:, None] * (O_ - Wc)[None, :]
            if not c["closed"]:
                for e, row in ((c["start"], 0), (c["end"], n - 1)):
                    g = e["edge"]
                    if e["kind"] == "boundary":
                        # layer end edge on the adjacent boundary curve, at arc distances d_k from the wall point
                        a, b, D, li, ei = E[g]; ed = self.lay.absd[li]["edges"][ei]
                        P = np.asarray(loops[li][ed["seg"]]["pts"]); s_w, _ = closest(P, Wp[row])
                        sgn = -1 if s_w > 0.5 * polylen(P)[-1] else 1
                        lay[row] = at_s(P, s_w + sgn * self.dk)
                        self.remap_edge(X, g, loops, Mid[row], lay[row, -1])
                    else:
                        self.remap_edge(X, g, loops, Mid[row], Wp[row])
                        e["Wend"], e["Nend"] = Wp[row], Nw[row]
            X[ids] = lay[:, -1]
            LAY[ci] = lay
            if not c["closed"]:
                for e in (c["start"], c["end"]):
                    if e["kind"] == "strip":
                        STRIP[ci, "start" if e is c["start"] else "end"] = self.strip(e, X, Xs, loops)
        return X, LAY, STRIP, Xs

    def strip(self, e, X, Xs, loops):
        m = self.m; fp = e["fine"]; k0, k1 = e["keys"]
        Y = np.array([Xs[self.m0.keys[(I, J)]] for I, J in fp])
        Wend, Nend = e["Wend"], e["Nend"]; Oend = Wend + DELTA * Nend
        dx, dy = dir_vec(e["ray"][1])
        perp = np.array([1.0, 0.0]) if dx == 0 else np.array([0.0, 1.0])     # increasing side index
        wall_side = 0 if Nend @ perp > 0 else 1
        if "wall_side" in e:
            assert e["wall_side"] == wall_side
        e["wall_side"] = wall_side
        kw, ko = (k0, k1) if wall_side == 0 else (k1, k0)
        lam = polylen(Y); lam /= lam[-1]
        foot = Y[-1]
        gf = next(g for g, (li, ei, ids, tau) in enumerate(m.bnd) if kw[-1] in ids or ko[-1] in ids)
        a, b, D, li, ei = self.E[gf]; ed = self.lay.absd[li]["edges"][ei]
        Pf = np.asarray(loops[li][ed["seg"]]["pts"]); s_f, _ = closest(Pf, foot); cP = polylen(Pf)
        tf = at_s(Pf, np.array([min(s_f + 1e-3, cP[-1])]))[0] - at_s(Pf, np.array([max(s_f - 1e-3, 0)]))[0]
        sgn = 1 if tf @ Nend > 0 else -1
        foot_o = at_s(Pf, np.array([s_f + sgn * DELTA]))[0]
        # boundary nodes of the foot edge on the outer side: compress into [s_f + DELTA, end] (before setting L/R)
        ids = np.array(m.bnd[gf][2]); sc = np.array([closest(Pf, p)[0] for p in X[ids]])
        side = ((sc - s_f) * sgn > 1e-9) & ~np.isin(ids, np.r_[kw, ko])
        if sgn > 0:
            sc_new = s_f + DELTA + (sc - s_f) * (cP[-1] - s_f - DELTA) / (cP[-1] - s_f)
        else:
            sc_new = s_f - DELTA - (s_f - sc) * (s_f - DELTA) / s_f
        X[ids[side]] = at_s(Pf, sc_new[side])
        L = Y + (1 - lam)[:, None] * (Wend - Y[0])
        R = Y + (1 - lam)[:, None] * (Oend - Y[0]) + lam[:, None] * (foot_o - foot)
        X[kw] = L; X[ko] = R
        S = L[:, None, :] + (self.dk / DELTA)[None, :, None] * (R - L)[:, None, :]
        return dict(S=S, kw=kw, ko=ko)

    # ------------------------------------------------------------ assembly
    def assemble(self, X, LAY, STRIP):
        """positions of all nodes + quads (4 ids, block id) + structured 2D blocks (ids arrays)"""
        m = self.m
        P = np.zeros((self.NT, 2)); P[:m.N] = X
        quads = [m.cells[:, :5].copy()]
        blocks = []                                   # (kind, ids array (a, b)) ; outside blocks handled separately
        nb = self.nb2
        for ci, c in enumerate(self.chains):
            Lid = c["lay_ids"].copy()
            if not c["closed"] and c["end"]["kind"] == "strip":
                Lid[-1, 0] = STRIP[ci, "end"]["kw"][0]       # wall node at the strip end = outside node
            if not c["closed"] and c["start"]["kind"] == "strip":
                Lid[0, 0] = STRIP[ci, "start"]["kw"][0]
            c["Lid"] = Lid
            lay = LAY[ci]
            P[Lid[:, :-1].ravel()] = lay[:, :-1].reshape(-1, 2)
            n = len(Lid); adj = c["adj_block"]; ne = len(adj)
            # split points: where the adjacent outside block changes (+ chain ends)
            cuts = [i for i in range(1, ne) if adj[i] != adj[i - 1]]
            if c["closed"]:
                if adj[0] != adj[-1]:
                    cuts = [0] + cuts
                segs = [(cuts[k], cuts[(k + 1) % len(cuts)]) for k in range(len(cuts))]
            else:
                cuts = [0] + cuts + [ne]
                segs = list(zip(cuts[:-1], cuts[1:]))
            for a, b in segs:
                rng = list(range(a, b + 1)) if b > a else list(range(a, ne)) + list(range(0, b + 1))
                ids = Lid[[i % n for i in rng]]
                blocks.append(("layer", ids, ci)); nb += 1
        for (ci, end), S in STRIP.items():
            c = self.chains[ci]; e = c[end]
            Sid = e["strip_ids"].copy(); Sid[:, 0] = S["kw"]; Sid[:, -1] = S["ko"]
            Sid[0] = c["Lid"][-1 if end == "end" else 0]
            P[Sid[1:, 1:-1].ravel()] = S["S"][1:, 1:-1].reshape(-1, 2)
            adj = e["adj_block"]; ne = len(adj)
            cuts = [0] + [j for j in range(1, ne) if adj[j] != adj[j - 1]] + [ne]
            for a, b in zip(cuts[:-1], cuts[1:]):
                blocks.append(("strip", Sid[a:b + 1], ci))
        # orient every extra block CCW and emit quads
        out_blocks = []
        bid = self.nb2
        if not hasattr(self, "_flip"):
            self._flip = {}
        for bi_, (kind, ids, ci) in enumerate(blocks):
            if bi_ not in self._flip:                 # orientation fixed once (summed signed cell area)
                Q = P[ids]; a, b, c_, d = Q[:-1, :-1], Q[1:, :-1], Q[1:, 1:], Q[:-1, 1:]
                A = 0.5 * ((c_ - a)[..., 0] * (d - b)[..., 1] - (c_ - a)[..., 1] * (d - b)[..., 0]).sum()
                self._flip[bi_] = A < 0
            if self._flip[bi_]:
                ids = ids[:, ::-1]
            out_blocks.append(dict(kind=kind, ids=ids, chain=ci, bid=bid))
            q = np.stack([ids[:-1, :-1], ids[1:, :-1], ids[1:, 1:], ids[:-1, 1:]], -1).reshape(-1, 4)
            quads.append(np.c_[q, np.full(len(q), bid)]); bid += 1
        Qa = np.vstack(quads)
        if UNTANGLE:
            free = np.zeros(self.NT, bool); free[self.m.inner] = True
            P = untangle(P, Qa, free)
        return P, Qa, out_blocks


BLUNT = 0.0          # knife edges: wall cut where the local thickness < BLUNT * nominal thickness


def blunt(W, closed, t):
    """remove knife-edge slivers: wall points whose distance to the nearest opposite-facing wall point is below
    BLUNT*t form runs around an apex; each run is replaced by a straight cut. Returns the polyline and the
    removed length (max over runs)"""
    from scipy.spatial import cKDTree
    P = W[:-1] if closed and np.linalg.norm(W[0] - W[-1]) < 1e-9 else W
    T = np.gradient(P, axis=0); T /= np.linalg.norm(T, axis=1)[:, None] + 1e-300
    N = np.c_[T[:, 1], -T[:, 0]]
    tree = cKDTree(P); thin = np.zeros(len(P), bool)
    for i, nb in enumerate(tree.query_ball_point(P, BLUNT * t)):
        nb = np.array(nb)
        if len(nb) and np.any((N[nb] @ N[i] < -0.7) & (np.abs(nb - i) > 3)):
            thin[i] = True
    if not thin.any():
        return W, 0.0
    keep = ~thin
    Q = P[keep]
    removed = np.sum(np.linalg.norm(np.diff(P, axis=0), axis=1)[thin[:-1]])
    if closed:
        Q = np.vstack([Q, Q[:1]])
    return Q, float(removed)


def corner_between(W, s0, s1, t, closed):
    """arc position of the sharpest wall corner (turn > 45 deg over +-0.15 t) between s0 and s1, else None"""
    Lw = polylen(W)[-1]; w = 0.15 * t
    ss = np.linspace(s0, s1, 200)
    f = (lambda x: x % Lw) if closed else (lambda x: np.clip(x, 0, Lw))
    a = at_s(W, f(ss - w)); b = at_s(W, f(ss)); c = at_s(W, f(ss + w))
    t1 = b - a; t2 = c - b
    ang = np.degrees(np.abs(np.arctan2(t1[:, 0] * t2[:, 1] - t1[:, 1] * t2[:, 0], (t1 * t2).sum(1))))
    k = np.argmax(ang)
    return ss[k] if ang[k] > 45 else None


def quad_sj(P, q):
    Q = P[q[:, :4]]; sj = np.full(len(q), np.inf)
    for c in range(4):
        p0, p1, p3 = Q[:, c], Q[:, (c + 1) % 4], Q[:, (c + 3) % 4]
        e1, e2 = p1 - p0, p3 - p0
        sj = np.minimum(sj, (e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0]) / (np.linalg.norm(e1, axis=1) * np.linalg.norm(e2, axis=1) + 1e-300))
    return sj


def winslow(m, X, iters=8000, tol=1e-9):
    inn, nb = m.inner, m.nb
    scale = np.ptp(X[m.is_b], axis=0).max()
    for it in range(iters):
        E_, W_, N_, S_, NE, NW, SE, SW = (X[nb[:, c]] for c in range(8))
        xi = 0.5 * (E_ - W_); eta = 0.5 * (N_ - S_)
        a = (eta ** 2).sum(1)[:, None]; g = (xi ** 2).sum(1)[:, None]; b = (xi * eta).sum(1)[:, None]
        new = (a * (E_ + W_) + g * (N_ + S_) - 0.5 * b * (NE - NW - SE + SW)) / (2 * (a + g) + 1e-300)
        d = np.abs(new - X[inn]).max(); X[inn] = new
        if d < tol * scale:
            break
    return X, it + 1


def untangle(P, Q, free, thr=0.1, rounds=6):
    """local optimisation: for free nodes of cells with sJ < thr, move the node to maximise the minimum corner sJ
    of its adjacent cells (Nelder-Mead, 2 unknowns per node); repeated a few rounds"""
    from scipy.optimize import minimize
    n2c = collections.defaultdict(list)
    for ci, q in enumerate(Q[:, :4]):
        for v in q:
            n2c[v].append(ci)
    for _ in range(rounds):
        sj = quad_sj(P, Q); bad = np.where(sj < thr)[0]
        if not len(bad):
            break
        nodes = {v for c in bad for v in Q[c, :4] if free[v]}
        nodes |= {w for v in list(nodes) for c in n2c[v] for w in Q[c, :4] if free[w]}
        for v in nodes:
            cs = np.array(n2c[v]); x0 = P[v].copy()
            h = np.median([np.linalg.norm(P[a] - P[b]) for c in cs for a, b in zip(Q[c, :4], np.roll(Q[c, :4], 1))])
            def f(x):
                P[v] = x; return -quad_sj(P, Q[cs]).min()
            r = minimize(f, x0, method="Nelder-Mead", options=dict(xatol=1e-4 * h, fatol=1e-6, initial_simplex=[x0, x0 + [0.3 * h, 0], x0 + [0, 0.3 * h]], maxiter=200))
            P[v] = r.x if r.fun <= f(x0) else x0
    return P
