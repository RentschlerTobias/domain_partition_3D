"""Coarse quad layouts of 2D domains (rectilinear abstraction + separatrices), the building block for 2D
blocking: section layouts of sweeps, the hub surface of a runner, profile sections.

    dom = Domain2D([Loop(outer), Loop(hole, pins=[...])], period=None)   # domain on the left of every loop
    lay = layout(dom)                     # -> Layout2D: index-space blocks, vertices, rays
    lay = collapse(dom, h=0.05, guard=0.2) # the same, after a guarded greedy chord collapse
    mesh = lay.mesh(dom, h=0.05)          # fine structured mesh per block, Winslow re-fit to the boundary

Method (generic, as in the r-z section layout of the sweep strategy, layout2d.py there):
  1. Every loop is split into parameter edges with a direction in {0, 90, 180, 270} deg: at pins (points that
     must become layout vertices, e.g. block corners of an attached structure), at sharp corners and where the
     tangent turns more than 45 deg (+ hysteresis) away from the current direction (virtual corners).
     Valence of a vertex v = 2 - turn/90 (turn > 0 = left, convex for the domain on the left).
  2. Parameter coordinates: the two ends of a horizontal edge share a y group, of a vertical edge an x group
     (union-find); a group's coordinate is the mean physical coordinate of its vertices.
  3. Separatrices: every vertex with v >= 2 sends v-1 rays into the domain until they hit the boundary; the
     faces of this arrangement in the index grid of the group coordinates are the blocks (T-junction free).
  4. Chord collapse: merging two neighbouring groups removes a chord of blocks; accepted greedily if the
     parameter polygon stays valid, the block count drops and the re-fitted fine mesh keeps min sJ >= guard.
A loop without any vertex (e.g. the straight end line of a periodic band) is a free loop: one group, no rays.
Periodic domains (``period`` in x): x is unwrapped, the index space wraps, rays and blocks may cross the seam.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


# ------------------------------------------------------------------ input
@dataclass
class Loop:
    points: np.ndarray                       # (n, 2) closed polyline (last point != first), domain on the left
    pins: tuple = ()                         # point indices that must be layout vertices
    shift: float = 0.0                       # x added to point 0 when closing the loop (periodic band lines)

    def __post_init__(self):
        self.points = np.asarray(self.points, float).reshape(-1, 2)
        self.pins = tuple(sorted(int(p) for p in self.pins))

    def closed(self):
        """Points with the closing point appended (shifted for a band line)."""
        return np.vstack([self.points, self.points[:1] + [self.shift, 0.0]])


@dataclass
class Domain2D:
    loops: list
    period: float | None = None              # x period (a cylinder unrolled), None = plane


class LayoutError(RuntimeError):
    pass


def best_shear(dom: Domain2D, s_max: float = 3.0, n: int = 241) -> float:
    """Shear s of x' = x - s*y that makes the boundary most axis-parallel (length-weighted squared deviation of
    the segment directions from the nearest axis). A shear keeps x-periodicity, unlike a rotation; for a
    blade footprint on a hub it turns the stagger direction vertical."""
    segs = [np.diff(lp.closed(), axis=0) for lp in dom.loops]
    d = np.vstack(segs)
    L = np.linalg.norm(d, axis=1)
    best, arg = np.inf, 0.0
    for sv in np.linspace(-s_max, s_max, n):
        dd = d.copy()
        dd[:, 0] -= sv * dd[:, 1]
        a = np.arctan2(dd[:, 1], dd[:, 0])
        dev = (a + np.pi / 4) % (np.pi / 2) - np.pi / 4
        c = float(np.sum(np.linalg.norm(dd, axis=1) * dev ** 2))
        if c < best - 1e-12:
            best, arg = c, float(sv)
    return arg


def sheared(dom: Domain2D, s: float) -> Domain2D:
    """The domain in sheared coordinates x' = x - s*y (period and band shifts unchanged)."""
    return Domain2D([Loop(lp.points - np.c_[s * lp.points[:, 1], np.zeros(len(lp.points))], lp.pins, lp.shift)
                     for lp in dom.loops], dom.period)


# ------------------------------------------------------------------ 1. rectilinear abstraction
def _angles(P):
    d = np.diff(P, axis=0)
    return np.degrees(np.arctan2(d[:, 1], d[:, 0]))


def _wrap180(a):
    return (a + 180.0) % 360.0 - 180.0


def abstract(loop: Loop, hyst=10.0, sharp=60.0, snap=0.04):
    """Parameter edges of a loop: list of dict(i0, i1 (point indices on the closed polyline, i1 may be n),
    D (deg, multiple of 90, unwrapped), v (valence of the start vertex)). Empty list = free loop.
    A virtual corner within ``snap`` x loop length of a pin is moved onto the pin (no tiny parameter edges)."""
    P = loop.closed()
    n = len(loop.points)
    A = _angles(P)                                  # segment k: point k -> k+1
    turn = np.r_[_wrap180(A[0] - A[-1]), _wrap180(np.diff(A))]   # turn at point k (into segment k)
    pins = set(loop.pins) | {k for k in range(n) if abs(turn[k]) > sharp}
    artificial = not pins
    if not pins:
        # no fixed vertex: a free loop if it never turns away from one direction, else start a virtual one
        D0 = 90.0 * round(A[0] / 90.0)
        dev = np.cumsum(np.r_[0.0, turn[1:]]) + A[0] - D0
        if np.all(np.abs(dev) <= 45 + hyst):
            return []
        pins = {0}
    start = min(pins)
    order = [(start + k) % n for k in range(n)]
    Acont = np.empty(n)                             # unwrapped tangent angle of the segments, from start
    Acont[0] = A[start]
    for m in range(1, n):
        Acont[m] = Acont[m - 1] + turn[order[m]]
    seglen = np.linalg.norm(np.diff(P, axis=0), axis=1)[order]
    arc = np.r_[0.0, np.cumsum(seglen)]             # arc length at walk step m (point order[m])
    reach = snap * arc[-1]
    pin_m = sorted(m for m in range(n) if order[m] in pins) + [n]
    next_pin = {}
    for m in range(n):
        next_pin[m] = next(q for q in pin_m if q > m)

    def turned(D, m):
        """Direction after the turns the tangent at step m asks for, starting from D."""
        dev = Acont[m] - D
        while abs(dev) > 45 + hyst:
            D += 90 * np.sign(dev)
            dev = Acont[m] - D
        return D
    edges = []
    D = 90.0 * round(Acont[0] / 90.0)
    m2 = 0                                          # look ahead from the start pin
    while m2 + 1 < n and arc[m2 + 1] - arc[0] <= reach and m2 + 1 < next_pin[0]:
        m2 += 1
    D = turned(D, m2)
    edges.append(dict(i0=start, D=D))
    for m in range(1, n):
        k = order[m]
        if k in pins:
            Dn = 90.0 * round(Acont[m] / 90.0)
            m2 = m                                  # take a turn that follows shortly after at the pin
            while m2 + 1 < n and arc[m2 + 1] - arc[m] <= reach and m2 + 1 < next_pin[m]:
                m2 += 1
            Dn = turned(Dn, m2)
            edges[-1]["i1"] = k
            edges.append(dict(i0=k, D=Dn, turn_in=Dn - D))
            D = Dn
            continue
        if abs(Acont[m] - D) > 45 + hyst and arc[next_pin[m]] - arc[m] <= reach and next_pin[m] < n:
            continue                                # the next pin is close: it takes this turn
        dev = Acont[m] - D
        while abs(dev) > 45 + hyst:                 # virtual corner(s) at point k
            sgn = np.sign(dev)
            Dn = D + 90 * sgn
            edges[-1]["i1"] = k
            edges.append(dict(i0=k, D=Dn, turn_in=Dn - D))
            D = Dn
            dev = Acont[m] - D
    edges[-1]["i1"] = start + n if start == 0 else start
    # closing turn into the first edge
    total = Acont[-1] + turn[start] - Acont[0]      # total turning of the polyline
    first_in = D                                    # direction arriving at the start vertex
    D0 = edges[0]["D"]
    # the start direction seen after one turn around the loop: D0 + round(total/90)*90 (360 for a ccw loop)
    D0_after = D0 + 90.0 * round(total / 90.0)
    edges[0]["turn_in"] = D0_after - first_in
    for e in edges:                                 # point range i0..i1 along the loop, i1 > i0 (may pass n)
        if e["i1"] <= e["i0"]:
            e["i1"] += n
    if artificial and len(edges) > 1 and abs(edges[0]["turn_in"]) < 1e-9:
        # the start point was only a seed for the walk: no turn there, so it is not a vertex
        last = edges.pop()
        first = edges[0]
        edges[0] = dict(last, i1=last["i1"] + (first["i1"] - first["i0"]))
    for e in edges:
        e["v"] = int(round(2 - e["turn_in"] / 90.0))
        if e["v"] < 1:
            raise LayoutError(f"vertex at point {e['i0']}: turn {e['turn_in']} deg (a cusp; give it as a pin)")
    for e in edges:
        if len(edges) > 1 and (e["i1"] - e["i0"]) % n == 0:
            raise LayoutError(f"empty parameter edge at point {e['i0']}")
    return edges


# ------------------------------------------------------------------ 2.-3. layout in index space
class _UF:
    def __init__(self, n):
        self.p = list(range(n))

    def f(self, x):
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def u(self, a, b):
        ra, rb = self.f(a), self.f(b)
        if ra != rb:
            self.p[max(ra, rb)] = min(ra, rb)


def _dirvec(D):
    return [(1, 0), (0, 1), (-1, 0), (0, -1)][int(round(D / 90.0)) % 4]


@dataclass
class Layout2D:
    """Blocks in index space: ``blocks`` (i0, i1, j0, j1) with i unwrapped (i1 may exceed nx for a block across
    the seam of a periodic domain); ``X``/``Y`` the group coordinates of the index lines."""
    X: np.ndarray
    Y: np.ndarray
    blocks: list
    verts: list                              # (loop, edge) per layout vertex of the boundary
    vij: list                                # index position per vertex
    edges: list                              # (a, b, D, loop, edge) parameter edges
    rays: list
    cells: set
    period: float | None
    merges: tuple = ()
    absd: list = field(default_factory=list)

    @property
    def n_blocks(self):
        return len(self.blocks)


def _vertices(dom: Domain2D, absd):
    verts, vpos, edges = [], [], []
    for li, (lp, ed) in enumerate(zip(dom.loops, absd)):
        P = lp.closed()
        base = len(verts)
        for ei, e in enumerate(ed):
            verts.append((li, ei))
            vpos.append(P[e["i0"] % len(lp.points)] + ([lp.shift, 0] if e["i0"] >= len(lp.points) else 0))
        for ei, e in enumerate(ed):
            edges.append((base + ei, base + (ei + 1) % len(ed), e["D"], li, ei))
    return verts, np.array(vpos).reshape(-1, 2), edges


def _facing_pairs(dom, absd, verts, edges, ax):
    """Pairs of parallel parameter edges across axis ``ax`` whose physical extents overlap along the other axis:
    (a_vertex, b_vertex, shift, gap) with the edge of a physically before the edge of b (gap > 0; shift = period
    multiple added to b's coordinate). Their parameter order must stay the physical order."""
    per = dom.period if ax == 0 else None
    E = []
    for a, b, D, li, ei in edges:
        if _dirvec(D)[ax] != 0:
            continue                              # an edge along ax does not face across ax
        e = absd[li][ei]
        Q = _piece(dom.loops[li], e["i0"], e["i1"])
        E.append((a, float(Q[:, ax].mean()), float(Q[:, 1 - ax].min()), float(Q[:, 1 - ax].max())))
    out = []
    for i in range(len(E)):
        for j in range(len(E)):
            if i == j:
                continue
            a, xa, lo_a, hi_a = E[i]
            b, xb, lo_b, hi_b = E[j]
            for w in ((-1, 0, 1) if per else (0,)):
                gap = xb + w * (per or 0) - xa
                if gap <= 0:
                    continue
                if per and gap >= per:
                    continue
                if min(hi_a, hi_b) - max(lo_a, lo_b) > 0:
                    out.append((a, b, w * (per or 0), gap))
    return out


def _groups(dom, verts, vpos, edges, free, merges, absd=None):
    n = len(verts)
    ux, uy = _UF(n + len(free)), _UF(n + len(free))
    for a, b, D, li, ei in edges:
        dx, dy = _dirvec(D)
        (uy if dy == 0 else ux).u(a, b)
    for ax, a, b in merges:
        (ux if ax == 0 else uy).u(a, b)
    U = np.zeros((n, 2))
    for ax, uf in ((0, ux), (1, uy)):
        g = {}
        for i in range(n):
            g.setdefault(uf.f(i), []).append(i)
        roots = list(g)
        gi = {r: k for k, r in enumerate(roots)}
        mean = np.array([vpos[g[r], ax].mean() for r in roots])
        # order constraints: an edge along this axis has positive length in its direction
        cons = []
        for a, b, D, li, ei in edges:
            d = _dirvec(D)[ax]
            if d == 0:
                continue
            ga, gb = gi[uf.f(a)], gi[uf.f(b)]
            span = abs(vpos[b] - vpos[a]).sum()
            cons.append((ga, gb, d, 0.05 * span + 1e-9))
        if absd is not None:
            for a, b, w, gap in _facing_pairs(dom, absd, verts, edges, ax):
                ga, gb = gi[uf.f(a)], gi[uf.f(b)]
                if ga != gb:
                    cons.append((ga, gb, 1, 0.25 * gap - w))     # x[gb] + w - x[ga] >= gap/4
        val = mean.copy()
        if any(d * (val[gb] - val[ga]) < m - 1e-12 for ga, gb, d, m in cons):
            val = _project_order(mean, cons)
        for r, k in gi.items():
            U[g[r], ax] = val[k]
    return np.round(U, 10) + 0.0, ux, uy


def _project_order(mean, cons):
    """Closest coordinates (least squares) to ``mean`` with d*(x[b]-x[a]) >= m for every (a, b, d, m)."""
    from scipy.optimize import minimize

    A = np.zeros((len(cons), len(mean)))
    lb = np.zeros(len(cons))
    for r, (ga, gb, d, m) in enumerate(cons):
        A[r, gb] += d
        A[r, ga] -= d
        lb[r] = m
    res = minimize(lambda x: 0.5 * np.sum((x - mean) ** 2), mean, jac=lambda x: x - mean, method="SLSQP",
                   constraints=[{"type": "ineq", "fun": lambda x: A @ x - lb, "jac": lambda x: A}],
                   options={"maxiter": 500, "ftol": 1e-14})
    if not res.success or np.any(A @ res.x - lb < -1e-9):
        raise LayoutError(f"no parameter coordinates with consistent edge directions ({res.message})")
    return res.x


def layout(dom: Domain2D, merges=(), absd=None) -> Layout2D:
    absd = absd if absd is not None else [abstract(lp) for lp in dom.loops]
    free = [li for li, a in enumerate(absd) if not a]
    verts, vpos, edges = _vertices(dom, absd)
    if not verts:
        raise LayoutError("no layout vertex at all (a domain without corners needs a pin)")
    U, ux, uy = _groups(dom, verts, vpos, edges, free, merges, absd)
    per = dom.period
    # parameter edges must have positive length along their direction and no offset across it
    for a, b, D, li, ei in edges:
        dx, dy = _dirvec(D)
        d = U[b] - U[a]
        if per and dy == 0 and dx != 0:
            d[0] = (d[0] * dx) % per * dx if (d[0] * dx) % per > 1e-12 else d[0]
        L = d @ np.array([dx, dy])
        off = d @ np.array([-dy, dx])
        if L <= 1e-9 or abs(off) > 1e-9:
            raise LayoutError(f"parameter edge {a}->{b} (loop {li}) has length {L:.4g}, offset {off:.2g}")
    # index lines
    xs = np.mod(U[:, 0], per) if per else U[:, 0]
    X = np.unique(np.round(xs, 9))
    ys = [U[:, 1]] + [np.array([dom.loops[li].points[:, 1].mean()]) for li in free]
    Y = np.unique(np.round(np.concatenate(ys), 9))
    nx = len(X)
    ix = {x: k for k, x in enumerate(X)}
    iy = {y: k for k, y in enumerate(Y)}
    vij = [(ix[round(float(x), 9)], iy[round(float(u[1]), 9)]) for x, u in zip(xs, U)]
    fy = {li: iy[round(float(dom.loops[li].points[:, 1].mean()), 9)] for li in free}
    ncx = nx if per else nx - 1

    def wrap(i):
        return i % nx if per else i

    def cx(i):                                   # x extent of index cell column i
        i = wrap(i)
        lo = X[i]
        hi = X[i + 1] if i + 1 < nx else X[0] + per
        return lo, hi

    # boundary edges (index-grid unit edges) of the parameter loops
    bedges, bnodes = set(), set()

    def ek(p, dx, dy):
        """Unit edge leaving node p in direction (dx, dy): ('h', column, row) or ('v', column line, row)."""
        if dy == 0:
            return ("h", p[0] if dx > 0 else wrap(p[0] - 1), p[1])
        return ("v", wrap(p[0]), p[1] if dy > 0 else p[1] - 1)
    for a, b, D, li, ei in edges:
        dx, dy = _dirvec(D)
        p = vij[a]
        steps = 0
        while p != vij[b] or steps == 0:
            q = (wrap(p[0] + dx), p[1] + dy)
            bedges.add(ek(p, dx, dy))
            bnodes.add(p)
            p = q
            steps += 1
            if steps > 10 * (nx + len(Y)):
                raise LayoutError(f"parameter edge {a}->{b} does not reach its end vertex")
        bnodes.add(p)
    for li, j in fy.items():                        # free loops: a full row of boundary edges
        for i in range(ncx):
            bedges.add(("h", i, j))
            bnodes.add((wrap(i), j))
    # domain cells: parity of crossings of an upward ray from the cell centre with the parameter loops
    segs = []
    for a, b, D, li, ei in edges:
        dx, dy = _dirvec(D)
        if dx != 0:                                  # horizontal parameter edge
            x0 = U[a, 0]
            L = (U[b, 0] - U[a, 0]) * dx
            if per:
                L = L % per if L % per > 1e-12 else L
            segs.append((min(x0, x0 + dx * L), max(x0, x0 + dx * L), U[a, 1]))
    for li in free:
        segs.append((-np.inf, np.inf, Y[fy[li]]))
    cells = set()
    for i in range(ncx):
        lo, hi = cx(i)
        xc = 0.5 * (lo + hi)
        for j in range(len(Y) - 1):
            yc = 0.5 * (Y[j] + Y[j + 1])
            cnt = 0
            for x0, x1, y in segs:
                if y <= yc:
                    continue
                if per and np.isfinite(x0):
                    k = np.floor((xc - x0) / per)
                    xx = xc - k * per                # bring xc into [x0, x0 + per)
                    hit = x0 < xx < x1
                else:
                    hit = x0 < xc < x1
                cnt += hit
            if cnt % 2:
                cells.add((i, j))

    def inside_edge(p, dx, dy):
        if dy == 0:
            c = p[0] if dx > 0 else wrap(p[0] - 1)
            return (c, p[1] - 1) in cells and (c, p[1]) in cells
        j = p[1] if dy > 0 else p[1] - 1
        return (wrap(p[0] - 1), j) in cells and (wrap(p[0]), j) in cells
    rays, redges = [], set()
    for vi, (li, ei) in enumerate(verts):
        e = absd[li][ei]
        for m in range(1, e["v"]):
            dx, dy = _dirvec(e["D"] + 90 * m)
            p, path = vij[vi], [vij[vi]]
            for _ in range(10 * (nx + len(Y))):
                q = (wrap(p[0] + dx), p[1] + dy)
                if ek(p, dx, dy) in bedges or not inside_edge(p, dx, dy):
                    raise LayoutError(f"ray from vertex {vi} leaves the domain at {p}->{q}")
                redges.add(ek(p, dx, dy))
                path.append(q)
                p = q
                if q in bnodes:
                    break
            else:
                raise LayoutError(f"ray from vertex {vi} does not end")
            rays.append((vi, e["D"] + 90 * m, path))
    walls = bedges | redges
    comp, blocks = {}, []
    for c in sorted(cells):
        if c in comp:
            continue
        b = len(blocks)
        comp[c] = (b, 0)
        stack, mem = [(c, 0)], [(c[0], c[1])]
        while stack:
            (i, j), o = stack.pop()
            for di, dj, wall in ((1, 0, ("v", wrap(i + 1), j)), (-1, 0, ("v", i, j)),
                                 (0, 1, ("h", i, j + 1)), (0, -1, ("h", i, j))):
                ni, nj = i + di, j + dj
                no = o + (1 if per and di == 1 and ni >= nx else -1 if per and di == -1 and ni < 0 else 0)
                n2 = (wrap(ni), nj)
                if n2 not in cells or wall in walls:
                    continue
                if n2 in comp:
                    if comp[n2] != (b, no):
                        raise LayoutError("a block wraps around the periodic domain (needs a vertex on the seam)")
                    continue
                comp[n2] = (b, no)
                stack.append((n2, no))
                mem.append((n2[0] + no * nx, n2[1]))
        I = [m[0] for m in mem]
        J = [m[1] for m in mem]
        i0, i1, j0, j1 = min(I), max(I) + 1, min(J), max(J) + 1
        if (i1 - i0) * (j1 - j0) != len(mem):
            raise LayoutError(f"block {b} is not a rectangle in index space")
        blocks.append((i0, i1, j0, j1))
    return Layout2D(X=X, Y=Y, blocks=blocks, verts=verts, vij=vij, edges=edges, rays=rays, cells=cells,
                    period=per, merges=tuple(merges), absd=absd)


# ------------------------------------------------------------------ 4. fine mesh + Winslow re-fit
def _piece(loop: Loop, i0, i1):
    """Polyline of the loop from point index i0 to i1 (i1 > i0, may pass n: then the loop's shift applies)."""
    n = len(loop.points)
    idx = np.arange(i0, i1 + 1)
    return loop.points[idx % n] + np.c_[(idx // n) * loop.shift, np.zeros(len(idx))]


def _at_fraction(P, t):
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
    s = s / s[-1]
    return np.c_[np.interp(t, s, P[:, 0]), np.interp(t, s, P[:, 1])]


@dataclass
class FineMesh:
    """Structured fine mesh of a layout: node positions ``X`` (unwrapped x), per block a node-id grid
    ``grids[b]`` of shape (ni+1, nj+1), cells (n, 4) counter-clockwise, quality = min scaled Jacobian."""
    X: np.ndarray
    grids: list
    cells: np.ndarray
    min_sj: float
    sj: np.ndarray


def mesh(lay: Layout2D, dom: Domain2D, h: float, nmin: int = 1, iters: int = 3000, tol: float = 1e-9) -> FineMesh:
    per = dom.period
    X, Y = lay.X, lay.Y
    nx = len(X)
    ncx = nx if per else nx - 1
    widths = [(X[c + 1] if c + 1 < nx else X[0] + per) - X[c] for c in range(ncx)]
    ncol = [max(nmin, int(round(w / h))) for w in widths]
    nrow = [max(nmin, int(round((Y[r + 1] - Y[r]) / h))) for r in range(len(Y) - 1)]
    oi, oj = np.r_[0, np.cumsum(ncol)], np.r_[0, np.cumsum(nrow)]
    NI = int(oi[-1])
    ux = np.concatenate([np.linspace(X[c], X[c] + widths[c], ncol[c] + 1)[:-1] for c in range(ncx)]
                        + ([] if per else [[X[-1]]]))
    uy = np.concatenate([np.linspace(Y[r], Y[r + 1], nrow[r] + 1)[:-1] for r in range(len(Y) - 1)] + [[Y[-1]]])

    def wI(I):
        return I % NI if per else I
    keys = {}
    fine_cells = []
    for (ci, cj) in sorted(lay.cells):
        for I in range(oi[ci], oi[ci + 1]):
            for J in range(oj[cj], oj[cj + 1]):
                q = []
                for A, B in ((I, J), (I + 1, J), (I + 1, J + 1), (I, J + 1)):
                    q.append(keys.setdefault((wI(A), B), len(keys)))
                fine_cells.append(q)
    N = len(keys)
    kl = list(keys)
    U = np.array([[ux[I] if I < len(ux) else ux[0] + per, uy[J]] for I, J in kl])   # parameter positions
    Xf = U.copy()
    fixed = np.zeros(N, bool)

    def fine_path(path, dx, dy):
        out = []
        for (i0, j0), (i1, j1) in zip(path[:-1], path[1:]):
            I0 = int(oi[i0])
            if dy != 0:
                J0, J1 = int(oj[j0]), int(oj[j1])
                seg = [(I0, J) for J in range(J0, J1 + dy, dy)]
            else:
                c = i0 if dx > 0 else (i0 - 1) % ncx
                seg = [(I0 + dx * k, int(oj[j0])) for k in range(ncol[c] + 1)]
            out += seg if not out else seg[1:]
        return out
    # boundary nodes of the parameter edges: arc-length fraction = parameter fraction along the edge
    for a, b, D, li, ei in lay.edges:
        e = lay.absd[li][ei]
        loop = dom.loops[li]
        dx, dy = _dirvec(D)
        p = lay.vij[a]
        path = [p]
        while True:
            q = ((p[0] + dx) % nx if per else p[0] + dx, p[1] + dy)
            path.append(q)
            p = q
            if q == lay.vij[b]:
                break
        fp = fine_path(path, dx, dy)
        uu = np.array([[ux[I % NI] if per else ux[I], uy[J]] for I, J in fp])
        d = np.r_[0, np.cumsum(np.abs(np.diff(uu, axis=0)).sum(1) % (per if per else np.inf))]
        tau = d / d[-1]
        pos = _at_fraction(_piece(loop, e["i0"], e["i1"]), tau)
        for (I, J), xy in zip(fp, pos):
            k = keys[(wI(I), J)]
            if per:
                xy = xy.copy()
                xy[0] -= per * np.round((xy[0] - U[k, 0]) / per)
            Xf[k] = xy
            fixed[k] = True
    free = [li for li, a in enumerate(lay.absd) if not a]
    for li in free:                                  # free loops: x from the parameter, y on the line
        loop = dom.loops[li]
        yl = loop.points[:, 1].mean()
        J = int(oj[list(Y).index(min(Y, key=lambda y: abs(y - yl)))])
        Pl = loop.closed()
        for (I, JJ), k in keys.items():
            if JJ == J:
                xq = U[k, 0]
                if per:
                    xs = Pl[:, 0]
                    lo = xs.min()
                    xq = lo + (xq - lo) % per
                o = np.argsort(Pl[:, 0])
                Xf[k] = [U[k, 0], np.interp(xq, Pl[o, 0], Pl[o, 1])]
                fixed[k] = True
    # Winslow on the free nodes (8-neighbour stencil in index space, periodic shift in x)
    inner = np.flatnonzero(~fixed)
    offs = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, 1), (1, -1), (-1, -1)]
    nb = np.zeros((len(inner), 8), int)
    sh = np.zeros((len(inner), 8))
    for r, k in enumerate(inner):
        I, J = kl[k]
        for c, (di, dj) in enumerate(offs):
            I2 = I + di
            nb[r, c] = keys[(wI(I2), J + dj)]
            sh[r, c] = (per * (I2 // NI)) if per else 0.0
    # start: the parameter positions shifted to the boundary (TFI-like through a Laplace solve would be better)
    for it in range(iters):
        E, W, Nn, S, NE, NW, SE, SW = (Xf[nb[:, c]] + np.c_[sh[:, c], np.zeros(len(inner))] for c in range(8))
        xi, eta = 0.5 * (E - W), 0.5 * (Nn - S)
        a = (eta ** 2).sum(1)[:, None]
        g = (xi ** 2).sum(1)[:, None]
        bb = (xi * eta).sum(1)[:, None]
        new = (a * (E + W) + g * (Nn + S) - 0.5 * bb * (NE - NW - SE + SW)) / (2 * (a + g) + 1e-300)
        d = np.abs(new - Xf[inner]).max() if len(inner) else 0.0
        Xf[inner] = new
        if d < tol:
            break
    C = np.array(fine_cells)
    P = Xf[C]
    if per:                                          # cells across the seam: unwrap to the first corner
        P[..., 0] -= per * np.round((P[..., 0] - P[:, :1, 0]) / per)
    sj = np.full(len(C), np.inf)
    for c in range(4):
        e1, e2 = P[:, (c + 1) % 4] - P[:, c], P[:, (c + 3) % 4] - P[:, c]
        cr = e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0]
        sj = np.minimum(sj, cr / (np.linalg.norm(e1, axis=1) * np.linalg.norm(e2, axis=1) + 1e-300))
    grids = []
    for i0, i1, j0, j1 in lay.blocks:
        I0, I1 = oi[i0 % ncx] + (i0 // ncx) * NI if per else oi[i0], None
        Iend = (oi[i1 % ncx] + (i1 // ncx) * NI) if per else oi[i1]
        if per and i1 % ncx == 0 and i1 > 0:
            Iend = (i1 // ncx) * NI
        G = np.array([[keys[(wI(I), J)] for J in range(oj[j0], oj[j1] + 1)] for I in range(I0, Iend + 1)])
        grids.append(G)
    return FineMesh(X=Xf, grids=grids, cells=C, min_sj=float(sj.min()) if len(sj) else 1.0, sj=sj)
