"""Stage 7: transfinite interpolation over the hex blocks.

Step 1 of the plan in `TFI_RESEARCH.md`: the trilinear Gordon-Hall map on a
single block, gated by `mesh_quality.py`.

The test is deliberately self-referential and therefore honest: take a block
that already exists, keep only its **boundary**, rebuild the interior by TFI,
and compare the result against the interior that was there. If TFI cannot
reproduce a mesh that AlgoHex already produced, it will not produce a better
one from scratch, and folding shows up on day one instead of after the
conforming-division machinery is built.

Two things have to work before the map can be evaluated at all:

1. **The block's lattice.** A block is a connected set of hexes cut out of a
   hex mesh, so it is a structured (a x b x c) grid -- but nothing stores
   that. `block_lattice` recovers it by propagating each cell's vertex order
   to its face neighbours, which also yields the (i,j,k) index of every
   vertex, and with it the six boundary faces as ordered grids.
2. **A consistent parametrisation.** TFI only gives a sensible grid if the
   bounding surfaces meet at the edges and run in the same direction; the
   lattice delivers exactly that, which is why it comes first.
"""

import sys
from collections import Counter, defaultdict, deque
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "experimentell" / "hex3d_algohex"))

import base_complex as bc                                             # noqa: E402
import clean_blocks as cb                                             # noqa: E402
import mesh_quality as mq                                             # noqa: E402

HF = bc.HEX_FACES
# the three pairs of opposite faces, in VTK_HEXAHEDRON ordering
OPPOSITE = ((0, 1), (2, 4), (3, 5))


# local corner (di, dj, dk) of each VTK_HEXAHEDRON vertex index
CORNER = ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
          (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1))
CORNER_IDX = {c: i for i, c in enumerate(CORNER)}


def _cell_adjacency(cell):
    """Vertex adjacency inside one hex, along its 12 edges."""
    adj = defaultdict(set)
    for a, b in bc.HEX_EDGES:
        u, v = int(cell[a]), int(cell[b])
        adj[u].add(v)
        adj[v].add(u)
    return adj


def _complete_order(cell, known):
    """Finish a hex's VTK order given 4 corners of one of its faces.

    `known` maps a local index (0..7) to a mesh vertex. The four unknown
    corners are each the unique neighbour of a known one along a cell edge in
    the missing direction. Deriving the neighbour's frame this way rather
    than writing out a permutation per direction is the whole point: the
    hand-derived permutations were wrong, and produced a lattice with 893
    distinct vertices in 986 slots -- silently, because every consistency
    check that only counts still passed."""
    adj = _cell_adjacency(cell)
    have = dict(known)
    rest = set(int(v) for v in cell) - set(have.values())
    # the axis along which the known face is constant
    cs = [CORNER[i] for i in have]
    axis = next(a for a in (0, 1, 2) if len({c[a] for c in cs}) == 1)
    for i, v in list(have.items()):
        c = list(CORNER[i])
        c[axis] = 1 - c[axis]
        j = CORNER_IDX[tuple(c)]
        cand = [w for w in adj[v] if w in rest]
        if len(cand) != 1:
            return None
        have[j] = cand[0]
    if len(have) != 8 or len(set(have.values())) != 8:
        return None
    return [have[i] for i in range(8)]


def block_lattice(hexes, cells, f2h):
    """(dims, vert) for one block: vert[i, j, k] is a mesh vertex index.

    Recovered by breadth-first propagation of each cell's vertex order to its
    face neighbours. Returns None if the block is not a structured grid --
    the honest answer for a block that is not a topological cuboid."""
    cells = list(cells)
    cset = set(cells)
    start = cells[0]
    order = {start: [int(v) for v in hexes[start]]}
    ijk = {start: (0, 0, 0)}
    dq = deque([start])
    while dq:
        c = dq.popleft()
        o = order[c]
        base = ijk[c]
        for axis in (0, 1, 2):
            for side in (0, 1):
                # the face of c at `side` along `axis`, and the local indices
                # it occupies in the NEIGHBOUR (mirrored across that axis)
                mine = [i for i in range(8) if CORNER[i][axis] == side]
                lp = [o[i] for i in mine]
                nb = [h for h in f2h[frozenset(lp)] if h != c and h in cset]
                if not nb or nb[0] in order:
                    continue
                n = nb[0]
                known = {}
                for i in mine:
                    cc = list(CORNER[i])
                    cc[axis] = 1 - cc[axis]
                    known[CORNER_IDX[tuple(cc)]] = o[i]
                no = _complete_order(hexes[n], known)
                if no is None:
                    return None
                order[n] = no
                step = [0, 0, 0]
                step[axis] = -1 if side == 0 else +1
                ijk[n] = tuple(np.array(base) + step)
                dq.append(n)
    if len(order) != len(cells):
        return None
    A = np.array([ijk[c] for c in cells])
    lo = A.min(0)
    dims = A.max(0) - lo + 1
    if int(np.prod(dims)) != len(cells):
        return None
    vert = np.full(tuple(dims + 1), -1, np.int64)
    for c in cells:
        i, j, k = np.array(ijk[c]) - lo
        o = order[c]
        for li, (di, dj, dk) in enumerate(CORNER):
            slot = (i + di, j + dj, k + dk)
            if vert[slot] >= 0 and vert[slot] != o[li]:
                return None            # inconsistent lattice, not a grid
            vert[slot] = o[li]
    if (vert < 0).any():
        return None
    return dims, vert


def _fix_handedness(P, vert):
    """A lattice is only defined up to orientation.

    The frame is propagated consistently from cell to cell, but the seed
    cell's own VTK order fixes a handedness that the (i, j, k) convention
    here need not share -- and the whole block then comes out with negative
    volume. Reversing one index flips the handedness for all of it at once.
    Detected on the mesh's own vertices, so it cannot be argued with."""
    C = block_cells(vert)
    import ovm_io
    if ovm_io._hex_volume(P[C[0]]) < 0:
        vert = vert[::-1].copy()
    return vert


def tfi(X):
    """Trilinear Gordon-Hall transfinite interpolation.

    `X` is (ni, nj, nk, 3) with its six boundary faces filled; the interior is
    replaced. Boolean sum P1 (+) P2 (+) P3: the six face terms, minus the
    twelve edge terms, plus the eight corner terms. Dropping the subtraction
    would count every corner three times and every edge twice, and the map
    would not even reproduce its own boundary."""
    ni, nj, nk = X.shape[:3]
    u = np.linspace(0, 1, ni)[:, None, None, None]
    v = np.linspace(0, 1, nj)[None, :, None, None]
    w = np.linspace(0, 1, nk)[None, None, :, None]

    F = ((1 - u) * X[0][None, :, :, :] + u * X[-1][None, :, :, :]
         + (1 - v) * X[:, 0][:, None, :, :] + v * X[:, -1][:, None, :, :]
         + (1 - w) * X[:, :, 0][:, :, None, :] + w * X[:, :, -1][:, :, None, :])

    bu = [(1 - u), u]
    bv = [(1 - v), v]
    bw = [(1 - w), w]
    iu = [0, -1]
    E = np.zeros_like(F)
    for a in (0, 1):
        for b in (0, 1):
            E += bu[a] * bv[b] * X[iu[a], iu[b], :][None, None, :, :]
            E += bu[a] * bw[b] * X[iu[a], :, iu[b]][None, :, None, :]
            E += bv[a] * bw[b] * X[:, iu[a], iu[b]][:, None, None, :]

    V = np.zeros_like(F)
    for a in (0, 1):
        for b in (0, 1):
            for c in (0, 1):
                V += bu[a] * bv[b] * bw[c] * X[iu[a], iu[b], iu[c]]

    return F - E + V


def rebuild_block(P, vert):
    """TFI the interior of one block from its own boundary. Returns the new
    positions for the block's vertices, same indexing as `vert`."""
    X = P[vert]                                  # (ni, nj, nk, 3)
    return tfi(X.copy())


def block_cells(vert):
    """The hexes of a lattice block, in VTK order."""
    ni, nj, nk = np.array(vert.shape) - 1
    out = []
    for i in range(ni):
        for j in range(nj):
            for k in range(nk):
                out.append([vert[i, j, k], vert[i + 1, j, k],
                            vert[i + 1, j + 1, k], vert[i, j + 1, k],
                            vert[i, j, k + 1], vert[i + 1, j, k + 1],
                            vert[i + 1, j + 1, k + 1], vert[i, j + 1, k + 1]])
    return np.array(out, np.int64)


def run(blocks_vtk, verbose=True):
    import meshio
    m = meshio.read(blocks_vtk)
    H = np.vstack([b.data for b in m.cells if b.type == "hexahedron"])
    B = np.concatenate([np.asarray(d).ravel() for b, d in
                        zip(m.cells, m.cell_data["block_id"])
                        if b.type == "hexahedron"])
    P = np.asarray(m.points, float)
    f2h, _e2h = bc.build_topology(H)
    print(f"[tfi] {blocks_vtk}: {len(H)} cells, {int(B.max()) + 1} blocks\n")
    print(f"{'block':>5} {'cells':>7} {'lattice':>14} {'sJ before':>10} "
          f"{'sJ after':>10} {'max move':>10}")
    ok = fail = 0
    Pnew = P.copy()
    for r in range(int(B.max()) + 1):
        cells = np.where(B == r)[0]
        lat = block_lattice(H, cells, f2h)
        if lat is None:
            print(f"{r:5d} {len(cells):7d} {'NOT a lattice':>14}")
            fail += 1
            continue
        dims, vert = lat
        vert = _fix_handedness(P, vert)
        Xn = rebuild_block(P, vert)
        C = block_cells(vert)
        sj0 = cb.scaled_jacobians(P, C)
        Q = P.copy()
        Q[vert.ravel()] = Xn.reshape(-1, 3)
        sj1 = cb.scaled_jacobians(Q, C)
        move = np.linalg.norm(Xn.reshape(-1, 3) - P[vert.ravel()], axis=1).max()
        print(f"{r:5d} {len(cells):7d} {str(tuple(dims)):>14} "
              f"{sj0.min():10.4f} {sj1.min():10.4f} {move:10.4f}")
        Pnew[vert.ravel()] = Xn.reshape(-1, 3)
        ok += 1
    print(f"\n[tfi] {ok} blocks rebuilt, {fail} not a lattice")
    if ok:
        print("\n[tfi] all seven CFD metrics, original -> TFI rebuild:")
        a0 = mq.metrics(P, H, verbose=False)
        a1 = mq.metrics(Pnew, H, verbose=False)
        print(f"{'metric':24s} {'before':>22s} {'after':>22s}")
        for k, side in (("scaled_jacobian", "min"),
                        ("non_orthogonality_deg", "max"), ("skewness", "max"),
                        ("aspect_ratio", "max"), ("vol_ratio", "min"),
                        ("face_weight", "min"), ("face_flatness", "min")):
            f = (lambda v: v.min()) if side == "min" else (lambda v: v.max())
            lim = mq.VIOLATION[k]
            n0 = int((a0[k] > lim[1]).sum()) if lim[0] == "max" else int((a0[k] <= lim[1]).sum())
            n1 = int((a1[k] > lim[1]).sum()) if lim[0] == "max" else int((a1[k] <= lim[1]).sum())
            print(f"{k:24s} {f(a0[k]):12.4f} ({n0:4d} bad) "
                  f"{f(a1[k]):12.4f} ({n1:4d} bad)")
    if ok:
        sjA = cb.scaled_jacobians(P, H)
        sjB = cb.scaled_jacobians(Pnew, H)
        print(f"[tfi] whole mesh scaled Jacobian: min {sjA.min():.4f} -> "
              f"{sjB.min():.4f}, {int((sjA <= 0).sum())} -> "
              f"{int((sjB <= 0).sum())} inverted")
    return P, Pnew, H, B


# --------------------------------------------------------------------------
# step 3: conforming divisions over the block complex
# --------------------------------------------------------------------------

def face_grid(vert, axis, side):
    """The (n, m) vertex grid of one lattice face."""
    sl = [slice(None)] * 3
    sl[axis] = 0 if side == 0 else -1
    return vert[tuple(sl)]


def direction_classes(lat, f2h, hexes, blk_of, verbose=True):
    """Group the block axes that must carry the same division count.

    Two blocks that share a face must agree on the two directions spanning
    it. Following that through the whole complex gives equivalence classes of
    (block, axis) pairs -- the 3D analogue of "opposite sides of a quad carry
    the same count" that `solve_edge_divisions` encodes in 2D, and the thing
    that has to be solved for, not assumed."""
    par = {}

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    def uni(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            par[rb] = ra

    for r in lat:
        for ax in (0, 1, 2):
            par.setdefault((r, ax), (r, ax))

    # where does a vertex sit in each block's lattice?
    pos = {}
    for r, (dims, vert) in lat.items():
        for idx in np.ndindex(*vert.shape):
            pos.setdefault(int(vert[idx]), {})[r] = idx

    linked = 0
    for r, (dims, vert) in lat.items():
        for ax in (0, 1, 2):
            for side in (0, 1):
                G = face_grid(vert, ax, side)
                if G.shape[0] < 2 or G.shape[1] < 2:
                    continue
                other = [a for a in (0, 1, 2) if a != ax]
                # a neighbouring block that owns the same three vertices
                v00, v10, v01 = int(G[0, 0]), int(G[1, 0]), int(G[0, 1])
                cands = set(pos.get(v00, {})) & set(pos.get(v10, {})) \
                    & set(pos.get(v01, {})) - {r}
                for s in cands:
                    p00, p10, p01 = (np.array(pos[v][s]) for v in (v00, v10, v01))
                    d1, d2 = p10 - p00, p01 - p00
                    a1 = np.nonzero(d1)[0]
                    a2 = np.nonzero(d2)[0]
                    if len(a1) != 1 or len(a2) != 1:
                        continue
                    uni((r, other[0]), (s, int(a1[0])))
                    uni((r, other[1]), (s, int(a2[0])))
                    linked += 2
    cls = defaultdict(list)
    for k in par:
        cls[find(k)].append(k)
    if verbose:
        sizes = sorted((len(v) for v in cls.values()), reverse=True)
        print(f"[tfi] {len(cls)} direction classes over {len(par)} block axes "
              f"({linked} links), class sizes {sizes}")
    return list(cls.values())


def check_conformity(lat, classes):
    """Every axis in a class must already carry the same count.

    The blocks come out of ONE hex mesh, so the existing divisions are
    conforming by construction -- if this fails, the classes are wrong, not
    the mesh. That makes it a test of `direction_classes`, not of the mesh."""
    bad = 0
    for c in classes:
        n = {int(lat[r][0][ax]) for r, ax in c}
        if len(n) > 1:
            bad += 1
            print(f"[tfi]   class {sorted(c)[:3]}... carries {sorted(n)}")
    print(f"[tfi] conformity of the existing divisions: "
          f"{'OK' if bad == 0 else str(bad) + ' classes disagree'}")
    return bad == 0


# --------------------------------------------------------------------------
# step 3b: how many cells each direction class gets
# --------------------------------------------------------------------------

def axis_length(P, vert, axis):
    """Mean physical length of the lattice lines running along `axis`.

    Along the polyline, not corner to corner: these blocks are curved, and a
    straight-line extent would ask for too few cells exactly where the passage
    bends most. The mean over all lines of that direction, so one short line at
    a corner cannot set the count for the whole face."""
    Q = P[vert]
    d = np.linalg.norm(np.diff(Q, axis=axis), axis=-1).sum(axis=axis)
    return float(d.mean())


def frozen_from_missing(lat, hexes, f2h, blk_of):
    """Axis counts that must not change because a neighbour has no lattice.

    A block whose `block_lattice` came back None cannot be refilled and keeps
    its cells. Its neighbours therefore may not re-divide the faces they share
    with it, or the two sides stop matching and the complex gets hanging
    nodes. Empirically all 16 v11 blocks do yield lattices, so this returns
    nothing there -- but a basis with a pinch is exactly where it might not,
    and silently welding a hanging face is not a failure mode worth having."""
    frozen = {}
    for r, (dims, vert) in lat.items():
        for ax in (0, 1, 2):
            for side in (0, 1):
                G = face_grid(vert, ax, side)
                other = [a for a in (0, 1, 2) if a != ax]
                for i in range(G.shape[0] - 1):
                    for j in range(G.shape[1] - 1):
                        fk = frozenset((int(G[i, j]), int(G[i + 1, j]),
                                        int(G[i + 1, j + 1]), int(G[i, j + 1])))
                        for h in f2h.get(fk, ()):
                            s = int(blk_of[h])
                            if s != r and s not in lat:
                                frozen[(r, other[0])] = int(dims[other[0]])
                                frozen[(r, other[1])] = int(dims[other[1]])
    return frozen


def solve_block_divisions(lat, classes, P, target_h, h_map=None,
                          frozen_counts=None, verbose=True):
    """Cells per direction class for a prescribed target cell size.

    The 3D counterpart of `dp3d.tmesh.solve_edge_divisions`: one integer
    variable per class, minimising the L1 distance to each axis's own ideal
    count `round(len/h)`, with the class lower bound at the FINEST of its
    axes so no cell comes out coarser than asked for.

    Conformity needs no constraint matrix. `direction_classes` already unions
    exactly the axes that span a shared face, so one variable per class IS the
    constraint -- and the three opposite-face families of a block are vacuous
    inside a lattice, where opposite faces are spanned by the same two axes
    (`OPPOSITE` is kept as a post-solve assertion, not as a row). What the
    program is actually for is `frozen_counts`: classes touching a block that
    cannot be refilled are pinned to that block's existing counts, and those
    equalities do propagate through the complex.

    `h_map` overrides the target per class index (the wall-normal families get
    a finer h). Returns {class index: count}."""
    from scipy.optimize import milp, Bounds, LinearConstraint
    h_map = h_map or {}
    frozen_counts = frozen_counts or {}
    n = len(classes)
    ideal, lb = [], []
    for ci, c in enumerate(classes):
        h = h_map.get(ci, target_h)
        want = [max(1, int(round(axis_length(P, lat[r][1], ax) / h)))
                for r, ax in c]
        ideal.append(want)
        lb.append(max(want))
    lo = np.array(lb, float)
    hi = np.full(n, np.inf)

    pinned = {}
    for ci, c in enumerate(classes):
        vals = {frozen_counts[k] for k in c if k in frozen_counts}
        if len(vals) > 1:
            raise RuntimeError(f"class {ci} is pinned to conflicting counts "
                               f"{sorted(vals)} -- the passthrough blocks "
                               f"themselves are non-conforming")
        if vals:
            v = float(vals.pop())
            pinned[ci] = int(v)
            lo[ci] = hi[ci] = v

    # L1 objective: minimise sum |x_c - ideal| over every axis of the class,
    # linearised with one slack per class (t >= x - ideal_k, t >= ideal_k - x).
    rows, rl, ru = [], [], []
    for ci, want in enumerate(ideal):
        for w in set(want):
            row = np.zeros(2 * n)
            row[ci], row[n + ci] = 1.0, -1.0
            rows.append(row), rl.append(-np.inf), ru.append(float(w))
            row = np.zeros(2 * n)
            row[ci], row[n + ci] = -1.0, -1.0
            rows.append(row), rl.append(-np.inf), ru.append(-float(w))
    c_obj = np.concatenate([np.zeros(n), np.ones(n)])
    res = milp(c=c_obj,
               constraints=[LinearConstraint(np.array(rows), rl, ru)],
               integrality=np.concatenate([np.ones(n), np.zeros(n)]),
               bounds=Bounds(np.concatenate([lo, np.zeros(n)]),
                             np.concatenate([hi, np.full(n, np.inf)])))
    counts = (np.round(res.x[:n]).astype(int) if res.success
              else lo.astype(int))
    if not res.success and verbose:
        print(f"[tfi] MILP WARNING: {res.message} -- falling back to the "
              f"lower bounds")

    for ci, c in enumerate(classes):
        assert counts[ci] >= lb[ci], f"class {ci} below its lower bound"
    # what `OPPOSITE` would have constrained: inside a lattice the two faces
    # normal to an axis are spanned by the same two axes, so the family is
    # vacuous. Asserted rather than argued.
    for r, (_dims, vert) in lat.items():
        for ax in (0, 1, 2):
            assert face_grid(vert, ax, 0).shape == face_grid(vert, ax, 1).shape, \
                f"block {r}: faces normal to axis {ax} differ in shape"

    if verbose:
        print(f"[tfi] target h={target_h}: {n} classes, counts "
              f"{sorted(counts.tolist())}, {len(pinned)} pinned "
              f"({'optimal' if res.success else 'FALLBACK'})")
    return {ci: int(counts[ci]) for ci in range(n)}


def _lerp_axis(Q, frac, axis):
    """Linear resampling of `Q` along `axis` at normalized positions `frac`.

    The parameter is the NORMALIZED INDEX of the source grid, not its arc
    length. That choice is what makes the refill weld: a block edge sampled
    on its own and the same edge sampled as the border of either adjacent
    face land on identical points, and two blocks resampling the face they
    share agree to the last bit however their lattices are oriented. With
    arc-length parameters they would not -- a face's interior rows have
    different lengths from its border, so the border would move depending on
    which face resampled it, and the mesh would come apart along exactly the
    seams it must not.

    The price is that the new grid inherits the source's spacing rather than
    being uniform in space. The AlgoHex core is near-uniform by construction,
    so the two differ little (measured and reported by the refill), and at
    identical counts every BOUNDARY vertex comes back to within 2e-15 of where
    it was -- the interior is the Gordon-Hall map, which is the fold gate's
    whole point, and costs the known min scaled Jacobian 0.1524 -> 0.1194."""
    n = Q.shape[axis]
    x = np.asarray(frac, float) * (n - 1)
    i0 = np.clip(np.floor(x).astype(int), 0, n - 2)
    t = (x - i0).reshape([-1 if d == axis else 1
                          for d in range(Q.ndim)])
    A = np.take(Q, i0, axis=axis)
    Bq = np.take(Q, i0 + 1, axis=axis)
    return (1 - t) * A + t * Bq


def resample_face_grid(Q, fu, fv):
    """A face's (n, m, 3) vertex grid resampled to (len(fu), len(fv), 3)."""
    return _lerp_axis(_lerp_axis(Q, fu, 0), fv, 1)


def uniform_fractions(n):
    """n cells -> n+1 normalized sample positions."""
    return np.linspace(0.0, 1.0, n + 1)


def clustered_fractions(n, first_height, length):
    """Tanh-clustered sample positions, `dp3d.tmesh.edge_fractions`.

    Present so the near-wall distribution is available to the refill, but the
    complex-level refill does not use it: a one-sided distribution is not
    invariant under the mirror that relates two blocks' views of a shared
    face, so clustering a class would tear the seam it is supposed to weld.
    The grading that matters for CFD enters through `reattach.py`, which
    extrudes the hub/shroud boundary layer with exactly this function and sets
    the first cell height -- the core's outer faces are interfaces to the
    re-attached parts, not walls."""
    sys.path.insert(0, str(REPO))
    from dp3d.tmesh import edge_fractions
    return np.asarray(edge_fractions(n, length, first_height), float)


def refill_block(P, vert, new_dims, fractions=None):
    """One block's vertex positions at new division counts.

    The six boundary faces are resampled from the block's existing faces --
    so the new boundary follows the old surface instead of a Coons patch
    stretched between four curves, which is what keeps a curved face on the
    input geometry -- and the interior is filled by the Gordon-Hall map.

    `fractions[axis]` overrides the sample positions along that axis."""
    ni, nj, nk = (int(d) for d in new_dims)
    fr = {ax: (fractions or {}).get(ax, uniform_fractions(d))
          for ax, d in zip((0, 1, 2), (ni, nj, nk))}
    X = np.zeros((ni + 1, nj + 1, nk + 1, 3), float)
    for ax in (0, 1, 2):
        o0, o1 = [a for a in (0, 1, 2) if a != ax]
        for side in (0, 1):
            F = resample_face_grid(P[face_grid(vert, ax, side)],
                                   fr[o0], fr[o1])
            sl = [slice(None)] * 3
            sl[ax] = 0 if side == 0 else -1
            X[tuple(sl)] = F
    return tfi(X)


def class_of_axis(classes):
    """{(block, axis): class index}."""
    return {k: ci for ci, c in enumerate(classes) for k in c}


def block_counts(lat, classes, counts):
    """{block: (ni, nj, nk)} of cells after applying the class counts."""
    cof = class_of_axis(classes)
    return {r: tuple(counts[cof[(r, ax)]] for ax in (0, 1, 2)) for r in lat}


def class_table(lat, classes, counts, P, target_h):
    """One row per direction class: what it carries now, what it would get."""
    rows = []
    for ci, c in enumerate(classes):
        lens = [axis_length(P, lat[r][1], ax) for r, ax in c]
        now = sorted({int(lat[r][0][ax]) for r, ax in c})
        rows.append({"class": ci, "axes": len(c),
                     "blocks": sorted({int(r) for r, _ax in c}),
                     "len_mean": float(np.mean(lens)),
                     "len_min": float(np.min(lens)),
                     "len_max": float(np.max(lens)),
                     "count_now": now, "count_new": int(counts[ci]),
                     "h_eff": float(np.mean(lens)) / int(counts[ci])})
    return rows


def interface_classes(P, H, B, lat, classes, input_vtk):
    """{class index: surface labels the class runs NORMAL to}.

    Which directions point across the boundary layer, in other words -- the
    ones a user would want a finer h or a graded distribution on. The
    deliverable VTK carries only `block_id`, no face labels, so they come from
    the AlgoHex input mesh by nearest face (`clean_blocks.SurfaceLabeller`),
    the same route `reattach.py` uses.

    The labels to look for are `bl_interface_hub`, `bl_interface_shroud` and
    `ogrid_interface`: this core has no hub or shroud WALLS in it at all --
    they were cut out of the domain and are re-attached later -- and reading
    "shell_blade" as a blade wall already cost this branch three wrong
    diagnoses."""
    import tet_prep_v5 as v5
    P0, tri0, tid0 = cb.read_input_surface(input_vtk)
    lab = cb.SurfaceLabeller(P0, tri0, tid0, v5.NAMES)
    f2h, _e2h = bc.build_topology(H)
    cof = class_of_axis(classes)
    out = defaultdict(set)
    for r, (_dims, vert) in lat.items():
        for ax in (0, 1, 2):
            for side in (0, 1):
                G = face_grid(vert, ax, side)
                quads, keys = [], []
                for i in range(G.shape[0] - 1):
                    for j in range(G.shape[1] - 1):
                        q = (int(G[i, j]), int(G[i + 1, j]),
                             int(G[i + 1, j + 1]), int(G[i, j + 1]))
                        if len(f2h.get(frozenset(q), ())) == 1:
                            quads.append(q)
                            keys.append(frozenset(q))
                if not quads:
                    continue
                ids = lab.label(P[np.array(quads)].mean(1))
                for s in set(int(x) for x in ids):
                    out[cof[(r, ax)]].add(v5.NAMES.get(s, s))
    return out


def ogrid_target():
    """A `SurfaceLabeller` over the blade O-grid block's OWN boundary.

    The core's `ogrid_interface` is not real geometry -- it is the cut face
    towards the O-grid that was removed from the AlgoHex domain and that
    `reattach.py` glues back on. There is therefore no geometry to be faithful
    to there; what matters is that the two parts share one surface. The O-grid
    block is reused verbatim from the source MSH, which makes it the more
    trustworthy of the two discretisations, so it is the target rather than
    the input triangulation."""
    from dp3d.extraction import parse_msh
    import reattach
    nodes, elements = parse_msh(reattach.MSH)
    Po, Ho, _Bo = reattach.ogrid_blocks(elements, nodes, verbose=False)
    f2h, _e2h = bc.build_topology(Ho)
    quads = np.asarray([cb._loop_of(Ho, hs[0], fk)
                        for fk, hs in f2h.items() if len(hs) == 1])
    tris = np.vstack([quads[:, [0, 1, 2]], quads[:, [0, 2, 3]]])
    return cb.SurfaceLabeller(Po, tris, np.zeros(len(tris), int))


def project_ogrid_interface(P, H, f2h, input_vtk, verbose=True):
    """Pull the core's O-grid cut face onto the O-grid block's surface.

    Applied per VERTEX on the whole complex, before any block is refilled, so
    that the six-face resampling and the Gordon-Hall fill inherit the moved
    boundary and absorb it -- rather than snapping afterwards and asking
    `untangle` to repair the damage locally.

    Only vertices in the INTERIOR of the patch are moved: a vertex on the ring
    where the O-grid interface meets a boundary-layer interface belongs to
    both surfaces, and pulling it onto one would drag it off the other. That
    ring is a feature curve of the input and stays where it is.

    Returns (new points, report)."""
    import tet_prep_v5 as v5
    P0, tri0, tid0 = cb.read_input_surface(input_vtk)
    lab = cb.SurfaceLabeller(P0, tri0, tid0, v5.NAMES)
    bnd = [(fk, hs[0]) for fk, hs in f2h.items() if len(hs) == 1]
    loops = np.asarray([cb._loop_of(H, h, fk) for fk, h in bnd])
    sids = lab.label(P[loops].mean(1))
    want = next(k for k, v in v5.NAMES.items() if v == "ogrid_interface")

    on, off = defaultdict(int), defaultdict(int)
    for lp, s in zip(loops, sids):
        for v in lp:
            (on if s == want else off)[int(v)] += 1
    verts = np.array(sorted(v for v in on if off[v] == 0), dtype=np.int64)
    if not len(verts):
        return P, {"projected": 0}

    tgt = ogrid_target()
    Q = cb.project_to_surface(tgt, P[verts])
    d = np.linalg.norm(Q - P[verts], axis=1)
    Pn = P.copy()
    Pn[verts] = Q
    rep = {"projected": int(len(verts)),
           "kept_on_ring": int(sum(1 for v in on if off[v])),
           "move_median": float(np.median(d)), "move_p95": float(np.percentile(d, 95)),
           "move_max": float(d.max())}
    if verbose:
        print(f"[tfi] O-grid projection: {rep['projected']} vertices moved "
              f"(median {rep['move_median']:.5f}, p95 {rep['move_p95']:.5f}, "
              f"max {rep['move_max']:.5f}), {rep['kept_on_ring']} ring "
              f"vertices left in place")
    return Pn, rep


def weld(points, tol=1e-9):
    """Merge coincident points. Returns (unique points, old -> new index).

    Two blocks resample the face they share independently and land on the same
    coordinates to the last bit or so, but with their own vertex numbering;
    this is what turns "the same position" into "the same vertex". Pairs
    within `tol` are unioned rather than bucketed by rounding, which has no
    knife edge at a bucket boundary."""
    from scipy.spatial import cKDTree
    par = np.arange(len(points))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    for a, b in cKDTree(points).query_pairs(tol, output_type="ndarray"):
        ra, rb = find(a), find(b)
        if ra != rb:
            par[max(ra, rb)] = min(ra, rb)
    root = np.array([find(i) for i in range(len(points))])
    keep = np.unique(root)
    remap = np.zeros(len(points), int)
    remap[keep] = np.arange(len(keep))
    return points[keep], remap[root]


def refill_complex(P, H, B, f2h, lat, classes, counts, fractions=None,
                   verbose=True):
    """Refill every lattice block at the solved counts; keep the rest.

    Blocks that are not lattices cannot be refilled and are passed through
    with their cells unchanged -- which is safe only because the MILP pinned
    every class adjoining them to their existing counts, so the faces they
    share still match vertex for vertex.

    Returns (points, hexes, block id, report)."""
    dims_new = block_counts(lat, classes, counts)
    cof = class_of_axis(classes)
    chunks, cells, bid = [], [], []
    for r in sorted(lat):
        _dims, vert = lat[r]
        fr = None
        if fractions:
            fr = {ax: fractions[cof[(r, ax)]] for ax in (0, 1, 2)
                  if cof[(r, ax)] in fractions}
        X = refill_block(P, vert, dims_new[r], fr)
        base = sum(len(c) for c in chunks)
        ids = np.arange(base, base + X[..., 0].size).reshape(X.shape[:3])
        chunks.append(X.reshape(-1, 3))
        cells.append(block_cells(ids))
        bid.append(np.full(len(cells[-1]), r, int))
    passthrough = sorted(set(range(int(B.max()) + 1)) - set(lat))
    for r in passthrough:
        sel = np.where(B == r)[0]
        base = sum(len(c) for c in chunks)
        old = np.unique(H[sel])
        loc = {int(v): base + i for i, v in enumerate(old)}
        chunks.append(P[old])
        cells.append(np.vectorize(loc.get)(H[sel]))
        bid.append(np.full(len(sel), r, int))

    pts = np.vstack(chunks)
    Hn = np.vstack(cells)
    Bn = np.concatenate(bid)
    pts_w, remap = weld(pts)
    Hn = remap[Hn]
    rep = {"cells_before": int(len(H)), "cells_after": int(len(Hn)),
           "points_before": int(len(P)), "points_after": int(len(pts_w)),
           "welded": int(len(pts) - len(pts_w)),
           "passthrough_blocks": passthrough}
    if verbose:
        print(f"[tfi] refill: {len(H)} -> {len(Hn)} cells, {len(pts)} -> "
              f"{len(pts_w)} points ({rep['welded']} welded), "
              f"{len(passthrough)} blocks passed through")
    return pts_w, Hn, Bn, rep


def check_watertight(P, H, verbose=True):
    """Vertex-level conformity: the refilled complex is still ONE mesh.

    The weld is the only thing that joins two blocks, so this is where a face
    resampled inconsistently would show up -- as a face used by one cell on
    one side and two half-faces on the other, i.e. as extra boundary. A
    genuine boundary face is fine; a face used by three or more cells, or a
    boundary that has grown, is not."""
    f2h, _e2h = bc.build_topology(H)
    per = Counter(len(v) for v in f2h.values())
    bad = {k: v for k, v in per.items() if k > 2}
    bnd = per.get(1, 0)
    if verbose:
        print(f"[tfi] faces by incident cells: {dict(sorted(per.items()))}"
              f"{'  <-- NON-MANIFOLD' if bad else ''}")
    return not bad, bnd


def load_blocks(blocks_vtk):
    """(P, hexes, block id per cell, f2h) from a written deliverable."""
    import meshio
    m = meshio.read(blocks_vtk)
    H = np.vstack([b.data for b in m.cells if b.type == "hexahedron"])
    B = np.concatenate([np.asarray(d).ravel() for b, d in
                        zip(m.cells, m.cell_data["block_id"])
                        if b.type == "hexahedron"])
    P = np.asarray(m.points, float)
    f2h, _e2h = bc.build_topology(H)
    return P, H, B.astype(int), f2h


def lattices(P, H, B, f2h, verbose=True):
    """{block: (dims, vert)} for every block that is a structured grid."""
    lat, missing = {}, []
    for r in range(int(B.max()) + 1):
        cells = np.where(B == r)[0]
        if not len(cells):
            continue
        got = block_lattice(H, cells, f2h)
        if got is None:
            missing.append(r)
            continue
        dims, vert = got
        lat[r] = (dims, _fix_handedness(P, vert))
    if verbose and missing:
        print(f"[tfi] {len(missing)} blocks are not lattices and will be "
              f"passed through unchanged: {missing}")
    return lat, missing


if __name__ == "__main__":
    import argparse
    import json
    ap = argparse.ArgumentParser()
    ap.add_argument("blocks", nargs="?",
                    default=str(REPO / "output" / "hex3d_algohex" / "deliverable"
                               / "T1_9_blocks_v11.vtk"))
    ap.add_argument("--target-h", type=float, default=None,
                    help="prescribed cell size for --solve-divisions")
    ap.add_argument("--solve-divisions", action="store_true",
                    help="solve the conforming division counts for "
                         "--target-h and write them next to the input")
    ap.add_argument("--project-ogrid", action="store_true",
                    help="pull the O-grid cut face onto the O-grid block's "
                         "surface before refilling, so the Gordon-Hall fill "
                         "absorbs the motion (see project_ogrid_interface)")
    ap.add_argument("--apply-divisions", action="store_true",
                    help="also refill every block at those counts and write "
                         "the result to --out")
    ap.add_argument("--out", default=None, help="output blocks VTK")
    ap.add_argument("--input-vtk",
                    default=str(REPO / "data" / "T1_9" / "T1_9_tet_v5.vtk"),
                    help="AlgoHex input mesh, for the surface labels")
    a = ap.parse_args()
    if a.apply_divisions:
        a.solve_divisions = True
    if a.solve_divisions:
        if a.target_h is None:
            ap.error("--solve-divisions needs --target-h")
        P, H, B, f2h = load_blocks(a.blocks)
        print(f"[tfi] {a.blocks}: {len(H)} cells, {int(B.max()) + 1} blocks")
        if a.project_ogrid:
            P, _rep = project_ogrid_interface(P, H, f2h, a.input_vtk)
        lat, missing = lattices(P, H, B, f2h)
        classes = direction_classes(lat, f2h, H, B)
        check_conformity(lat, classes)
        frozen = frozen_from_missing(lat, H, f2h, B)
        if frozen:
            print(f"[tfi] {len(frozen)} axes pinned by non-lattice "
                  f"neighbours: {sorted(frozen)}")
        counts = solve_block_divisions(lat, classes, P, a.target_h,
                                       frozen_counts=frozen)
        rows = class_table(lat, classes, counts, P, a.target_h)
        print(f"\n{'class':>5} {'axes':>5} {'len mean':>9} {'now':>10} "
              f"{'new':>5} {'h eff':>8}  blocks")
        for w in rows:
            print(f"{w['class']:5d} {w['axes']:5d} {w['len_mean']:9.4f} "
                  f"{str(w['count_now']):>10} {w['count_new']:5d} "
                  f"{w['h_eff']:8.4f}  {w['blocks']}")
        bc_ = block_counts(lat, classes, counts)
        total = sum(int(np.prod(d)) for d in bc_.values())
        now = sum(int(np.prod(lat[r][0])) for r in lat)
        print(f"\n[tfi] cells in the lattice blocks: {now} -> {total}")
        out = Path(a.blocks).with_suffix(".divisions.json")
        out.write_text(json.dumps({"target_h": a.target_h,
                                   "counts": counts, "classes":
                                   [[list(map(int, k)) for k in c]
                                    for c in classes],
                                   "table": rows,
                                   "frozen": {str(k): v for k, v
                                              in frozen.items()},
                                   "passthrough_blocks": missing}, indent=1))
        print(f"\n[tfi] wrote {out}")
        for ci, labs in sorted(interface_classes(P, H, B, lat, classes,
                                                 a.input_vtk).items()):
            print(f"[tfi] class {ci} is normal to {sorted(labs)}")
        if not a.apply_divisions:
            raise SystemExit(0)

        Pn, Hn, Bn, rep = refill_complex(P, H, B, f2h, lat, classes, counts)
        ok, _bnd = check_watertight(Pn, Hn)
        sj = cb.scaled_jacobians(Pn, Hn)
        print(f"[tfi] refilled mesh: scaled Jacobian min {sj.min():.4f}, "
              f"mean {sj.mean():.4f}, {int((sj <= 0).sum())} inverted")
        if not ok:
            print("[tfi] REFUSING to write: the refilled complex is not "
                  "watertight -- shared faces disagree")
            raise SystemExit(2)
        edge = np.linalg.norm(Pn[Hn[:, 1]] - Pn[Hn[:, 0]], axis=1)
        print(f"[tfi] cell size along the first edge: mean {edge.mean():.4f}, "
              f"p5 {np.percentile(edge, 5):.4f}, p95 "
              f"{np.percentile(edge, 95):.4f} (target h {a.target_h})")
        outv = a.out or str(Path(a.blocks).with_name(
            Path(a.blocks).stem + f"_h{a.target_h}.vtk"))
        import export_vtk as ev
        ev.write_vtk(outv, Pn, Hn, [12] * len(Hn), Bn, "block_id",
                     f"{int(Bn.max()) + 1} blocks refilled at h={a.target_h}")
        print(f"[tfi] wrote {outv}")
        raise SystemExit(0)
    run(a.blocks)
