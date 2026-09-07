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
    a = ap.parse_args()
    if a.solve_divisions:
        if a.target_h is None:
            ap.error("--solve-divisions needs --target-h")
        P, H, B, f2h = load_blocks(a.blocks)
        print(f"[tfi] {a.blocks}: {len(H)} cells, {int(B.max()) + 1} blocks")
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
        raise SystemExit(0)
    run(a.blocks)
