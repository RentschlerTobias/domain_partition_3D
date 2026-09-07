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


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("blocks", nargs="?",
                    default=str(REPO / "output" / "hex3d_algohex" / "deliverable"
                               / "T1_9_blocks_v11.vtk"))
    a = ap.parse_args()
    run(a.blocks)


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
