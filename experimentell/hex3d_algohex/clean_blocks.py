"""Stage 5: postprocessing of the 3D block structure (the 3D analogue of
`dp3d/clean_separatrix.py`).

The raw base complex of an AlgoHex mesh is never directly usable, exactly as
a raw separatrix network is never directly usable in 2D. `dp3d` spends most
of its code cleaning the 2D network (drop stub separatrices, snap endpoints,
collapse seam wedges); this module is the same pass one dimension up, on
sheets and blocks instead of separatrices and quads.

Three things are worth knowing before reading the code:

**1. There are two layers, and most of the module only touches the upper
one.** The cleanup is mostly a simplification of the *cut set* -- the set of
sheet faces that separate blocks. Dropping a face from the cut merges the two
blocks it separated; adding one splits a block; the hex mesh is untouched, so
those steps cannot invert a cell or move the boundary at all.

Exactly two steps change the mesh itself: `fill_cavities` (on by default;
adds the cells missing from internal holes, and every vertex it uses already
exists) and `collapse_mesh_sheets` (`--collapse-rounds`, off by default;
Gao et al.'s sheet collapse, which welds vertices and removes a layer). Both
are rejected outright if they would add an inverted cell or push the boundary
further from the input surface.

**2. Classification is fixed before anything is collapsed.** Two of the
"defects" in the raw report are counting artefacts, not geometry:

  * interior block faces were labelled by *sheet id*. One interface between
    the same two blocks can be covered by two different sheets, which splits
    a perfectly good face in two. Labelling by the *neighbouring block*
    instead is both the correct base-complex definition and what TFI needs.
  * a block may legitimately touch the same surface on two opposite sides.
    In the v9 mesh `bl_interface_hub` appears as ONE connected shell wrapping hub
    side *and* shroud side, so a passage block touches it top and bottom.
    That is a cuboid, and counting it as a defect was the artefact this
    module had to rule out first. (The shell being connected at all was
    later traced to a labelling bug in `tet_prep_v5`, fixed there; with the
    fix hub and shroud are two separate surfaces. The rule stands either
    way: only *adjacent* same-surface patches are a real defect.)

Both were measured before any collapse ran -- see `report_structure`.

Surface labels are transferred from the AlgoHex input mesh by nearest-face
lookup, never derived from a coordinate threshold: a threshold cuts across
the triangulation and has already cost this project twice (zigzag feature
curves in `tet_prep_v5`, and a drop from 74 % to 48 % cuboids in the
`block_faces.physical_of` classifier).

**3. Block merging, which the plan expected to do the work, does almost
nothing -- for a structural reason, not a tuning one.**

  * `collapse_small_sheets` -- all three sliver sheets rejected, each at
    exactly +3 excess faces. The separating sheets are bounded by singular
    edges, so dropping one merges two blocks without merging their side
    neighbours, and each of the four side faces of the union stays split in
    two.
  * `absorb_tiny_blocks` -- 0 of 8 accepted, same cause. Merging a small
    cuboid into a neighbouring cuboid cost +3 in every one of the 40
    candidate pairs examined.
  * `merge_non_cuboids` / `split_non_cuboid` -- 2 and 1 accepted.

`collapse_mesh_sheets` is what gets past this: welding a whole mesh sheet
merges the side neighbours in the same step, which is precisely what the
cut-set merge cannot express. It is also the expensive step (~13 min per
round, every sheet re-evaluated), so it is off by default.
"""

import sys
from collections import Counter, defaultdict, deque
from contextlib import contextmanager
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "experimentell" / "hex3d_algohex"))

import base_complex as bc                                             # noqa: E402
import block_faces as bfx                                             # noqa: E402
import ovm_io                                                         # noqa: E402

HF = bc.HEX_FACES
OUT = REPO / "output" / "hex3d_algohex"


# --------------------------------------------------------------------------
# surface labels: exact transfer from the AlgoHex input mesh
# --------------------------------------------------------------------------

def read_input_surface(path):
    """(points, triangles, surface id) from an AlgoHex input VTK.

    The input written by `tet_prep.write_algohex_vtk_analytic` carries one
    int CELL_DATA `color` per cell; for triangle cells that colour IS the
    surface id, so no reclassification is needed or wanted here."""
    txt = Path(path).read_text().split("\n")

    def seek(prefix, i):
        while not txt[i].startswith(prefix):
            i += 1
        return i

    i = seek("POINTS", 0)
    npts = int(txt[i].split()[1])
    P = np.array([[float(x) for x in txt[i + 1 + k].split()]
                  for k in range(npts)])
    i = seek("CELLS", i + npts)
    ncell = int(txt[i].split()[1])
    conn = [txt[i + 1 + k].split()[1:] for k in range(ncell)]
    i = seek("CELL_TYPES", i + ncell)
    types = np.array([int(txt[i + 1 + k]) for k in range(ncell)])
    i = seek("LOOKUP_TABLE", i + ncell)
    col = np.array([int(txt[i + 1 + k]) for k in range(ncell)])
    keep = np.where(types == 5)[0]                       # VTK_TRIANGLE
    tri = np.array([[int(x) for x in conn[k]] for k in keep])
    return P, tri, col[keep]


def _closest_point_dist2(q, A, B, C):
    """Squared distance from points q (n,3) to triangles (n,k,3) each."""
    return _closest_point_on_tris(q, A, B, C)[0]


def project_to_surface(labeller, Q, k=16):
    """Nearest point on the input surface for each of the points Q.

    Used to keep boundary vertices ON the geometry while they are allowed to
    slide along it -- the same constraint that made the untangling safe."""
    k = min(k, len(labeller.tris))
    _d, cand = labeller.tree.query(Q, k=k)
    cand = np.atleast_2d(cand)
    T = labeller.tris[cand]
    d2, p = _closest_point_on_tris(Q, labeller.P[T[..., 0]],
                                   labeller.P[T[..., 1]], labeller.P[T[..., 2]])
    best = np.argmin(d2, axis=1)
    return p[np.arange(len(Q)), best]


def _closest_point_on_tris(q, A, B, C):
    """(squared distance, closest point) from q (n,3) to triangles (n,k,3).

    Ericson, *Real-Time Collision Detection*, closest point on triangle,
    vectorised over the candidate axis. Using the true face distance rather
    than the centroid distance matters right at a label boundary, where the
    nearest centroid can belong to the triangle on the wrong side."""
    q = q[:, None, :]
    ab, ac, aq = B - A, C - A, q - A
    d1 = (ab * aq).sum(-1)
    d2 = (ac * aq).sum(-1)
    bq = q - B
    d3 = (ab * bq).sum(-1)
    d4 = (ac * bq).sum(-1)
    cq = q - C
    d5 = (ab * cq).sum(-1)
    d6 = (ac * cq).sum(-1)
    vc = d1 * d4 - d3 * d2
    vb = d5 * d2 - d1 * d6
    va = d3 * d6 - d5 * d4
    den = np.where(np.abs(va + vb + vc) < 1e-30, 1e-30, va + vb + vc)
    v = np.clip(vb / den, 0, 1)
    w = np.clip(vc / den, 0, 1)
    # region tests, applied in the same order as the reference algorithm
    p = A + v[..., None] * ab + w[..., None] * ac
    on_ab = (vc <= 0) & (d1 > 0) & (d3 < 0)
    t = np.clip(np.where(d1 - d3 == 0, 0.0, d1 / np.where(d1 - d3 == 0, 1, d1 - d3)), 0, 1)
    p = np.where(on_ab[..., None], A + t[..., None] * ab, p)
    on_ac = (vb <= 0) & (d2 > 0) & (d6 < 0)
    t = np.clip(np.where(d2 - d6 == 0, 0.0, d2 / np.where(d2 - d6 == 0, 1, d2 - d6)), 0, 1)
    p = np.where(on_ac[..., None], A + t[..., None] * ac, p)
    on_bc = (va <= 0) & (d4 - d3 > 0) & (d5 - d6 > 0)
    dd = (d4 - d3) + (d5 - d6)
    t = np.clip(np.where(dd == 0, 0.0, (d4 - d3) / np.where(dd == 0, 1, dd)), 0, 1)
    p = np.where(on_bc[..., None], B + t[..., None] * (C - B), p)
    p = np.where(((d1 <= 0) & (d2 <= 0))[..., None], A, p)
    p = np.where(((d3 >= 0) & (d4 <= d3))[..., None], B, p)
    p = np.where(((d6 >= 0) & (d5 <= d6))[..., None], C, p)
    return ((q - p) ** 2).sum(-1), p


class SurfaceLabeller:
    """Exact surface label for a hex-mesh boundary quad, by nearest input
    face. Never a coordinate threshold -- see the module docstring."""

    def __init__(self, points, tris, ids, names=None, k=8):
        from scipy.spatial import cKDTree
        self.P, self.tris, self.ids = points, tris, np.asarray(ids)
        self.k = min(k, len(tris))
        self.tree = cKDTree(points[tris].mean(1))
        self.names = names or {int(s): f"surf{int(s)}"
                               for s in np.unique(self.ids)}
        self.max_dist = 0.0

    def label(self, Q):
        """Q: (n,3) quad centroids -> (n,) surface ids."""
        _, cand = self.tree.query(Q, k=self.k)
        cand = np.atleast_2d(cand)
        T = self.tris[cand]
        d2 = _closest_point_dist2(Q, self.P[T[..., 0]], self.P[T[..., 1]],
                                  self.P[T[..., 2]])
        best = np.argmin(d2, axis=1)
        self.max_dist = float(np.sqrt(d2[np.arange(len(Q)), best]).max())
        return self.ids[cand[np.arange(len(Q)), best]]


def boundary_components(loops, lab=None):
    """Edge-connected components of the boundary quads, optionally required
    to share a label. Returns (component id per quad, quad adjacency)."""
    e2f = defaultdict(list)
    for i, lp in enumerate(loops):
        for k in range(4):
            u, v = lp[k], lp[(k + 1) % 4]
            e2f[(min(u, v), max(u, v))].append(i)
    adj = defaultdict(list)
    for fs in e2f.values():
        for a in range(len(fs)):
            for b in range(a + 1, len(fs)):
                adj[fs[a]].append(fs[b])
                adj[fs[b]].append(fs[a])
    comp = -np.ones(len(loops), int)
    n = 0
    for s in range(len(loops)):
        if comp[s] >= 0:
            continue
        st = [s]
        comp[s] = n
        while st:
            u = st.pop()
            for w in adj[u]:
                if comp[w] < 0 and (lab is None or lab[w] == lab[u]):
                    comp[w] = n
                    st.append(w)
        n += 1
    return comp, adj


def _rot(loop, v):
    i = loop.index(v)
    return loop[i:] + loop[:i]


def _peel_cavity(quads, vadj, max_cells=500):
    """Reconstruct the hexes missing from a cavity, by corner peeling.

    At a corner of the missing region exactly three of its faces meet at one
    vertex, and those three faces determine the cell: seven of its eight
    vertices are on them, and the eighth is the only vertex joined to at
    least two of the three face diagonals. Emit that cell, XOR its six faces
    into the surface (a face already present becomes interior and drops out),
    and repeat. The surface shrinks to nothing exactly when the region is
    filled.

    Two more obvious reconstructions do NOT work here and were tried first:
    propagating lattice coordinates by the parallelogram rule stalls after
    eight vertices, because two adjacent quads share only two; and matching
    opposite quads fails because the through-thickness edges in the middle of
    a one-cell-thick slab belong only to missing cells, so they are on
    neither the cavity surface nor the mesh.

    Returns the cells, or None if the region is not hex-fillable this way."""
    S = {tuple(q) for q in quads}
    cells = []
    for _ in range(max_cells):
        if not S:
            return cells
        by_key = {frozenset(q): q for q in S}
        v2q = defaultdict(list)
        for q in S:
            for x in q:
                v2q[x].append(q)
        found = None
        for v, qs in v2q.items():
            if len(qs) != 3:
                continue
            r = [_rot(list(q), v) for q in qs]
            for i in range(3):
                q1 = r[i]
                a, d1, b = q1[1], q1[2], q1[3]
                rest = [r[j] for j in range(3) if j != i]
                q2 = next((q for q in rest if b in (q[1], q[3])), None)
                if q2 is None:
                    continue
                if q2[3] == b:                       # orient so q2 = (v,b,d2,c)
                    q2 = [q2[0], q2[3], q2[2], q2[1]]
                d2, c = q2[2], q2[3]
                q3 = [q for q in rest if frozenset(q) != frozenset(q2)][0]
                if q3[3] == c:                       # orient so q3 = (v,c,d3,a)
                    q3 = [q3[0], q3[3], q3[2], q3[1]]
                if q3[1] != c or q3[3] != a:
                    continue
                d3 = q3[2]
                known = (v, a, b, c, d1, d2, d3)
                cand = [w for w in vadj if w not in known
                        and sum(w in vadj[d] for d in (d1, d2, d3)) >= 2]
                if len(cand) != 1:
                    continue
                found = [v, a, d1, b, c, d3, cand[0], d2]
                break
            if found:
                break
        if found is None:
            return None
        cell = found
        for f in ((cell[0], cell[1], cell[2], cell[3]),
                  (cell[4], cell[5], cell[6], cell[7]),
                  (cell[0], cell[1], cell[5], cell[4]),
                  (cell[1], cell[2], cell[6], cell[5]),
                  (cell[2], cell[3], cell[7], cell[6]),
                  (cell[3], cell[0], cell[4], cell[7])):
            k = frozenset(f)
            if k in by_key:
                S.discard(by_key[k])
                by_key.pop(k)
            else:
                S.add(tuple(f))
                by_key[k] = tuple(f)
        cells.append(cell)
    return None


def fill_cavities(P, hexes, f2h, e2h, verbose=True):
    """Refill the internal cavities of a hex mesh. Returns the new cells.

    This is the one operation in this module that changes the hex mesh, and
    it is the smallest possible such change: every vertex it uses already
    exists, so no geometry is invented and nothing is moved or welded -- the
    missing cells are simply reconnected. It is rejected outright if any new
    cell is inverted, which is the hard constraint on the whole stage."""
    loops, owner = [], []
    for fk, hs in f2h.items():
        if len(hs) != 1:
            continue
        c = hexes[hs[0]]
        for fc in HF:
            lp = [int(c[i]) for i in fc]
            if frozenset(lp) == fk:
                loops.append(lp)
                owner.append(hs[0])
                break
    comp, _adj = boundary_components(loops)
    sizes = Counter(comp.tolist())
    outer = sizes.most_common(1)[0][0]
    new = []
    for cid, n in sorted(sizes.items(), key=lambda kv: -kv[1]):
        if cid == outer:
            continue
        quads = [loops[i] for i in np.where(comp == cid)[0]]
        vs = {x for q in quads for x in q}
        vadj = defaultdict(set)
        for q in quads:
            for k in range(4):
                vadj[q[k]].add(q[(k + 1) % 4])
                vadj[q[(k + 1) % 4]].add(q[k])
        for g in e2h:                     # mesh edges between cavity vertices
            if g[0] in vs and g[1] in vs:
                vadj[g[0]].add(g[1])
                vadj[g[1]].add(g[0])
        cells = _peel_cavity(quads, vadj)
        if cells is None:
            if verbose:
                print(f"[clean_blocks] cavity of {n} quads: NOT fillable by "
                      f"corner peeling -- left open")
            continue
        if verbose:
            print(f"[clean_blocks] cavity of {n} quads -> {len(cells)} cells")
        new += cells
    if not new:
        return np.zeros((0, 8), np.int64)
    new = np.array([c if ovm_io._hex_volume(P[c]) > 0 else c[4:] + c[:4]
                    for c in new], np.int64)
    sj = np.array([ovm_io.scaled_jacobian(P, c) for c in new])
    if verbose:
        print(f"[clean_blocks] fill_cavities: {len(new)} cells, scaled "
              f"Jacobian min {sj.min():.4f} mean {sj.mean():.4f}")
    if sj.min() <= 0:
        print(f"[clean_blocks] REJECTED: fill would add {int((sj <= 0).sum())} "
              f"inverted cells")
        return np.zeros((0, 8), np.int64)
    return new


# --------------------------------------------------------------------------
# mesh-level sheet collapse (Gao et al. 2017)
# --------------------------------------------------------------------------

# the three parallel edge classes of a hex in VTK_HEXAHEDRON ordering
_PARALLEL = ([(0, 1), (3, 2), (7, 6), (4, 5)],
             [(0, 3), (1, 2), (5, 6), (4, 7)],
             [(0, 4), (1, 5), (2, 6), (3, 7)])


def scaled_jacobians(P, hexes):
    """Vectorised `ovm_io.scaled_jacobian` over a whole mesh. Same corner
    tetrahedra; the scalar version is ~100x too slow to screen sheets."""
    idx = [(0, 1, 3, 4), (1, 2, 0, 5), (2, 3, 1, 6), (3, 0, 2, 7),
           (4, 7, 5, 0), (5, 4, 6, 1), (6, 5, 7, 2), (7, 6, 4, 3)]
    Q = P[hexes]
    out = np.empty((len(hexes), 8))
    for k, (a, b, c, d) in enumerate(idx):
        u, v, w = Q[:, b] - Q[:, a], Q[:, c] - Q[:, a], Q[:, d] - Q[:, a]
        n = (np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1)
             * np.linalg.norm(w, axis=1))
        out[:, k] = np.einsum('ij,ij->i', u, np.cross(v, w)) / np.where(n > 0, n, 1.0)
    return out.min(1)


def mesh_sheets(hexes):
    """The sheets of the mesh, as equivalence classes of parallel edges.

    A hex mesh is a stack of intersecting sheets; two edges belong to the
    same sheet if they are parallel inside some hex. This is the dual view
    Gao et al. simplify, and it is a different object from the *separating*
    sheets of the base complex in `base_complex.sheet_faces` -- those are
    bounded by singular edges, these run through the whole mesh."""
    par = {}

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    for c in hexes:
        for grp in _PARALLEL:
            es = [(min(int(c[i]), int(c[j])), max(int(c[i]), int(c[j])))
                  for i, j in grp]
            for e in es:
                par.setdefault(e, e)
            for e in es[1:]:
                ra, rb = find(es[0]), find(e)
                if ra != rb:
                    par[rb] = ra
    cls = defaultdict(list)
    for e in par:
        cls[find(e)].append(e)
    return list(cls.values())


def collapse_sheet(P, hexes, edges):
    """Remove one sheet: weld the two endpoints of each of its edges.

    Every hex of the sheet degenerates and drops out; every other hex keeps
    eight distinct vertices. Welded vertices go to the centroid of their
    class. Returns (points, hexes, n_dropped)."""
    vp = {}

    def find(x):
        vp.setdefault(x, x)
        while vp[x] != x:
            vp[x] = vp[vp[x]]
            x = vp[x]
        return x

    for u, v in edges:
        ru, rv = find(u), find(v)
        if ru != rv:
            vp[rv] = ru
    grp = defaultdict(list)
    for v in range(len(P)):
        grp[find(v)].append(v)
    Q = P.copy()
    for _r, vs in grp.items():
        if len(vs) > 1:
            Q[vs] = P[vs].mean(0)
    keep, dropped = [], 0
    for c in hexes:
        m = [find(int(v)) for v in c]
        if len(set(m)) == 8:
            keep.append(m)
        else:
            dropped += 1
    return Q, np.array(keep, np.int64), dropped


def _structure_stats(P, hexes, labeller):
    """Block-structure summary of a candidate mesh, without building a full
    BlockStructure (this runs once per sheet per round)."""
    import contextlib
    import io
    with contextlib.redirect_stdout(io.StringIO()):
        f2h, e2h = bc.build_topology(hexes)
        sing = bc.singular_edges(hexes, P, f2h, e2h)
        sheet_of, _n = bfx.label_sheets(hexes, f2h, sing)
        blocks = bc.blocks_from_cut(hexes, f2h, set(sheet_of))
    blk = np.zeros(len(hexes), int)
    for i, b in enumerate(blocks):
        for h in b:
            blk[h] = i
    loops, keys = [], []
    for fk, hs in f2h.items():
        if len(hs) != 1:
            continue
        c = hexes[hs[0]]
        for fc in HF:
            lp = [int(c[i]) for i in fc]
            if frozenset(lp) == fk:
                loops.append(lp)
                keys.append(fk)
                break
    ids = labeller.label(P[np.asarray(loops)].mean(1))
    surf = {k: int(s) for k, s in zip(keys, ids)}
    cub = exc = 0
    sizes = []
    for bi, b in enumerate(blocks):
        rec = []
        for h in b:
            c = hexes[h]
            for fc in HF:
                lp = [int(c[i]) for i in fc]
                k = frozenset(lp)
                nb = [x for x in f2h[k] if x != h]
                if not nb:
                    rec.append((("P", surf[k]), lp))
                elif blk[nb[0]] != bi:
                    rec.append((("B", int(blk[nb[0]])), lp))
        p = _patch_split(rec)
        cub += cuboid_status(p)[0] == "cuboid"
        exc += _excess(p)
        sizes.append(len(b))
    sj = scaled_jacobians(P, hexes)
    return {"cells": len(hexes), "blocks": len(blocks), "cuboids": cub,
            "excess": exc, "sing": len(sing),
            "tiny": int(sum(1 for s in sizes if s < 10)),
            "inverted": int((sj <= 0).sum()),
            "min_sj": float(sj.min()),
            "hausdorff": float(labeller.max_dist)}


def collapse_mesh_sheets(P, hexes, labeller, max_rounds=4, verbose=True):
    """Greedy sheet collapse -- the operation the merge-only cleanup cannot
    express (see the module docstring).

    Accept a sheet only if it adds NO inverted cell, does not push the
    boundary further from the input surface, and strictly improves
    (excess faces, then block count). Returns (P, hexes, log)."""
    base = _structure_stats(P, hexes, labeller)
    if verbose:
        print(f"[clean_blocks] sheet collapse, start: {base}")
    log = [dict(base, sheet=None)]
    for rnd in range(max_rounds):
        sheets = mesh_sheets(hexes)
        best, best_state = None, None
        for i, es in enumerate(sheets):
            Q, H, _drop = collapse_sheet(P, hexes, es)
            if not len(H):
                continue
            if scaled_jacobians(Q, H).min() <= 0 and base["inverted"] == 0:
                continue
            st = _structure_stats(Q, H, labeller)
            if st["inverted"] > base["inverted"]:
                continue
            # and never let the worst cell get worse: the count alone can
            # stay flat while a cell folds further
            if st["min_sj"] < base["min_sj"] - 1e-9:
                continue
            if st["hausdorff"] > base["hausdorff"] + 1e-9:
                continue
            key = (st["excess"], st["blocks"])
            if key < (base["excess"], base["blocks"]) and \
                    (best is None or key < best):
                best, best_state = key, (i, Q, H, st)
        if best_state is None:
            if verbose:
                print(f"[clean_blocks] round {rnd + 1}: no sheet improves, stop")
            break
        i, P, hexes, base = best_state
        if verbose:
            print(f"[clean_blocks] round {rnd + 1}: collapsed sheet {i} "
                  f"({len(sheets)} candidates) -> {base}")
        log.append(dict(base, sheet=i))
    return P, hexes, log


def untangle(P, hexes, f2h, labeller=None, rings=1, verbose=True):
    """Gao et al.'s second half: repair element quality after simplification.

    Only INTERIOR vertices are allowed to move, which is what makes this
    safe to bolt onto the rest of the module: the domain boundary is
    untouched to the last bit, so the Hausdorff distance to the input surface
    cannot change, and connectivity is untouched, so the block topology
    cannot change either. All this can do is move inner vertices of a few
    cells around the inverted ones.

    Note what this is NOT here: in Gao's pipeline the smoothing repairs
    damage done by the collapses. Our collapses did no damage -- none of the
    145 sheets added an inverted cell. The only thing left to repair is a
    defect inherited from HexEx.

    The objective is a one-sided barrier, sum of max(0, eps - J)^2 over the
    affected cells: it pushes every cell above eps and leaves cells that are
    already good alone. Maximising the soft MINIMUM instead -- the obvious
    first choice -- is actively harmful here: it raised the worst cell from
    -0.0582 to -0.0453 while dragging three of its neighbours below zero, so
    the mesh went from 1 inverted cell to 4. The acceptance test therefore
    ranks (number of inverted cells, worst value), never the worst value on
    its own."""
    sj = scaled_jacobians(P, hexes)
    bad = np.where(sj <= 0)[0]
    if not len(bad):
        if verbose:
            print("[clean_blocks] untangle: no inverted cells, nothing to do")
        return P, 0
    bnd_v = {v for fk, hs in f2h.items() if len(hs) == 1 for v in fk}
    v2c = defaultdict(list)
    for i, c in enumerate(hexes):
        for v in c:
            v2c[int(v)].append(i)

    seed = {int(v) for b in bad for v in hexes[b]}
    for _ in range(rings - 1):
        seed = {int(v) for x in seed for c in v2c[x] for v in hexes[c]}
    free = sorted(v for v in seed if v not in bnd_v)
    if not free:
        if verbose:
            print("[clean_blocks] untangle: inverted cells have no interior "
                  "vertex to move -- would need boundary motion, skipped")
        return P, 0
    cells = sorted({c for v in free for c in v2c[v]})
    loc = {v: i for i, v in enumerate(sorted({int(x) for c in cells
                                              for x in hexes[c]}))}
    Pl = np.array([P[v] for v in loc])
    Hl = np.array([[loc[int(x)] for x in hexes[c]] for c in cells])
    fidx = np.array([loc[v] for v in free])
    g0 = (int((sj <= 0).sum()), -sj.min())
    best = _smooth(P, hexes, Pl, Hl, free, fidx, cells, g0, None, verbose)
    # Escalate whenever cells are still inverted, NOT merely when the
    # interior sweep failed outright: that sweep does improve the worst value
    # (-0.0582 -> -0.0451) while leaving the cell inverted, and keying the
    # escalation on "did it improve" let that count as success and stopped
    # the search one step short of the actual repair.
    if (best is None or best[0][0] > 0) and labeller is not None:
        if best is not None:
            P, g0 = best[1], best[0]
        # Escalation. Interior motion alone cannot untangle a cell whose
        # inversion is already present in its boundary quads -- measured:
        # it lifts the worst cell from -0.0582 only to -0.0349. Letting the
        # boundary vertices SLIDE, penalised by their distance to the input
        # surface, does untangle it, and in practice they stay on the
        # surface to 5 decimals, so the domain is not deformed.
        free = sorted(seed)
        cells = sorted({c for v in free for c in v2c[v]})
        loc = {v: i for i, v in enumerate(sorted({int(x) for c in cells
                                                  for x in hexes[c]}))}
        Pl = np.array([P[v] for v in loc])
        Hl = np.array([[loc[int(x)] for x in hexes[c]] for c in cells])
        fidx = np.array([loc[v] for v in free])
        on_b = np.where([v in bnd_v for v in free])[0]
        if verbose:
            print(f"[clean_blocks]   interior-only failed, letting "
                  f"{len(on_b)} boundary vertices slide on the input surface")
        slid = _smooth(P, hexes, Pl, Hl, free, fidx, cells, g0,
                       (labeller, on_b, np.array([P[v] for v in free])[on_b]),
                       verbose)
        if slid is not None:
            best = slid
    if best is None:
        if verbose:
            print(f"[clean_blocks] untangle: nothing improves "
                  f"({g0[0]} inverted, min {-g0[1]:.4f}), rejected")
        return P, 0
    g, P, drift = best
    if verbose:
        print(f"[clean_blocks] untangle: {len(free)} vertices moved over "
              f"{len(cells)} cells; inverted {g0[0]} -> {g[0]}, min scaled "
              f"Jacobian {-g0[1]:.4f} -> {-g[1]:.4f}, boundary drift {drift:.6f}")
    return P, len(free)


def _smooth(P, hexes, Pl, Hl, free, fidx, cells, g0, slide, verbose,
            max_drift=1e-3):
    """One smoothing sweep. `slide` is None for interior-only, else
    (labeller, indices into `free` that are on the boundary, their positions).

    The objective is one-sided: sum of max(0, eps - J)^2, plus, when sliding,
    a penalty on the squared distance of each moved boundary vertex to the
    input surface."""
    from scipy.optimize import minimize
    A = B = C = None
    if slide is not None:
        labeller, on_b, bpos0 = slide
        _d, cand = labeller.tree.query(bpos0, k=min(24, len(labeller.tris)))
        T3 = labeller.tris[cand]
        A, B, C = (labeller.P[T3[..., 0]], labeller.P[T3[..., 1]],
                   labeller.P[T3[..., 2]])

    def drift_of(x):
        if slide is None:
            return np.zeros(1)
        q = x.reshape(-1, 3)[on_b]
        return np.sqrt(_closest_point_dist2(q, A, B, C).min(1))

    best = None
    grid = ((0.05, 0.0), (0.15, 0.0), (0.30, 0.0)) if slide is None else \
           ((0.15, 50.0), (0.30, 50.0), (0.30, 200.0))
    for eps, w in grid:
        def obj(x, eps=eps, w=w):
            Q = Pl.copy()
            Q[fidx] = x.reshape(-1, 3)
            s = scaled_jacobians(Q, Hl)
            f = float((np.maximum(0.0, eps - s) ** 2).sum())
            if w:
                q = x.reshape(-1, 3)[on_b]
                f += w * float(_closest_point_dist2(q, A, B, C).min(1).sum())
            return f

        res = minimize(obj, Pl[fidx].ravel(), method="L-BFGS-B",
                       options={"maxiter": 800})
        Q = Pl.copy()
        Q[fidx] = res.x.reshape(-1, 3)
        cand_P = P.copy()
        for v, i in zip(free, fidx):
            cand_P[v] = Q[i]
        s = scaled_jacobians(cand_P, hexes)
        g = (int((s <= 0).sum()), -s.min())
        d = float(drift_of(res.x).max())
        if verbose:
            print(f"[clean_blocks]   eps={eps} w={w}: {g[0]} inverted, "
                  f"min {-g[1]:.4f}, boundary drift {d:.6f}")
        if d > max_drift:
            continue
        if g < g0 and (best is None or g < best[0]):
            best = (g, cand_P, d)
    return best


def label_cavities(S, verbose=True):
    """Give every non-outer boundary component its own label.

    Measured on v9 and NOT documented anywhere before: the hex mesh has two
    internal cavities, closed quad surfaces of 28 and 22 faces enclosing
    0.00170 and 0.00145 volume. That is where the "99.8 % coverage" goes.
    Both are structured boxes -- 8 valence-3 corners, all other vertices
    valence 4, so 28 quads = a 1x2x4 grid and 22 quads = 1x1x5 -- and all of
    their vertices already exist in the mesh, so an exact refill is possible
    in principle. It is NOT done here: refilling changes the hex mesh, which
    is outside this module's remit (nothing else here moves or adds a cell).

    Labelling them explicitly matters because nearest-face lookup otherwise
    gives a cavity quad the label of whatever input surface happens to be
    closest, which is meaningless and hides the defect."""
    comp, _adj = boundary_components(S.bnd_loops)
    sizes = Counter(comp.tolist())
    outer = sizes.most_common(1)[0][0]
    nxt = max(S.surf_names) + 1 if S.surf_names else 1
    cav = []
    for c, n in sorted(sizes.items(), key=lambda kv: -kv[1]):
        if c == outer:
            continue
        S.surf_names[nxt] = f"void{len(cav)}"
        for i in np.where(comp == c)[0]:
            S.surf_of[frozenset(S.bnd_loops[i])] = nxt
        cav.append((nxt, n))
        nxt += 1
    S.cavities = cav
    S._cache.clear()
    if verbose and cav:
        print(f"[clean_blocks] INTERNAL CAVITIES: {len(cav)} "
              f"({', '.join(f'{n} quads' for _s, n in cav)}) -- the hex mesh "
              f"is not solid; these blocks cannot be cuboids until refilled")
    return cav


def remove_label_islands(S, max_quads=64, verbose=True):
    """Clean the transferred labels: a same-label component that is entirely
    surrounded by ONE strictly larger label is an island and is relabelled.

    The hex boundary and the input triangulation are two different
    discretisations of the same geometry, so the transferred label boundary
    can wander by a quad and pinch off an island, which then splits a
    block's contact with a surface into two patches and reads as a
    non-cuboid. Purely topological -- an island is defined by the mesh's own
    quad adjacency, never by position.

    The "strictly larger" rule is load-bearing: without it two mutually
    bordering components (which is what a cavity's two label patches are)
    swap labels forever."""
    loops = S.bnd_loops
    total = 0
    for _sweep in range(10):
        lab = np.array([S.surf_of[frozenset(lp)] for lp in loops])
        comp, adj = boundary_components(loops, lab)
        sizes = Counter(comp.tolist())
        moved = 0
        for c, n in sizes.items():
            if n > max_quads:
                continue
            m = np.where(comp == c)[0]
            ring = Counter()
            for i in m:
                for w in adj[i]:
                    if comp[w] != c:
                        ring[comp[w]] += 1
            if len(ring) != 1:
                continue                     # ambiguous, or a closed surface
            other = next(iter(ring))
            if sizes[other] <= n:
                continue
            tgt = int(lab[np.where(comp == other)[0][0]])
            for i in m:
                S.surf_of[frozenset(loops[i])] = tgt
            if verbose:
                print(f"[clean_blocks] island: {n} quads "
                      f"{S.surf_names.get(int(lab[m[0]]))} -> "
                      f"{S.surf_names.get(tgt)}")
            moved += n
        total += moved
        if not moved:
            break
    S._cache.clear()
    if verbose:
        print(f"[clean_blocks] remove_label_islands: {total} quads relabelled")
    return total


def feature_angle_labels(P, loops, sharp_deg=40.0):
    """Fallback when no input mesh is given: segment the boundary into
    patches at sharp edges. Still a mesh-intrinsic criterion (a dihedral
    angle between two existing faces), not a coordinate cut. Used for the
    AlgoHex cylinder demo, which has no labelled input here."""
    def nrm(lp):
        Q = P[lp]
        n = np.cross(Q[2] - Q[0], Q[3] - Q[1])
        ln = np.linalg.norm(n)
        return n / ln if ln else n

    nn = [nrm(lp) for lp in loops]
    e2f = defaultdict(list)
    for i, lp in enumerate(loops):
        for k in range(4):
            u, v = lp[k], lp[(k + 1) % 4]
            e2f[(min(u, v), max(u, v))].append(i)
    cthr = np.cos(np.radians(sharp_deg))
    adj = defaultdict(list)
    for fs in e2f.values():
        for a in range(len(fs)):
            for b in range(a + 1, len(fs)):
                if float(np.dot(nn[fs[a]], nn[fs[b]])) >= cthr:
                    adj[fs[a]].append(fs[b])
                    adj[fs[b]].append(fs[a])
    lab = -np.ones(len(loops), int)
    n = 0
    for s in range(len(loops)):
        if lab[s] >= 0:
            continue
        st = [s]
        lab[s] = n
        while st:
            u = st.pop()
            for w in adj[u]:
                if lab[w] < 0:
                    lab[w] = n
                    st.append(w)
        n += 1
    return lab, {i: f"patch{i}" for i in range(n)}


# --------------------------------------------------------------------------
# the block structure
# --------------------------------------------------------------------------

def _patch_split(items):
    """Group (key, loop) face records into maximal edge-connected patches of
    equal key. Returns [(key, [loops])]."""
    e2f = defaultdict(list)
    for fi, (_k, lp) in enumerate(items):
        for a in range(4):
            u, v = lp[a], lp[(a + 1) % 4]
            e2f[(min(u, v), max(u, v))].append(fi)
    adj = defaultdict(list)
    for fl in e2f.values():
        for x in range(len(fl)):
            for y in range(x + 1, len(fl)):
                if items[fl[x]][0] == items[fl[y]][0]:
                    adj[fl[x]].append(fl[y])
                    adj[fl[y]].append(fl[x])
    seen, out = set(), []
    for s in range(len(items)):
        if s in seen:
            continue
        comp, st = [], [s]
        seen.add(s)
        while st:
            u = st.pop()
            comp.append(u)
            for w in adj[u]:
                if w not in seen:
                    seen.add(w)
                    st.append(w)
        out.append((items[s][0], [items[i][1] for i in comp]))
    return out


class BlockStructure:
    """Base complex of a hex mesh, plus the merge operations of the cleanup.

    A merge is a union-find union of two blocks: the interface between them
    stops separating. The hex mesh is never modified."""

    def __init__(self, ovm_path, labeller=None, sharp_deg=40.0, fill=True,
                 collapse_rounds=0, untangle_mesh=False):
        P, e, f, poly = ovm_io.read_ovm(ovm_path)
        self.P = P
        self.hexes, skipped = ovm_io.ovm_to_cells(P, e, f, poly)
        self.f2h, self.e2h = bc.build_topology(self.hexes)
        self.n_filled = 0
        if fill:
            add = fill_cavities(P, self.hexes, self.f2h, self.e2h)
            if len(add):
                self.hexes = np.vstack([self.hexes, add])
                self.n_filled = len(add)
                self.f2h, self.e2h = bc.build_topology(self.hexes)
        self.n_untangled = 0
        # Untangle BEFORE collapsing, not only after. The collapse guard
        # compares against the starting inverted-cell count, so starting from
        # 96 makes it toothless: a collapse that keeps 96 is waved through.
        # Repairing first gives the guard something to protect.
        if untangle_mesh:
            self.P, n = untangle(self.P, self.hexes, self.f2h, labeller)
            self.n_untangled += n
            P = self.P
        self.collapse_log = []
        if collapse_rounds and labeller is not None:
            self.P, self.hexes, self.collapse_log = collapse_mesh_sheets(
                self.P, self.hexes, labeller, max_rounds=collapse_rounds)
            P = self.P
            self.f2h, self.e2h = bc.build_topology(self.hexes)
            # A collapse welds vertices and can open a hole. fill_cavities
            # used to run only BEFORE it, so such a hole would never have
            # been repaired -- it did not happen on v11, but the ordering
            # was wrong.
            add = fill_cavities(P, self.hexes, self.f2h, self.e2h)
            if len(add):
                self.hexes = np.vstack([self.hexes, add])
                self.n_filled += len(add)
                self.f2h, self.e2h = bc.build_topology(self.hexes)
                print(f"[clean_blocks] {len(add)} cells filled AFTER the "
                      f"collapse (a collapse had opened a cavity)")
        if untangle_mesh:
            self.P, n = untangle(self.P, self.hexes, self.f2h, labeller)
            self.n_untangled += n
            P = self.P
        self.sing = bc.singular_edges(self.hexes, P, self.f2h, self.e2h)
        self.sheet_of, self.n_sheet = bfx.label_sheets(self.hexes, self.f2h,
                                                       self.sing)
        base = bc.blocks_from_cut(self.hexes, self.f2h, set(self.sheet_of))
        self.blk0 = np.zeros(len(self.hexes), int)
        for i, b in enumerate(base):
            for h in b:
                self.blk0[h] = i
        self.n_base = len(base)
        self.parent = list(range(self.n_base))
        self.cut = set(self.sheet_of)
        self.extra_cut = set()

        # boundary quads and their exact surface labels
        self.bnd = [(fk, hs[0]) for fk, hs in self.f2h.items() if len(hs) == 1]
        loops = []
        for fk, h in self.bnd:
            c = self.hexes[h]
            for fc in HF:
                lp = [int(c[i]) for i in fc]
                if frozenset(lp) == fk:
                    loops.append(lp)
                    break
        self.bnd_loops = loops
        if labeller is not None:
            ids = labeller.label(P[np.asarray(loops)].mean(1))
            self.surf_names = labeller.names
            self.label_source = f"input mesh, nearest face " \
                                f"(max dist {labeller.max_dist:.4f})"
        else:
            ids, self.surf_names = feature_angle_labels(P, loops, sharp_deg)
            self.label_source = f"feature angle {sharp_deg:.0f} deg (no input mesh)"
        self.surf_of = {fk: int(s) for (fk, _h), s in zip(self.bnd, ids)}
        self._cache = {}

    @classmethod
    def from_arrays(cls, P, hexes, blk_of, labeller):
        """Rebuild from an already-partitioned mesh (points, hexes and a
        block id per cell), skipping the OVM read and the base complex.

        Used to re-derive the block faces of a written deliverable without
        rerunning the pipeline -- the block partition is the expensive part
        and it is already stored as `block_id`."""
        S = cls.__new__(cls)
        S.P, S.hexes = np.asarray(P), np.asarray(hexes)
        S.f2h, S.e2h = bc.build_topology(S.hexes)
        S.n_filled = S.n_untangled = 0
        S.collapse_log = []
        S.sheet_of, S.n_sheet = {}, 0
        S.sing = set()
        S.blk0 = np.asarray(blk_of, int).ravel().copy()
        S.n_base = int(S.blk0.max()) + 1
        S.parent = list(range(S.n_base))
        S.cut = set()
        S.extra_cut = set()
        S.bnd = [(fk, hs[0]) for fk, hs in S.f2h.items() if len(hs) == 1]
        loops = []
        for fk, h in S.bnd:
            c = S.hexes[h]
            for fc in HF:
                lp = [int(c[i]) for i in fc]
                if frozenset(lp) == fk:
                    loops.append(lp)
                    break
        S.bnd_loops = loops
        ids = labeller.label(S.P[np.asarray(loops)].mean(1))
        S.surf_names = labeller.names
        S.label_source = "input mesh, nearest face"
        S.surf_of = {fk: int(s) for (fk, _h), s in zip(S.bnd, ids)}
        S.cavities = []
        S._cache = {}
        return S

    # -- union-find ---------------------------------------------------------
    def root(self, b):
        p = self.parent
        while p[b] != b:
            p[b] = p[p[b]]
            b = p[b]
        return b

    def _roots(self):
        return sorted({self.root(b) for b in range(self.n_base)})

    def cells_of(self):
        d = defaultdict(list)
        for h in range(len(self.hexes)):
            d[self.root(self.blk0[h])].append(h)
        return d

    def merge(self, a, b):
        ra, rb = self.root(a), self.root(b)
        if ra == rb:
            return False
        self.parent[rb] = ra
        self._cache.clear()      # neighbour keys change too, so drop all
        return True

    # -- splitting ----------------------------------------------------------
    def _face_keys(self, h):
        c = self.hexes[h]
        return [frozenset(int(c[i]) for i in fc) for fc in HF]

    def components_within(self, cells, extra):
        """Connected components of one block's cells once `extra` faces are
        made separating. Staying inside `cells` is what keeps the existing
        block boundary intact."""
        cs = set(cells)
        seen, out = set(), []
        for s in cells:
            if s in seen:
                continue
            st, comp = [s], []
            seen.add(s)
            while st:
                u = st.pop()
                comp.append(u)
                for k in self._face_keys(u):
                    if k in extra:
                        continue
                    nb = [x for x in self.f2h[k] if x != u]
                    if nb and nb[0] in cs and nb[0] not in seen:
                        seen.add(nb[0])
                        st.append(nb[0])
            out.append(comp)
        return out

    def split(self, r, extra):
        """Make `extra` separating inside block r. Returns the new roots."""
        cells = self.cells_of()[r]
        comps = self.components_within(cells, extra)
        if len(comps) < 2:
            return None
        self.extra_cut |= set(extra)
        self.cut |= set(extra)
        roots = [r]
        for comp in comps[1:]:
            nid = self.n_base
            self.n_base += 1
            self.parent.append(nid)
            for h in comp:
                self.blk0[h] = nid
            roots.append(nid)
        self._cache.clear()
        return roots

    # -- faces --------------------------------------------------------------
    def face_records(self, cells, me):
        """(key, loop) for every face of `cells` that bounds the block."""
        rec = []
        for h in cells:
            c = self.hexes[h]
            for fc in HF:
                lp = [int(c[i]) for i in fc]
                k = frozenset(lp)
                nb = [x for x in self.f2h[k] if x != h]
                if not nb:
                    rec.append((("P", self.surf_of[k]), lp))
                elif self.root(self.blk0[nb[0]]) != me:
                    rec.append((("B", self.root(self.blk0[nb[0]])), lp))
        return rec

    def patches(self, root=None, cells=None):
        if cells is None:
            if root in self._cache:
                return self._cache[root]
            cells = self.cells_of()[root]
            out = _patch_split(self.face_records(cells, root))
            self._cache[root] = out
            return out
        return _patch_split(self.face_records(cells, root))

    def all_patches(self):
        return {r: self.patches(r) for r, _ in self.cells_of().items()}

    def neighbours(self, root, cells):
        return {k[1] for k, _ in self.patches(root=root) if k[0] == "B"}


# --------------------------------------------------------------------------
# cuboid test -- classification, before any geometry is touched
# --------------------------------------------------------------------------

def _patch_adjacency(patches):
    """Which patches share an edge (the block's face adjacency graph)."""
    e2p = defaultdict(set)
    for pi, (_k, loops) in enumerate(patches):
        for lp in loops:
            for a in range(4):
                u, v = lp[a], lp[(a + 1) % 4]
                e2p[(min(u, v), max(u, v))].add(pi)
    adj = defaultdict(set)
    for ps in e2p.values():
        for x in ps:
            for y in ps:
                if x != y:
                    adj[x].add(y)
    return adj


def cuboid_status(patches):
    """('cuboid' | reason, detail).

    A block is a topological cuboid iff it has exactly 6 faces. Two faces
    may carry the same surface label as long as they are NOT adjacent -- a
    block that touches one surface on two opposite sides is a legitimate
    cuboid, and counting that as a defect was the counting artefact this
    module had to rule out first."""
    n = len(patches)
    if n != 6:
        return f"{n} faces", n
    adj = _patch_adjacency(patches)
    for i in range(6):
        if len(adj[i]) != 4:
            return "not cube-adjacent", len(adj[i])
    for i in range(6):
        for j in adj[i]:
            if patches[i][0] == patches[j][0]:
                return "adjacent same-surface faces", patches[i][0]
    return "cuboid", 6


def duplicate_surface_faces(patches):
    """Surface labels used by more than one face of the block, split into
    opposite (legitimate) and adjacent (a real defect)."""
    adj = _patch_adjacency(patches)
    by_lab = defaultdict(list)
    for i, (k, _l) in enumerate(patches):
        if k[0] == "P":
            by_lab[k[1]].append(i)
    opp, bad = {}, {}
    for lab, ps in by_lab.items():
        if len(ps) < 2:
            continue
        if any(b in adj[a] for a in ps for b in ps if a != b):
            bad[lab] = len(ps)
        else:
            opp[lab] = len(ps)
    return opp, bad


# --------------------------------------------------------------------------
# step 1: diagnose
# --------------------------------------------------------------------------

def touches_cavity(S, patches):
    """Faces of this block that lie on an internal cavity."""
    voids = {s for s, _n in getattr(S, "cavities", [])}
    return [i for i, (k, _l) in enumerate(patches)
            if k[0] == "P" and k[1] in voids]


def report_structure(S, title="structure", verbose=True):
    cells = S.cells_of()
    pat = {r: S.patches(r) for r in cells}
    sizes = sorted((len(c) for c in cells.values()), reverse=True)
    nf = Counter(len(p) for p in pat.values())
    stat = {r: cuboid_status(p) for r, p in pat.items()}
    cub = [r for r, (s, _d) in stat.items() if s == "cuboid"]
    ncell_cub = sum(len(cells[r]) for r in cub)
    n = len(cells)
    sj = np.array([ovm_io.scaled_jacobian(S.P, c) for c in S.hexes])
    # projection, clearly separated from the measurement: how many blocks
    # are non-cuboid ONLY because an internal cavity chips extra faces off
    # them. Not a result -- it says what a cavity refill would be worth.
    # `len(p) - len(cavity faces) == 6` alone is NOT enough: for a block with
    # 6 faces that fails the cube-adjacency test it reduces to 6 - 0 == 6 and
    # blamed an internal cavity that the block does not even touch. Require
    # that it actually touches one.
    cav_only = [r for r, p in pat.items()
                if stat[r][0] != "cuboid"
                and len(touches_cavity(S, p)) > 0
                and len(p) - len(touches_cavity(S, p)) == 6]
    rep = {
        "blocks": n,
        "cuboids": len(cub),
        "cuboids_if_cavities_filled": len(cub) + len(cav_only),
        "cuboid_pct": 100.0 * len(cub) / n,
        "cells_in_cuboids": int(ncell_cub),
        "cells_pct": 100.0 * ncell_cub / len(S.hexes),
        "tiny_blocks": int(sum(1 for s in sizes if s < 10)),
        "small_blocks": int(sum(1 for s in sizes if s < 100)),
        "faces_per_block": {int(k): int(v) for k, v in sorted(nf.items())},
        "min_cells": int(sizes[-1]),
        "inverted_cells": int((sj <= 0).sum()),
        "min_scaled_jacobian": float(sj.min()),
        "mean_scaled_jacobian": float(sj.mean()),
    }
    if verbose:
        print(f"\n===== {title} =====")
        print(f"  {n} blocks, {len(cub)} cuboids ({rep['cuboid_pct']:.0f} %), "
              f"covering {ncell_cub}/{len(S.hexes)} cells "
              f"({rep['cells_pct']:.0f} %)")
        print(f"  faces per block: {rep['faces_per_block']}")
        print(f"  block sizes: {sizes[:8]} ... {sizes[-8:]}")
        print(f"  {rep['tiny_blocks']} blocks < 10 cells, "
              f"{rep['small_blocks']} < 100 cells")
        print(f"  scaled Jacobian: {rep['inverted_cells']} cells <= 0, "
              f"min {sj.min():.4f}, mean {sj.mean():.4f}")
        reasons = Counter(s for s, _d in stat.values() if s != "cuboid")
        if reasons:
            print(f"  non-cuboid reasons: {dict(reasons)}")
        if cav_only:
            print(f"  of those, {len(cav_only)} are non-cuboid ONLY because "
                  f"of an internal cavity ({sum(len(cells[r]) for r in cav_only)} "
                  f"cells) -- a refill would make them cuboids")
        nopp = nbad = 0
        for r, p in pat.items():
            opp, bad = duplicate_surface_faces(p)
            nopp += bool(opp)
            nbad += bool(bad)
        print(f"  blocks touching one surface on two OPPOSITE faces: {nopp} "
              f"(legitimate) / on two ADJACENT faces: {nbad} (defect)")
    return rep, stat, pat, cells


# --------------------------------------------------------------------------
# steps 2-4: the merge operations
# --------------------------------------------------------------------------

@contextmanager
def _temp_merge(S, ra, rb):
    """Apply a merge, evaluate, roll it back. `parent` is a plain array and
    `root` path-compresses, so the rollback has to restore every entry that
    could have been rewritten -- a snapshot is the only safe way."""
    snap = list(S.parent)
    S.parent[rb] = ra
    S._cache.clear()
    try:
        yield
    finally:
        S.parent[:] = snap
        S._cache.clear()


def _excess(patches):
    """How far this block is from being a cuboid, in faces.

    Counting non-cuboid *blocks* is the wrong objective and produced a
    measurably worse structure: merging a 9-faced and a 7-faced block into
    one 12-faced block reduces the count 2 -> 1 and reads as progress, and
    iterating that built 17- and 27-faced monsters while the cuboid
    percentage rose (because the denominator shrank). Summing the excess
    faces instead is monotone -- a merge can only be accepted if the
    structure genuinely gets closer to all-cuboid."""
    if cuboid_status(patches)[0] == "cuboid":
        return 0
    return max(1, len(patches) - 6)


def _score(S, roots):
    """(total excess faces, -cuboid count) over a set of blocks.

    Lexicographic: excess first, cuboid count as the tie-break. A split that
    turns a 10-faced block into a 9-faced one plus a genuine cuboid leaves
    the excess unchanged but is clearly progress, and a pure excess test
    throws it away."""
    exc = cub = 0
    for r in roots:
        p = S.patches(r)
        exc += _excess(p)
        cub += cuboid_status(p)[0] == "cuboid"
    return (exc, -cub)


def _merge_gain(S, a, b, cells):
    """Effect of merging blocks a, b, measured as total excess faces over
    every block the merge can change. Returns (excess_before, excess_after)."""
    ra, rb = S.root(a), S.root(b)
    touched = {ra, rb} | S.neighbours(ra, None) | S.neighbours(rb, None)
    before = sum(_excess(S.patches(r)) for r in touched)
    union = cells[ra] + cells[rb]
    with _temp_merge(S, ra, rb):
        after = 0
        for r in touched:
            if r == rb:
                continue
            c = union if r == ra else cells[r]
            after += _excess(S.patches(root=r, cells=c))
    return before, after


def total_excess(S):
    """Excess faces summed over the whole block structure."""
    return sum(_excess(S.patches(r)) for r in S.cells_of())


def collapse_small_sheets(S, max_faces, verbose=True):
    """Step 2 -- the 3D analogue of `_drop_degenerate_corner_seps`.

    A sheet with very few faces separates almost nothing; removing it from
    the cut set merges the blocks on either side.

    This has to be ATOMIC over the whole sheet. Merging one block pair at a
    time never works: stacking two cuboids leaves each of their four side
    faces split in two, because the side neighbours are still separate
    blocks -- measured, every such pair-wise merge cost exactly +3 excess
    faces and was rejected, so a pair-wise version of this function does
    nothing at all. Dropping the entire sheet merges the side neighbours in
    the same step, and only then is the union a cuboid again."""
    n0 = len(S.cells_of())
    sheet_size = Counter(S.sheet_of.values())
    small = sorted((n, s) for s, n in sheet_size.items() if n <= max_faces)
    if verbose:
        print(f"[clean_blocks] sliver sheets (<= {max_faces} faces): "
              f"{[n for n, _s in small]} of {sorted(sheet_size.values(), reverse=True)}")
    done = 0
    for _n, s in small:
        cells = S.cells_of()
        pairs = set()
        for fk, hs in S.f2h.items():
            if len(hs) != 2 or S.sheet_of.get(fk) != s:
                continue
            ra, rb = S.root(S.blk0[hs[0]]), S.root(S.blk0[hs[1]])
            if ra != rb:
                pairs.add((min(ra, rb), max(ra, rb)))
        if not pairs:
            continue
        before = total_excess(S)
        snap = list(S.parent)
        for a, b in pairs:
            S.merge(a, b)
        after = total_excess(S)
        bad = sum(1 for r in S.cells_of()
                  if duplicate_surface_faces(S.patches(r))[1])
        if after <= before and not bad:
            done += 1
            if verbose:
                print(f"[clean_blocks]   sheet {s} ({_n} faces): "
                      f"{len(pairs)} pairs merged, excess {before} -> {after}")
        else:
            S.parent[:] = snap
            S._cache.clear()
            if verbose:
                print(f"[clean_blocks]   sheet {s} ({_n} faces): rejected "
                      f"(excess {before} -> {after}, {bad} folded blocks)")
    n1 = len(S.cells_of())
    if verbose:
        print(f"[clean_blocks] collapse_small_sheets: {done} sheets removed, "
              f"{n0} -> {n1} blocks")
    return done


def _surface_conflict(S, a, b, cells):
    """Guard from the plan: never merge across a sheet whose two sides carry
    different physical surfaces in a way that would put two ADJACENT faces of
    the merged block on the same surface (a block folding around a wall)."""
    ra, rb = S.root(a), S.root(b)
    union = cells[ra] + cells[rb]
    with _temp_merge(S, ra, rb):
        _opp, bad = duplicate_surface_faces(S.patches(root=ra, cells=union))
    return bool(bad)


def absorb_tiny_blocks(S, min_cells, verbose=True):
    """Step 3 -- blocks below a cell threshold are merged into the neighbour
    they share the largest interface with. Analogue of the stub filter."""
    n0 = len(S.cells_of())
    done = 0
    while True:
        cells = S.cells_of()
        tiny = sorted((r for r, c in cells.items() if len(c) < min_cells),
                      key=lambda r: len(cells[r]))
        hit = False
        for r in tiny:
            if S.root(r) != r:
                continue
            iface = Counter()
            for k, loops in S.patches(r):
                if k[0] == "B":
                    iface[k[1]] += len(loops)
            for nb, _sz in iface.most_common():
                if _surface_conflict(S, r, nb, cells):
                    continue
                bef, aft = _merge_gain(S, r, nb, cells)
                if aft <= bef:
                    S.merge(nb, r)
                    done += 1
                    hit = True
                    break
            if hit:
                break
        if not hit:
            break
    n1 = len(S.cells_of())
    if verbose:
        print(f"[clean_blocks] absorb_tiny_blocks (< {min_cells} cells): "
              f"{done} merges, {n0} -> {n1} blocks")
    return done


def merge_non_cuboids(S, max_cells=None, verbose=True):
    """Step 4b -- a non-cuboid block merged with a neighbour is often a
    cuboid again (the sheet cut a corner off it). Only STRICT reductions of
    total excess faces are accepted, so this can never trade two mediocre
    blocks for one bad one.

    Blocks whose only defect is an internal cavity are skipped: no merge can
    repair a hole, so attempting one just coarsens the structure blindly."""
    n0 = len(S.cells_of())
    done = 0
    while True:
        cells = S.cells_of()
        bad = []
        for r, c in cells.items():
            p = S.patches(r)
            if cuboid_status(p)[0] == "cuboid":
                continue
            tc = touches_cavity(S, p)
            if tc and len(p) - len(tc) == 6:
                continue                     # cavity-only defect, unfixable
            if max_cells is None or len(c) <= max_cells:
                bad.append(r)
        bad.sort(key=lambda r: len(cells[r]))
        hit = False
        for r in bad:
            if S.root(r) != r:
                continue
            iface = Counter()
            for k, loops in S.patches(r):
                if k[0] == "B":
                    iface[k[1]] += len(loops)
            best, gain = None, 0
            for nb, sz in iface.most_common():
                if _surface_conflict(S, r, nb, cells):
                    continue
                before, after = _merge_gain(S, r, nb, cells)
                if before - after > gain:
                    best, gain = nb, before - after
            if best is not None:
                S.merge(best, r)
                done += 1
                hit = True
                break
        if not hit:
            break
    n1 = len(S.cells_of())
    if verbose:
        print(f"[clean_blocks] merge_non_cuboids: {done} merges, "
              f"{n0} -> {n1} blocks")
    return done


def _loop_of(hexes, h, fk):
    c = hexes[h]
    for fc in HF:
        lp = [int(c[i]) for i in fc]
        if frozenset(lp) == fk:
            return lp
    return None


def _face_edges(loop):
    return [(min(loop[k], loop[(k + 1) % 4]), max(loop[k], loop[(k + 1) % 4]))
            for k in range(4)]


def _straight_across(S, fk, g, inside):
    """Rotate twice around edge g, staying inside the block. Same walk as
    `base_complex.sheet_faces`, which is what makes the result a surface
    rather than a chain of faces."""
    for h0 in S.f2h.get(fk, ()):
        if h0 not in inside:
            continue
        f1 = _other_face_with_edge(S, fk, g, h0)
        if f1 is None:
            continue
        h1 = [h for h in S.f2h.get(f1, ()) if h != h0 and h in inside]
        if not h1:
            continue
        f2 = _other_face_with_edge(S, f1, g, h1[0])
        if f2 is not None and f2 != fk:
            return f2
    return None


def _other_face_with_edge(S, fk, g, hi):
    c = S.hexes[hi]
    for fc in HF:
        lp = [int(c[i]) for i in fc]
        k2 = frozenset(lp)
        if k2 == fk:
            continue
        if g in _face_edges(lp):
            return k2
    return None


def extension_surface(S, cells, patch_loops):
    """Continue a block face into the block, as a separating surface.

    The 3D analogue of extending a separatrix from a dangling endpoint
    (`_snap_separatrix_endpoints`). A block with 7 faces is a cuboid with a
    notch: the sheet that would have cut it in two died at a singular edge
    on the notch's border. Restarting the same straight-across propagation
    from that border, restricted to the block's interior, reconstructs the
    missing piece and splits the block."""
    inside = set(cells)
    own = {frozenset(lp) for lp in patch_loops}
    # border edges of the patch: used by exactly one of its faces
    ec = Counter()
    for lp in patch_loops:
        for g in _face_edges(lp):
            ec[g] += 1
    border = {g for g, n in ec.items() if n == 1}
    seeds = set()
    for g in border:
        for h in S.e2h.get(g, ()):
            if h not in inside:
                continue
            for k in S._face_keys(h):
                if k in own or g not in _face_edges(_loop_of(S.hexes, h, k)):
                    continue
                hs = S.f2h[k]
                if len(hs) == 2 and hs[0] in inside and hs[1] in inside:
                    seeds.add(k)
    surf, dq = set(), deque(seeds)
    while dq:
        fk = dq.popleft()
        if fk in surf:
            continue
        surf.add(fk)
        for g in _face_edges(_loop_of(S.hexes, S.f2h[fk][0], fk)):
            if g in S.sing:
                continue
            nxt = _straight_across(S, fk, g, inside)
            if nxt is None or nxt in surf:
                continue
            hs = S.f2h[nxt]
            if len(hs) == 2 and hs[0] in inside and hs[1] in inside:
                dq.append(nxt)
    return surf


def split_non_cuboid(S, verbose=True):
    """Step 4b -- split a block whose face count is too high along the
    surface that continues one of its faces. Accepted only when the total
    excess strictly drops, so a split that merely trades one bad block for
    two can never be applied."""
    n0 = len(S.cells_of())
    done = 0
    for _round in range(40):
        cells = S.cells_of()
        bad = []
        for r, c in cells.items():
            p = S.patches(r)
            if cuboid_status(p)[0] == "cuboid":
                continue
            tc = touches_cavity(S, p)
            if tc and len(p) - len(tc) == 6:
                continue                      # cavity-only, unsplittable
            bad.append(r)
        bad.sort(key=lambda r: -len(S.patches(r)))
        hit = False
        for r in bad:
            p = S.patches(r)
            touched = {r} | {k[1] for k, _l in p if k[0] == "B"}
            before = _score(S, touched)
            best = None
            # try continuing each face, smallest first: a sliver face is the
            # one most likely to be the stump of a sheet that died
            for _sz, pi in sorted((len(l), i) for i, (_k, l) in enumerate(p)):
                surf = extension_surface(S, cells[r], p[pi][1])
                if not surf:
                    continue
                snap = (list(S.parent), S.blk0.copy(), set(S.cut),
                        set(S.extra_cut), S.n_base)
                roots = S.split(r, surf)
                if roots is not None:
                    after = _score(S, set(roots) | (touched - {r}))
                    folded = any(duplicate_surface_faces(S.patches(x))[1]
                                 for x in roots)
                    if after < before and not folded:
                        if verbose:
                            print(f"[clean_blocks]   split b{r} "
                                  f"({len(cells[r])} cells, {len(p)} faces) "
                                  f"along {len(surf)} faces -> {len(roots)} "
                                  f"blocks, (excess, -cuboids) {before} -> "
                                  f"{after}")
                        best = True
                        break
                (S.parent, S.blk0, S.cut, S.extra_cut, S.n_base) = snap
                S._cache.clear()
            if best:
                done += 1
                hit = True
                break
        if not hit:
            break
    n1 = len(S.cells_of())
    if verbose:
        print(f"[clean_blocks] split_non_cuboid: {done} splits, "
              f"{n0} -> {n1} blocks")
    return done


# --------------------------------------------------------------------------
# step 5: validation
# --------------------------------------------------------------------------

class HexBlockValidator:
    """3D counterpart of `dp3d.field.quad_partition_validator.
    QuadPartitionValidator`: same three parts -- hard `is_valid()` checks,
    a soft `quality_score()`, and human-readable `diagnostics()`."""

    def __init__(self, S):
        self.S = S
        self.cells = S.cells_of()
        self.pat = {r: S.patches(r) for r in self.cells}
        self.sj = np.array([ovm_io.scaled_jacobian(S.P, c) for c in S.hexes])

    # -- hard ---------------------------------------------------------------
    def is_valid(self):
        return not [m for m in self.hard_failures()]

    def hard_failures(self):
        out = []
        covered = sum(len(c) for c in self.cells.values())
        if covered != len(self.S.hexes):
            out.append(f"block complex covers {covered}/{len(self.S.hexes)} cells")
        n_inv = int((self.sj <= 0).sum())
        if n_inv:
            out.append(f"{n_inv} inverted cells (scaled Jacobian <= 0)")
        cav = getattr(self.S, "cavities", [])
        if cav:
            out.append(f"{len(cav)} internal cavities in the hex mesh "
                       f"({', '.join(str(n) for _s, n in cav)} quads) -- "
                       f"inherited from HexEx, not created here")
        nbad = 0
        for r, p in self.pat.items():
            _opp, bad = duplicate_surface_faces(p)
            if bad:
                nbad += 1
        if nbad:
            out.append(f"{nbad} blocks with two ADJACENT faces on the same "
                       f"surface")
        mixed = 0
        for r, p in self.pat.items():
            for k, loops in p:
                if k[0] != "P":
                    continue
                labs = {self.S.surf_of[frozenset(lp)] for lp in loops}
                if len(labs) > 1:
                    mixed += 1
        if mixed:
            out.append(f"{mixed} block faces spanning two physical surfaces")
        return out

    # -- soft ---------------------------------------------------------------
    def quality_score(self):
        n = len(self.cells)
        cub = sum(cuboid_status(p)[0] == "cuboid" for p in self.pat.values())
        sizes = np.array([len(c) for c in self.cells.values()])
        return {
            "blocks": n,
            "cuboid_fraction": cub / n,
            "tiny_blocks": int((sizes < 10).sum()),
            "min_cells_per_block": int(sizes.min()),
            "median_cells_per_block": float(np.median(sizes)),
            "min_scaled_jacobian": float(self.sj.min()),
            "mean_scaled_jacobian": float(self.sj.mean()),
        }

    def boundary_hausdorff(self, labeller):
        """Distance of the block-complex boundary to the input surface.
        With cut-set-only simplification this is invariant by construction;
        it is measured anyway so that any future geometric step cannot
        silently deform the domain."""
        if labeller is None:
            return None
        Q = self.S.P[np.asarray(self.S.bnd_loops)].mean(1)
        labeller.label(Q)
        return labeller.max_dist

    def diagnostics(self):
        lines = []
        fails = self.hard_failures()
        lines.append("VALID" if not fails else "INVALID")
        lines += [f"  ! {m}" for m in fails]
        for k, v in self.quality_score().items():
            lines.append(f"  {k}: {v}")
        reasons = Counter(cuboid_status(p)[0] for p in self.pat.values())
        lines.append(f"  block face counts: "
                     f"{dict(sorted(Counter(len(p) for p in self.pat.values()).items()))}")
        lines.append(f"  status: {dict(reasons)}")
        return "\n".join(lines)


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------

def postprocess(ovm_path, input_vtk=None, names=None, sheet_frac=0.30,
                min_cells=10, fix_non_cuboid=True, fill=True,
                collapse_rounds=0, untangle_mesh=False, verbose=True):
    """`sheet_frac` is a fraction of the MEDIAN sheet size. The plan asked
    for "~5 %, i.e. the 27/54/154-face sheets", but on v9 the median is 533,
    so 5 % is 26 and catches none of them -- the two numbers in the plan
    contradict each other. 30 % (=160) is what actually selects the three
    sheets the plan names."""
    """Mirrors `StreamlinePostProcessor`: diagnose, clean, validate."""
    labeller = None
    if input_vtk is not None:
        P, tri, tid = read_input_surface(input_vtk)
        labeller = SurfaceLabeller(P, tri, tid, names)
    S = BlockStructure(ovm_path, labeller, fill=fill,
                       collapse_rounds=collapse_rounds,
                       untangle_mesh=untangle_mesh)
    print(f"[clean_blocks] surface labels from: {S.label_source}")
    # step 0: fix the classification before any geometry/topology is touched
    label_cavities(S, verbose)
    remove_label_islands(S, verbose=verbose)
    print(f"[clean_blocks] surfaces: "
          f"{ {S.surf_names.get(s, s): n for s, n in sorted(Counter(S.surf_of.values()).items())} }")
    before, _st, _p, _c = report_structure(S, "BEFORE")

    med = np.median(list(Counter(S.sheet_of.values()).values()))
    collapse_small_sheets(S, max_faces=max(1, int(sheet_frac * med)),
                          verbose=verbose)
    absorb_tiny_blocks(S, min_cells=min_cells, verbose=verbose)
    if fix_non_cuboid:
        merge_non_cuboids(S, verbose=verbose)
        split_non_cuboid(S, verbose=verbose)

    after, _st, _p, _c = report_structure(S, "AFTER")
    v = HexBlockValidator(S)
    print("\n----- HexBlockValidator -----")
    print(v.diagnostics())
    hd = v.boundary_hausdorff(labeller)
    if hd is not None:
        print(f"  boundary Hausdorff to input surface: {hd:.5f}")
    return S, before, after, v


def block_edge_curves(S):
    """The 1D skeleton of the block structure: the curves where two faces of
    the same block meet.

    A mesh edge is on a block edge iff some block sees it on the border
    between two of its own faces. Block membership alone is not enough to
    find these -- an edge running along the corner between two boundary
    faces of a single block is incident to exactly one block, like every
    other edge of that block's surface -- which is why this needs the face
    partition, and hence the surface labels.

    Returns (segments, curve id per segment). Curves are split at block
    corners, i.e. wherever a vertex is not simply passed through, so each
    curve is one edge of the block topology -- the entity the conforming
    division MILP assigns a division count to later."""
    on_edge = set()
    for r in S.cells_of():
        e2p = defaultdict(set)
        for pi, (_k, loops) in enumerate(S.patches(r)):
            for lp in loops:
                for g in _face_edges(lp):
                    e2p[g].add(pi)
        for g, ps in e2p.items():
            if len(ps) > 1:
                on_edge.add(g)
    adj = defaultdict(list)
    for u, v in on_edge:
        adj[u].append(v)
        adj[v].append(u)
    # a corner is any vertex that is not a plain pass-through
    corner = {v for v, ns in adj.items() if len(ns) != 2}
    seen, segs, cid = set(), [], []
    curve = 0

    def walk(a, b):
        nonlocal curve
        path = [a, b]
        seen.add((min(a, b), max(a, b)))
        while path[-1] not in corner:
            nxt = [w for w in adj[path[-1]] if w != path[-2]]
            if len(nxt) != 1:
                break
            g = (min(path[-1], nxt[0]), max(path[-1], nxt[0]))
            if g in seen:
                break
            seen.add(g)
            path.append(nxt[0])
        for k in range(len(path) - 1):
            segs.append((path[k], path[k + 1]))
            cid.append(curve)
        curve += 1

    for c in sorted(corner):
        for w in adj[c]:
            if (min(c, w), max(c, w)) not in seen:
                walk(c, w)
    for u, v in sorted(on_edge):          # closed loops without any corner
        if (u, v) not in seen:
            walk(u, v)
    return np.array(segs, np.int64), np.array(cid, np.int64), curve


def edge_kink_stats(P, segs, cid):
    """Turn angle at each interior vertex of every block-edge curve."""
    byc = defaultdict(list)
    for (a, b), c in zip(segs, cid):
        byc[int(c)].append((int(a), int(b)))
    out = {}
    for c, es in byc.items():
        path = _polyline(es)
        A = []
        for i in range(1, len(path) - 1):
            u, w = P[path[i]] - P[path[i - 1]], P[path[i + 1]] - P[path[i]]
            nu, nw = np.linalg.norm(u), np.linalg.norm(w)
            if nu * nw:
                A.append(np.degrees(np.arccos(
                    np.clip(float(np.dot(u, w)) / (nu * nw), -1, 1))))
        out[c] = (path, np.array(A))
    return out


def _polyline(es):
    adj = defaultdict(list)
    for a, b in es:
        adj[a].append(b)
        adj[b].append(a)
    ends = [v for v, n in adj.items() if len(n) == 1]
    start = ends[0] if ends else es[0][0]
    path, prev, cur = [start], None, start
    while True:
        nxt = [w for w in adj[cur] if w != prev]
        if not nxt:
            break
        prev, cur = cur, nxt[0]
        path.append(cur)
        if cur == start:
            break
    return path


def smooth_block_edges(P, hexes, f2h, segs, cid, labeller, iters=30,
                       alpha=0.5, skip=(), verbose=True):
    """Straighten the block-edge curves by moving the mesh vertices on them.

    Curve endpoints (block corners) are pinned; interior curve vertices move
    toward the midpoint of their two curve neighbours. Vertices on the domain
    boundary are re-projected onto the input surface after every step, so
    they slide along the geometry instead of leaving it. Every step is
    reverted if it inverts a cell.

    `skip` excludes curves this cannot legitimately fix. On v9 that is the 14
    `bl_interface_hub | ogrid_interface` curves: they carry 45 of the 63 kinks above 30
    degrees, but the hex mesh is not aligned to that feature curve at all --
    NO hex boundary vertex lies within 0.010 of it, the staircases run 1-3
    cells beside it (median 0.074). There is no smooth edge chain to move
    onto, so smoothing them would only shrink the staircase amplitude while
    dragging the surface mesh off the geometry: cosmetics over a real
    alignment defect, which is why it is opt-in and not the default."""
    byc = defaultdict(list)
    for (a, b), c in zip(segs, cid):
        byc[int(c)].append((int(a), int(b)))
    deg = Counter()
    for a, b in segs:
        deg[int(a)] += 1
        deg[int(b)] += 1
    pinned = {v for v, d in deg.items() if d != 2}
    bnd_v = {v for fk, hs in f2h.items() if len(hs) == 1 for v in fk}
    v2c = defaultdict(list)
    for i, c in enumerate(hexes):
        for v in c:
            v2c[int(v)].append(i)

    moving = []
    for c, es in byc.items():
        if c in skip:
            continue
        path = _polyline(es)
        for i in range(1, len(path) - 1):
            if path[i] not in pinned:
                moving.append((path[i], path[i - 1], path[i + 1]))
    if not moving:
        return P, 0
    idx = np.array([m[0] for m in moving])
    lo = np.array([m[1] for m in moving])
    hi = np.array([m[2] for m in moving])
    on_b = np.array([v in bnd_v for v in idx])
    watch = sorted({c for v in idx for c in v2c[int(v)]})
    P = P.copy()
    ok0 = scaled_jacobians(P, hexes[watch]).min()
    # Per-VERTEX acceptance with backtracking. A single global damped step
    # is useless here: the touched cells start at a minimum scaled Jacobian
    # of only ~0.10, so one sweep at alpha=0.5 inverts a cell and every one
    # of 30 iterations gets rejected. Moving one vertex at a time, and
    # halving its step until its own cells stay healthy, converges instead.
    cell_of = {v: np.array(v2c[int(v)]) for v in idx}
    # Quality floor measured against the ORIGINAL mesh, not the running one:
    # a per-step "no worse than 80 % of current" rule compounds over the
    # sweeps and drove the touched cells from 0.0986 down to 0.0047.
    floor = {v: max(0.0, 0.8 * float(scaled_jacobians(P, hexes[cell_of[v]]).min()))
             for v in map(int, idx)}

    def kink_score(Q):
        # only over the curves this operation may actually move: the skipped
        # hub|blade staircases contribute 45 of the 63 kinks above 30 deg and
        # cannot improve, so including them makes every sweep look useless
        ks = edge_kink_stats(Q, segs, cid)
        A = np.concatenate([v[1] for c, v in ks.items()
                            if c not in skip and len(v[1])])
        return (int((A > 30).sum()), float(np.percentile(A, 95)),
                float(np.median(A)))

    score = kink_score(P)
    if verbose:
        print(f"[clean_blocks]   start (>30 deg, p95, median) = {score}")
    nmoved = 0
    for _it in range(iters):
        snapshot = P.copy()
        moved_this = 0
        for j in range(len(idx)):
            v = int(idx[j])
            cl = cell_of[v]
            tgt = 0.5 * (P[int(lo[j])] + P[int(hi[j])])
            step = alpha
            for _bt in range(4):
                q = P[v] + step * (tgt - P[v])
                if on_b[j] and labeller is not None:
                    q = project_to_surface(labeller, q[None, :])[0]
                old = P[v].copy()
                P[v] = q
                if float(scaled_jacobians(P, hexes[cl]).min()) > floor[v]:
                    moved_this += 1
                    break
                P[v] = old
                step *= 0.5
        new = kink_score(P)
        # accept the sweep only if the TAILS improve; the mean is not the
        # thing that matters and optimising it made the tails worse
        if new >= score:
            P = snapshot
            break
        score = new
        nmoved += moved_this
        if not moved_this:
            break
    if verbose:
        print(f"[clean_blocks]   end   (>30 deg, p95, median) = {score}")
    ok1 = scaled_jacobians(P, hexes[watch]).min()
    if verbose:
        print(f"[clean_blocks] smooth_block_edges: {len(idx)} vertices over "
              f"{len(byc) - len(set(skip))} curves, {nmoved} accepted moves; "
              f"min scaled Jacobian on the {len(watch)} touched cells "
              f"{ok0:.4f} -> {ok1:.4f}")
    return P, len(idx)


def _ogrid_interface_curves(S, segs, cid):
    """Curves where the two artificial cut surfaces meet: the interface to
    the removed hub/shroud boundary layer against the interface to the
    removed blade O-grid. Neither is a wall -- see the note in
    `tet_prep_v5.NAMES`. Excluded from smoothing, see `smooth_block_edges`."""
    e2k = defaultdict(set)
    for r in S.cells_of():
        for k, loops in S.patches(r):
            for lp in loops:
                for g in _face_edges(lp):
                    e2k[g].add(k)
    byc = defaultdict(set)
    for (a, b), c in zip(segs, cid):
        for k in e2k[(min(int(a), int(b)), max(int(a), int(b)))]:
            if k[0] == "P":
                byc[int(c)].add(S.surf_names.get(k[1]))
    return {c for c, n in byc.items()
            if {"bl_interface_hub", "ogrid_interface"} <= n}


def write_block_edges_vtk(path, P, segs, cid, title):
    used = np.unique(segs)
    remap = {int(v): i for i, v in enumerate(used)}
    with open(path, "w") as f:
        f.write(f"# vtk DataFile Version 2.0\n{title}\nASCII\n"
                f"DATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {len(used)} float\n")
        for v in used:
            f.write(f"{P[v][0]} {P[v][1]} {P[v][2]}\n")
        f.write(f"CELLS {len(segs)} {3 * len(segs)}\n")
        for a, b in segs:
            f.write(f"2 {remap[int(a)]} {remap[int(b)]}\n")
        f.write(f"CELL_TYPES {len(segs)}\n")
        for _ in segs:
            f.write("3\n")                                     # VTK_LINE
        f.write(f"CELL_DATA {len(segs)}\nSCALARS block_edge_id int 1\n"
                f"LOOKUP_TABLE default\n")
        for c in cid:
            f.write(f"{int(c)}\n")


def write_blocks(S, path, title="cleaned blocks"):
    import export_vtk as ev
    cells = S.cells_of()
    order = {r: i for i, r in enumerate(sorted(cells))}
    bid = np.zeros(len(S.hexes), int)
    nf = np.zeros(len(S.hexes), int)
    for r, c in cells.items():
        k = len(S.patches(r))
        for h in c:
            bid[h] = order[r]
            nf[h] = k
    segs, cid, ncurve = block_edge_curves(S)
    print(f"[clean_blocks] block edges: {len(segs)} segments in "
          f"{ncurve} curves")
    sj = scaled_jacobians(S.P, S.hexes)

    # Compact once, at the very end, and remap the block edges with the same
    # table -- they index the same point array.
    n0 = len(S.P)
    used = np.unique(S.hexes)
    remap = np.full(n0, -1, np.int64)
    remap[used] = np.arange(len(used))
    P, H = S.P[used], remap[S.hexes]
    segs = remap[segs] if len(segs) else segs
    if len(used) != n0:
        print(f"[clean_blocks] compacted points {n0} -> {len(used)} "
              f"({n0 - len(used)} unreferenced, left behind by the collapse)")

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    ev.write_vtk(str(path), P, H, [12] * len(H), bid,
                 "block_id", f"{len(cells)} blocks -- {title}")
    ev.write_vtk(str(path).replace(".vtk", "_nfaces.vtk"), P, H,
                 [12] * len(H), nf, "n_block_faces",
                 "faces per block (6 = cuboid)")
    ev.write_vtk(str(path).replace(".vtk", "_quality.vtk"), P, H,
                 [12] * len(H), np.round(sj * 1000).astype(int),
                 "scaled_jacobian_x1000",
                 f"min {sj.min():.4f}, mean {sj.mean():.4f}, "
                 f"{int((sj <= 0).sum())} inverted")
    print(f"wrote {str(path).replace('.vtk', '_quality.vtk')}")
    ovm_io.write_hex_msh(str(path).replace(".vtk", ".msh"), P, H, bid,
                         lines=segs, line_tags=cid + 1)
    write_block_edges_vtk(str(path).replace(".vtk", "_edges.vtk"), P, segs,
                          cid, f"block edges -- {ncurve} curves, {title}")
    print(f"wrote {str(path).replace('.vtk', '_edges.vtk')}")
    return bid


if __name__ == "__main__":
    import argparse
    import tet_prep_v5 as v5
    ap = argparse.ArgumentParser()
    ap.add_argument("ovm", nargs="?", default=str(OUT / "T1_9_hex_v9.ovm"))
    ap.add_argument("--input-vtk",
                    default=str(REPO / "data" / "T1_9" / "T1_9_tet_v5.vtk"))
    ap.add_argument("--no-input-vtk", action="store_true",
                    help="label the boundary by feature angle instead "
                         "(for meshes without a labelled input, e.g. cylinder)")
    ap.add_argument("--sheet-frac", type=float, default=0.30)
    ap.add_argument("--collapse-rounds", type=int, default=0,
                    help="greedy mesh-level sheet collapses (Gao et al.); "
                         "changes the hex mesh, slow (~15 min per round)")
    ap.add_argument("--smooth-edges", action="store_true",
                    help="straighten block-edge curves (with --edges-only); "
                         "skips the hub|blade staircases, which no amount of "
                         "smoothing can legitimately fix")
    ap.add_argument("--edges-only", metavar="BLOCKS_VTK",
                    help="do not run the pipeline: read an existing blocks "
                         "VTK, recompute its block edges and rewrite the .msh "
                         "with them as 1D elements plus a _edges.vtk")
    ap.add_argument("--untangle", action="store_true",
                    help="repair inverted cells by local smoothing (Gao's "
                         "quality-repair half); moves interior vertices, and "
                         "boundary vertices only along the input surface")
    ap.add_argument("--no-fill", action="store_true",
                    help="leave internal cavities open (they are refilled by "
                         "default; this is the only step that adds cells)")
    ap.add_argument("--min-cells", type=int, default=10)
    ap.add_argument("--out", default=str(OUT / "deliverable"
                                         / "T1_9_blocks_v9_clean.vtk"))
    a = ap.parse_args()
    if a.edges_only:
        import meshio
        m = meshio.read(a.edges_only)
        hx = np.vstack([b.data for b in m.cells if b.type == "hexahedron"])
        bid = np.concatenate([np.asarray(d).ravel() for b, d in
                              zip(m.cells, m.cell_data["block_id"])
                              if b.type == "hexahedron"])
        P0, tri0, tid0 = read_input_surface(a.input_vtk)
        lb = SurfaceLabeller(P0, tri0, tid0, v5.NAMES)
        St = BlockStructure.from_arrays(m.points, hx, bid, lb)
        segs, cid, ncurve = block_edge_curves(St)
        if a.smooth_edges:
            ks = edge_kink_stats(St.P, segs, cid)
            skip = _ogrid_interface_curves(St, segs, cid)
            allA = np.concatenate([v[1] for v in ks.values() if len(v[1])])
            print(f"[clean_blocks] kinks before: median {np.median(allA):.1f} "
                  f"deg, p95 {np.percentile(allA, 95):.1f}, "
                  f"{int((allA > 30).sum())} above 30 deg")
            St.P, _n = smooth_block_edges(St.P, St.hexes, St.f2h, segs, cid,
                                          lb, skip=skip)
            ks = edge_kink_stats(St.P, segs, cid)
            allB = np.concatenate([v[1] for v in ks.values() if len(v[1])])
            print(f"[clean_blocks] kinks after:  median {np.median(allB):.1f} "
                  f"deg, p95 {np.percentile(allB, 95):.1f}, "
                  f"{int((allB > 30).sum())} above 30 deg")
            import export_vtk as ev
            ev.write_vtk(a.edges_only, St.P, St.hexes, [12] * len(St.hexes),
                         bid, "block_id", f"{len(St.cells_of())} blocks")
        print(f"[clean_blocks] {len(St.cells_of())} blocks, {len(segs)} block-"
              f"edge segments in {ncurve} curves")
        ovm_io.write_hex_msh(a.edges_only.replace(".vtk", ".msh"), St.P,
                             St.hexes, bid, lines=segs, line_tags=cid + 1)
        write_block_edges_vtk(a.edges_only.replace(".vtk", "_edges.vtk"),
                              St.P, segs, cid,
                              f"block edges -- {ncurve} curves")
        print(f"wrote {a.edges_only.replace('.vtk', '.msh')} (with 1D block "
              f"edges) and {a.edges_only.replace('.vtk', '_edges.vtk')}")
        raise SystemExit(0)
    S, before, after, v = postprocess(
        a.ovm, None if a.no_input_vtk else a.input_vtk, v5.NAMES,
        a.sheet_frac, a.min_cells, fill=not a.no_fill,
        collapse_rounds=a.collapse_rounds, untangle_mesh=a.untangle)
    write_blocks(S, a.out)
    print("\nbefore/after:")
    for k in before:
        print(f"  {k:22s} {before[k]}  ->  {after[k]}")
