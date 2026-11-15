"""beam_collapse on the base-complex level: coarsen first, then the unchanged beam search.

The cell-level beam search (beam_collapse.py) evaluates every sheet of the fine hex mesh. On a fine mesh
(canada diffuser, n = 100k: 105837 cells, ~500 sheets per state) almost every sheet is just an inner
cell layer of a block: collapsing it changes nothing structural, and each evaluation runs over all
cells (~2 h per search depth on 8 cores). A block layer that is k cells thick also needs k collapse
steps, i.e. k search depths, before it can disappear.

Here the mesh is first reduced STRUCTURE-NEUTRALLY: a sheet is collapsed if every block it crosses keeps
at least one more cell layer in that direction (another sheet crosses the block without sharing a cell
with it). Per round many such sheets go at once, at most one per block; afterwards the structure stats
(blocks, cuboids, excess faces) are verified to be unchanged -- a batch that changes them is bisected.
Moved boundary vertices are projected back onto the input surface. The result has the same block
structure with every block one cell thick per direction (a few hundred cells), and the unchanged beam
search runs on it: every candidate sheet is then a structural one, and one step removes a whole block
layer.

Guard added for both modes: every input feature curve must stay on a block edge (--feature-tol, default
0.05, never below the start mesh's value).

Modes: default = coarsen + beam on the coarse mesh; --no-coarsen = cell-level beam; --macro = beam on the
fine mesh with block-layer-stack moves. The vectorised fast_stats.structure_stats replaces
_structure_stats in all modes (identical results, ~4x faster; --slow-stats for the original).

usage: same arguments as beam_collapse.py (+ --no-coarsen | --macro, --feature-tol T, --slow-stats)
"""
import sys, time
from collections import defaultdict
from pathlib import Path

HEX = Path("/opt/stack/domain_partition_3D/experimentell/hex3d_algohex")
sys.path.insert(0, str(HEX.parent.parent)); sys.path.insert(0, str(HEX))
import numpy as np                                            # noqa: E402
import clean_blocks as cb                                     # noqa: E402
import beam_collapse as bm                                    # noqa: E402
import base_complex as bc                                     # noqa: E402
import block_faces as bfx                                     # noqa: E402


def _blocks(P, hexes):
    import contextlib, io
    with contextlib.redirect_stdout(io.StringIO()):
        f2h, e2h = bc.build_topology(hexes)
        sing = bc.singular_edges(hexes, P, f2h, e2h)
        sheet_of, _n = bfx.label_sheets(hexes, f2h, sing)
        blocks = bc.blocks_from_cut(hexes, f2h, set(sheet_of))
    blk = np.empty(len(hexes), int)
    for i, b in enumerate(blocks):
        blk[list(b)] = i
    return blk, f2h


def _sig(st):
    return (st["blocks"], st["cuboids"], st["excess"])


def _hex_sheets(hexes):
    """sheet id of each of the three parallel edge groups of every hex"""
    sheets = cb.mesh_sheets(hexes)
    sid = {}
    for k, es in enumerate(sheets):
        for e in es:
            sid[e] = k
    hs = np.empty((len(hexes), 3), int)
    for h, c in enumerate(hexes):
        for g, grp in enumerate(cb._PARALLEL):
            i, j = grp[0]
            hs[h, g] = sid[(min(int(c[i]), int(c[j])), max(int(c[i]), int(c[j])))]
    return sheets, hs


def _neutral_batch(hexes, blk):
    """sheets whose collapse keeps >= 1 cell layer in every block they cross; at most one per block"""
    sheets, hs = _hex_sheets(hexes)
    in_blk = defaultdict(lambda: defaultdict(set))       # block -> sheet -> hexes
    for h in range(len(hexes)):
        for s in hs[h]:
            in_blk[blk[h]][int(s)].add(h)
    crosses = defaultdict(set)
    for b, d in in_blk.items():
        for s in d:
            crosses[s].add(b)
    def keeps_layer(s, b):
        mine = in_blk[b][s]
        return any(t != s and not (cells & mine) for t, cells in in_blk[b].items())
    cand = [s for s in crosses if all(keeps_layer(s, b) for b in crosses[s])]
    cand.sort(key=lambda s: len(sheets[s]))
    used, batch = set(), []
    for s in cand:
        if crosses[s] & used:
            continue
        used |= crosses[s]
        batch.append(s)
    return sheets, batch


def _on_features(Q, tol=None):
    """mask: point lies on an input feature curve (point-segment distance below 2e-4 of the model size;
    AlgoHex places boundary vertices on the feature curves only up to ~1e-4: tistos 67 vertices at 1e-6,
    180 at 1e-3, 184 at 1e-2)"""
    tol = _FEAT["on_tol"] if tol is None else tol
    from scipy.spatial import cKDTree
    if _FEAT.get("segs") is None:
        return np.zeros(len(Q), bool)
    A, B, tree = _FEAT["segs"]
    k = min(8, len(A))
    _d, j = tree.query(Q, k=k)
    j = np.atleast_2d(j)
    a, b = A[j], B[j]
    ab = b - a
    t = np.clip(((Q[:, None] - a) * ab).sum(-1) / np.maximum((ab * ab).sum(-1), 1e-30), 0, 1)
    d = np.linalg.norm(Q[:, None] - (a + t[..., None] * ab), axis=-1).min(1)
    return d < tol


def _collapse(P, hexes, sheets, batch, labeller):
    """collapse the sheets of the batch at once. Position of a welded vertex class: a member on an input
    feature curve keeps its place (block edges must stay on the geometric edges), else boundary members
    are averaged and projected onto the surface, else the plain average."""
    edges = [e for s in batch for e in sheets[s]]
    Q, H, _ = cb.collapse_sheet(P, hexes, edges)          # topology + plain averages
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
    for v in list(vp):
        grp[find(v)].append(v)
    f2h0, _ = bc.build_topology(hexes)
    bnd0 = {v for fk, hs_ in f2h0.items() if len(hs_) == 1 for v in fk}
    members = sorted({v for vs in grp.values() for v in vs})
    feat = dict(zip(members, _on_features(P[members]))) if members else {}
    Q = Q.copy()
    proj = []
    for r, vs in grp.items():
        if len(vs) < 2:
            continue
        fv = [v for v in vs if feat.get(v)]
        bv = [v for v in vs if v in bnd0]
        if fv:
            Q[vs] = P[fv[0]]
        elif bv:
            Q[vs] = P[bv].mean(0)
            proj += vs
    if proj and labeller is not None:
        proj = sorted(set(proj))
        Q[proj] = cb.project_to_surface(labeller, Q[proj])
    return Q, H


def coarsen(P, hexes, labeller, verbose=True):
    t0 = time.time()
    ref = cb._structure_stats(P, hexes, labeller)
    if verbose:
        print(f"[coarse] start {ref['cells']} cells, {ref['blocks']} blocks, {ref['cuboids']} cuboids, "
              f"excess {ref['excess']} ({time.time() - t0:.0f}s)", flush=True)
    rnd = 0
    while True:
        rnd += 1
        blk, _ = _blocks(P, hexes)
        sheets, batch = _neutral_batch(hexes, blk)
        if not batch:
            break
        todo, done = [batch], 0
        while todo:
            part = todo.pop()
            Q, H = _collapse(P, hexes, sheets, part, labeller)
            st = cb._structure_stats(Q, H, labeller) if len(H) else None
            if st is not None and _sig(st) == _sig(ref):
                P, hexes, done = Q, H, done + len(part)
                break          # sheet ids refer to the old mesh: recompute before the next part
            if len(part) > 1:
                todo += [part[len(part) // 2:], part[:len(part) // 2]]
        if verbose:
            print(f"[coarse] round {rnd}: {len(batch)} neutral candidates, {done} sheets collapsed -> "
                  f"{len(hexes)} cells ({time.time() - t0:.0f}s)", flush=True)
        if not done:
            break
    st = cb._structure_stats(P, hexes, labeller)
    if verbose:
        print(f"[coarse] result {st} ({time.time() - t0:.0f}s)", flush=True)
    return P, hexes


# ---- feature guard: every input feature curve must stay on a block edge --------------------------------
# On the coarse mesh the Hausdorff guard (boundary face centroids) is too weak: cells are large, a collapse
# that drops a block edge from a geometric edge hardly moves the face centroids (tistos_s1_0004: 10 blocks,
# 33 % of the feature curves > 0.05 from any block edge; cell level: 12 blocks, max 0.036).
_FEAT = {"pts": None, "tol": None, "labeller": None, "segs": None, "on_tol": 1e-6}


def set_feature_curves(input_vtk, step=None):
    import meshio
    m = meshio.read(input_vtk)
    X = np.asarray(m.points, float)
    E = np.vstack([c.data for c in m.cells if c.type == "line"])
    s = np.linspace(0, 1, 5)[:, None, None]
    _FEAT["pts"] = (X[E[:, 0]][None] * (1 - s) + X[E[:, 1]][None] * s).reshape(-1, 3)
    from scipy.spatial import cKDTree
    _FEAT["segs"] = (X[E[:, 0]], X[E[:, 1]], cKDTree(0.5 * (X[E[:, 0]] + X[E[:, 1]])))
    _FEAT["on_tol"] = 2e-4 * float(np.linalg.norm(X.max(0) - X.min(0)))


def feature_dist(P, hexes):
    """max distance of the input feature curve samples to the block edges of the mesh. A block edge is a
    mesh edge where the boundary of one block changes patch (other block or other surface label)."""
    from scipy.spatial import cKDTree
    if _FEAT["pts"] is None or not len(hexes):
        return 0.0
    blk, f2h = _blocks(P, hexes)
    edges = set()
    # patch key of every block-boundary face, per block
    ekeys = defaultdict(set)
    # surface label of every boundary face (same labeller as _structure_stats): a label change along the
    # boundary is a block edge too -- without it the raw mesh measured 0.73 instead of ~0.03
    bfaces = [(fk, hs_[0]) for fk, hs_ in f2h.items() if len(hs_) == 1]
    surf = {}
    if bfaces and _FEAT.get("labeller") is not None:
        cen = np.array([P[list(fk)].mean(0) for fk, _h in bfaces])
        ids = _FEAT["labeller"].label(cen)
        surf = {fk: int(i) for (fk, _h), i in zip(bfaces, ids)}
    for fk, hs_ in f2h.items():
        if len(hs_) == 2 and blk[hs_[0]] == blk[hs_[1]]:
            continue
        for h in hs_:
            other = [x for x in hs_ if x != h]
            key = ("B", int(blk[other[0]])) if other else ("P", surf.get(fk, -1))
            c = hexes[h]
            for fc in cb.HF:
                lp = [int(c[i]) for i in fc]
                if frozenset(lp) == fk:
                    for a_, b_ in zip(lp, lp[1:] + lp[:1]):
                        ekeys[(int(blk[h]), min(a_, b_), max(a_, b_))].add(key)
                    break
    for (b, u, v), ks in ekeys.items():
        if len(ks) > 1:
            edges.add((u, v))
    # boundary faces carry surface labels: a label change along the boundary is a block edge as well
    if not edges:
        return np.inf
    E = np.array(sorted(edges))
    s = np.linspace(0, 1, 9)[:, None, None]
    B = (P[E[:, 0]][None] * (1 - s) + P[E[:, 1]][None] * s).reshape(-1, 3)
    d, _ = cKDTree(B).query(_FEAT["pts"])
    return float(d.max())


_stats0 = cb._structure_stats


def _stats_with_feat(P, hexes, labeller):
    st = _stats0(P, hexes, labeller)
    st["feat"] = feature_dist(P, hexes)
    return st


_ok0 = bm._ok


def _ok_feat(st, bar, ceiling, guard, sj_floor):
    if st is not None and _FEAT["tol"] is not None and st.get("feat", 0.0) > _FEAT["tol"] + 1e-9:
        return False
    return _ok0(st, bar, ceiling, guard, sj_floor)


# ---- macro moves on the fine mesh ---------------------------------------------------------------------
def macro_sheets(hexes):
    """candidate moves = block layer stacks: sheets that cross the same set of blocks and are parallel
    inside them (share no cell) are collapsed together. A block layer k cells thick disappears in ONE
    move (cell level: k search depths, k-1 of them structure-neutral); geometry stays fine-mesh exact."""
    import fast_stats as fs
    sheets = _mesh_sheets0(hexes)
    if not len(hexes):
        return sheets
    T = fs.topology(hexes)
    cut = fs.separatrix_faces(T, fs.singular(T))
    _nb, blk = fs.blocks(T, cut)
    sid = {}
    for k, es in enumerate(sheets):
        for e in es:
            sid[e] = k
    cells = defaultdict(set)
    for h, c in enumerate(hexes):
        for grp in cb._PARALLEL:
            i, j = grp[0]
            cells[sid[(min(int(c[i]), int(c[j])), max(int(c[i]), int(c[j])))]].add(h)
    bykey = defaultdict(list)
    for k in range(len(sheets)):
        bykey[frozenset(int(blk[h]) for h in cells[k])].append(k)
    groups = []
    for key, ks in bykey.items():
        classes = []                       # parallel classes: members pairwise disjoint
        for k in ks:
            for cl in classes:
                if all(not (cells[k] & cells[m]) for m in cl):
                    cl.append(k)
                    break
            else:
                classes.append([k])
        groups += classes
    return [[e for k in g for e in sheets[k]] for g in groups]


_mesh_sheets0 = cb.mesh_sheets


if __name__ == "__main__":
    no_coarsen = "--no-coarsen" in sys.argv
    if no_coarsen:
        sys.argv.remove("--no-coarsen")
    use_macro = "--macro" in sys.argv          # macro moves on the fine mesh (implies no coarsening)
    if use_macro:
        sys.argv.remove("--macro"); no_coarsen = True
    slow_stats = "--slow-stats" in sys.argv    # original pure-Python _structure_stats (for comparison)
    if slow_stats:
        sys.argv.remove("--slow-stats")
    else:
        import fast_stats
        _stats0 = fast_stats.structure_stats
        cb._structure_stats = _stats0
    if use_macro:
        cb.mesh_sheets = macro_sheets
    # wrap the patch of beam_collapse: coarsen, then its beam search
    _patch = bm.patch
    feat_tol = None
    if "--feature-tol" in sys.argv:
        k = sys.argv.index("--feature-tol"); feat_tol = float(sys.argv[k + 1]); del sys.argv[k:k + 2]
    set_feature_curves(sys.argv[sys.argv.index("--input-vtk") + 1])

    def patch(*args, **kw):
        _patch(*args, **kw)
        beam = cb.collapse_mesh_sheets
        def collapse(P, hexes, labeller, **kw2):
            _FEAT["labeller"] = labeller
            # the feature guard: tolerance = max(given, value of the start mesh), like the Hausdorff ceiling
            f0 = feature_dist(P, hexes)
            if not no_coarsen:
                P, hexes = coarsen(P, hexes, labeller)
            f1 = feature_dist(P, hexes)
            _FEAT["tol"] = max(feat_tol if feat_tol is not None else 0.05, f0, f1)
            print(f"[coarse] feature curves -> block edges: fine start {f0:.4f}, coarse start {f1:.4f}, "
                  f"guard <= {_FEAT['tol']:.4f}", flush=True)
            cb._structure_stats = _stats_with_feat
            bm._ok = _ok_feat
            try:
                return beam(P, hexes, labeller, **kw2)
            finally:
                cb._structure_stats = _stats0
                bm._ok = _ok0
        cb.collapse_mesh_sheets = collapse
    bm.patch = patch
    src = open(HEX / "beam_collapse.py").read()
    main = src[src.index('if __name__ == "__main__":'):]
    exec(compile(main, str(HEX / "beam_collapse.py"), "exec"),
         {**vars(bm), "__name__": "__main__", "patch": patch})
