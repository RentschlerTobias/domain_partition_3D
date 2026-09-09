"""Optimal block merging as a weighted exact cover, solved as an ILP.

`clean_blocks.merge_non_cuboids` merges greedily: it takes the neighbour that
helps most right now, and a merge it makes can block a better one later. The
sheet collapse has the same shape of problem and stops the same way -- "no
sheet improves" -- because both accept only strict local improvement.

Merging is the half of this that can be solved exactly, and cheaply, because
**a block merge changes no geometry at all**. It is a union-find union on the
cut set: the interface between two blocks stops separating them, no vertex
moves, no cell changes. Inverted cells, scaled Jacobians and the Hausdorff
distance are therefore invariant under everything this module does, which
removes the whole class of "the metric moved the right way while the mesh got
worse" failures. The only thing at stake is the block topology.

That makes the problem a set-partitioning ILP:

    variables   one binary per candidate group of blocks
    constraint  every block is covered by exactly one selected group
    objective   minimise the number of resulting blocks, plus a penalty for
                every block that stays defective (non-cuboid, or too small)

Candidate groups are restricted to connected sets whose union is still a valid
block, and singletons are always candidates, so the cover stays feasible and
the solver can always answer "leave this one alone". The optimum is exact over
the enumerated candidates -- not over all conceivable block structures, since
minimum cuboid decomposition is NP-hard in general and the enumeration is
bounded by `--max-group`.

**What counts as "still a valid block" decides everything here.** Two criteria
are available and they disagree completely on this data:

* `--valid cuboid` -- `clean_blocks.cuboid_status`: exactly 6 patches, each
  adjacent to 4 others. Patches are keyed by what lies on the other side, so
  two faces of the merged region pointing at two different neighbour blocks
  stay two patches. On v16 not one of 1650 candidate groups passes: the
  unions come out with 8 to 11 patches. Merging is impossible under it.
* `--valid lattice` (default) -- `tfi.block_lattice`: the merged cells form a
  structured (a x b x c) grid. On v16 **all 44** seed-neighbour pairs pass.

The second is the one the deliverable actually needs. TFI fills a block from
its lattice; it does not care how many distinct neighbours sit on one side of
it. Judging merges by patch count answers a question nobody asked and reports
"no merge possible" on a structure where every candidate is a perfect grid.

    PY merge_ilp.py BLOCKS.vtk [--valid lattice|cuboid] [--max-group 4]
                    [--min-cells 10] [--apply OUT.vtk]
"""

import argparse
import sys
from collections import defaultdict, deque
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "experimentell" / "hex3d_algohex"))

import clean_blocks as cb                                              # noqa: E402
import tfi                                                             # noqa: E402


def block_graph(S, cells_of, root_of):
    """{block: {neighbouring blocks}} from the shared faces."""
    adj = defaultdict(set)
    for r in cells_of:
        adj[r] = set()
    for fk, hs in S.f2h.items():
        if len(hs) != 2:
            continue
        a, b = root_of[hs[0]], root_of[hs[1]]
        if a != b:
            adj[a].add(b)
            adj[b].add(a)
    return adj


def group_patches(S, group, cells_of, root_of):
    """The boundary patches of a merged set of blocks.

    Mirrors `BlockStructure.face_records` + `_patch_split`, but asks whether
    the neighbour is in the GROUP rather than whether it is one particular
    block -- which is exactly what a merge does to the face records."""
    rec = []
    for r in group:
        for h in cells_of[r]:
            c = S.hexes[h]
            for fc in cb.HF:
                lp = [int(c[i]) for i in fc]
                k = frozenset(lp)
                nb = [x for x in S.f2h[k] if x != h]
                if not nb:
                    rec.append((("P", S.surf_of[k]), lp))
                elif root_of[nb[0]] not in group:
                    rec.append((("B", root_of[nb[0]]), lp))
    return cb._patch_split(rec)


def enumerate_groups(adj, seeds, max_group, cap=200000):
    """Connected block sets up to `max_group` that contain at least one seed.

    Which blocks may seed a group decides the result far more than how large
    groups may get. Measured on v11:

        max-group  seeds       blocks
            4      defective   16 -> 12
            4      all         16 ->  6
            6      defective   16 -> 12
            6      all         16 ->  6

    Group size is irrelevant here; the seed set is everything. Seeding only
    from defective blocks expresses a REPAIR objective -- fix the tiny and
    non-cuboid ones, leave the rest alone -- and it costs six blocks when the
    goal is actually the smallest structure.

    The same numbers also show why this is solved globally rather than
    iteratively: 16 -> 6 in one solve beats 16 -> 12 -> 7 in two."""
    out, seen = [], set()
    for s in seeds:
        frontier = [frozenset((s,))]
        while frontier:
            nxt = []
            for g in frontier:
                if len(g) >= max_group:
                    continue
                for r in g:
                    for nb in adj[r]:
                        if nb in g:
                            continue
                        ng = g | {nb}
                        if ng in seen:
                            continue
                        seen.add(ng)
                        out.append(ng)
                        nxt.append(ng)
                        if len(out) >= cap:
                            print(f"[merge_ilp] candidate cap {cap} reached")
                            return out
            frontier = nxt
    return out


def feature_edges(S):
    """Mesh edges on the domain boundary where two surface labels meet.

    These are the feature curves of the input, seen from the hex mesh: the
    blade-root and blade-tip rings between the cut surfaces, the inlet and
    outlet rims, the periodic seams. A block edge that lies on one of them is
    load-bearing, because the geometry has a real kink there -- the hub ring's
    median dihedral angle is 81 degrees."""
    e2s = defaultdict(set)
    for (fk, _h), lp in zip(S.bnd, S.bnd_loops):
        s = S.surf_of[fk]
        for a in range(4):
            u, v = int(lp[a]), int(lp[(a + 1) % 4])
            e2s[(min(u, v), max(u, v))].add(s)
    return {e for e, s in e2s.items() if len(s) > 1}


def crosses_feature(S, g, cells_of, root_of, feats):
    """Would merging this group dissolve a block edge that sits on a feature?

    Merging removes the interfaces between the blocks of the group. If one of
    those interfaces meets the domain boundary along a feature curve, the
    merged block spans a geometric kink, and the Gordon-Hall fill then bridges
    that kink linearly through the block's interior.

    OFF BY DEFAULT, because measuring it refuted the idea it was built for.
    The hypothesis was that v11's [4, 9, 13, 15] merge folds on refill because
    it spans the blade-root ring. It does not cross a feature edge at all --
    `crosses_feature` returns False for it -- so this check never touched it.
    What it does reject is [4, 5, 7, 13] and [9, 11], which are precisely the
    two merges that make up the 12-block structure that refills cleanly. And
    with the check on, the 6-block refill went from 2 inverted cells to 40.

    So it blocks the good merges and lets the bad one through, on this data.
    Kept because "block edges belong on feature curves" is standard practice
    in block-structured meshing and may earn its keep on another geometry --
    but nothing here justifies switching it on."""
    inner = set()
    for r in g:
        for h in cells_of[r]:
            c = S.hexes[h]
            for fc in cb.HF:
                lp = [int(c[i]) for i in fc]
                k = frozenset(lp)
                nb = [x for x in S.f2h[k] if x != h]
                if nb and root_of[nb[0]] in g and root_of[nb[0]] != r:
                    for a in range(4):
                        u, v = lp[a], lp[(a + 1) % 4]
                        inner.add((min(u, v), max(u, v)))
    return bool(inner & feats)


def is_valid_merge(S, g, cells_of, root_of, mode, feats=None):
    """Would this group still be a usable block after merging?

    `lattice` also rejects a merge that would fold the block around a wall --
    two ADJACENT faces landing on the same physical surface -- reusing
    `duplicate_surface_faces`, the guard `merge_non_cuboids` already applies.
    Without it, minimising the block count would happily weld a block across
    the passage and hand TFI a block with the inlet on two touching sides."""
    if feats and crosses_feature(S, g, cells_of, root_of, feats):
        return False
    if mode == "cuboid":
        st, _d = cb.cuboid_status(group_patches(S, g, cells_of, root_of))
        return st == "cuboid"
    cells = [h for r in g for h in cells_of[r]]
    if tfi.block_lattice(S.hexes, cells, S.f2h) is None:
        return False
    _opp, bad = cb.duplicate_surface_faces(
        group_patches(S, g, cells_of, root_of))
    return not bad


def solve(S, max_group=4, min_cells=10, mode="lattice", seed="all",
          keep_features=False, banned=None, verbose=True):
    cells_of = S.cells_of()
    root_of = {h: S.root(S.blk0[h]) for h in range(len(S.hexes))}
    blocks = sorted(cells_of)
    size = {r: len(cells_of[r]) for r in blocks}
    adj = block_graph(S, cells_of, root_of)

    status = {r: cb.cuboid_status(S.patches(r))[0] for r in blocks}
    defective = [r for r in blocks
                 if status[r] != "cuboid" or size[r] < min_cells]
    seeds = blocks if seed == "all" else defective
    if verbose:
        print(f"[merge_ilp] {len(blocks)} blocks, {sum(len(v) for v in adj.values()) // 2} "
              f"adjacencies, {len(defective)} defective (non-cuboid or "
              f"< {min_cells} cells), seeding from {len(seeds)}")
    if not seeds:
        print("[merge_ilp] nothing to merge")
        return None

    cand = enumerate_groups(adj, seeds, max_group)
    if verbose:
        print(f"[merge_ilp] {len(cand)} connected candidate groups of size "
              f"2..{max_group}, testing them against '{mode}' ...")

    feats = feature_edges(S) if keep_features else None
    if verbose and feats:
        print(f"[merge_ilp] {len(feats)} boundary edges lie on a feature "
              f"curve; merges that would dissolve one are rejected")
    banned = banned or set()
    valid = [g for g in cand if g not in banned
             and is_valid_merge(S, g, cells_of, root_of, mode, feats)]
    if verbose:
        print(f"[merge_ilp] {len(valid)} of them stay a valid block")

    # singletons keep the status quo and make the cover always feasible
    groups = [frozenset((r,)) for r in blocks] + valid
    W_DEFECT, W_TINY = 10.0, 5.0
    cost = []
    for g in groups:
        n = sum(size[r] for r in g)
        if len(g) == 1:
            r = next(iter(g))
            c = 1.0 + (W_DEFECT if status[r] != "cuboid" else 0.0) \
                + (W_TINY if n < min_cells else 0.0)
        else:
            c = 1.0 + (W_TINY if n < min_cells else 0.0)
        cost.append(c)

    from scipy.optimize import milp, LinearConstraint, Bounds
    idx = {r: i for i, r in enumerate(blocks)}
    rows, cols = [], []
    for j, g in enumerate(groups):
        for r in g:
            rows.append(idx[r])
            cols.append(j)
    from scipy.sparse import coo_matrix
    A = coo_matrix((np.ones(len(rows)), (rows, cols)),
                   shape=(len(blocks), len(groups))).tocsc()
    res = milp(c=np.array(cost),
               constraints=[LinearConstraint(A, 1, 1)],
               integrality=np.ones(len(groups)),
               bounds=Bounds(0, 1))
    if not res.success:
        print(f"[merge_ilp] ILP failed: {res.message}")
        return None

    chosen = [groups[j] for j in np.where(np.round(res.x) > 0.5)[0]]
    merges = [g for g in chosen if len(g) > 1]
    left_defect = [next(iter(g)) for g in chosen
                   if len(g) == 1 and status[next(iter(g))] != "cuboid"]
    left_tiny = [g for g in chosen
                 if sum(size[r] for r in g) < min_cells]
    if verbose:
        print(f"[merge_ilp] optimum: {len(chosen)} blocks "
              f"({len(blocks)} before), {len(merges)} merges")
        for g in sorted(merges, key=lambda g: -len(g)):
            print(f"[merge_ilp]   merge {sorted(g)} -> "
                  f"{sum(size[r] for r in g)} cells")
        print(f"[merge_ilp] still non-cuboid: {len(left_defect)}; "
              f"still under {min_cells} cells: {len(left_tiny)}")
    return {"chosen": chosen, "blocks_before": len(blocks),
            "blocks_after": len(chosen), "merges": merges,
            "left_defect": left_defect, "left_tiny": len(left_tiny),
            "cells_of": cells_of}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("blocks")
    ap.add_argument("--input-vtk",
                    default=str(REPO / "data" / "T1_9" / "T1_9_tet_v5.vtk"))
    ap.add_argument("--valid", choices=("lattice", "cuboid"),
                    default="lattice",
                    help="what a merged group must still be; see the module "
                         "docstring -- 'cuboid' admits nothing on this data")
    ap.add_argument("--keep-features", action="store_true",
                    help="reject merges that dissolve a block edge sitting on "
                         "a feature curve. Off by default -- measured to block "
                         "the merges that work and not the one that folds; "
                         "see crosses_feature")
    ap.add_argument("--probe-interfaces", action="store_true",
                    help="test each block interface on its own: merge just "
                         "those two blocks, refill, and report whether the "
                         "boundary error grows. Linear in interfaces, and it "
                         "measures what block metrics cannot see")
    ap.add_argument("--probe-tol", type=float, default=1.5,
                    help="how much the boundary error may grow and still count "
                         "as safe (default 1.5x)")
    ap.add_argument("--probe-h", type=float, default=0.05,
                    help="target cell size for the probe refills")
    ap.add_argument("--seed", choices=("all", "defective"), default="all",
                    help="which blocks may seed a candidate group. 'all' "
                         "minimises the block count (v11: 16 -> 6), "
                         "'defective' only repairs tiny and non-cuboid "
                         "blocks and leaves the rest alone (16 -> 12)")
    ap.add_argument("--max-group", type=int, default=4)
    ap.add_argument("--min-cells", type=int, default=10)
    ap.add_argument("--apply", default=None, metavar="OUT_VTK")
    a = ap.parse_args()

    S, bid, lab = cb.read_blocks_vtk(a.blocks, a.input_vtk)
    if a.probe_interfaces:
        probe_interfaces(S, lab, a.probe_h, a.probe_tol)
        return
    out = solve(S, a.max_group, a.min_cells, a.valid, a.seed,
                keep_features=a.keep_features)
    if out is None or not a.apply:
        return
    # apply the merges to the structure itself and export through
    # clean_blocks.write_blocks, so the result carries the same artifacts as
    # every other deliverable: .vtk, .msh with the block edges as 1D cells,
    # _edges.vtk, _nfaces.vtk and _quality.vtk
    for g in out["merges"]:
        it = iter(sorted(g))
        first = next(it)
        for r in it:
            S.merge(first, r)
    cb.write_blocks(S, a.apply,
                    f"{out['blocks_after']} blocks after the merge ILP")



# --------------------------------------------------------------------------
# which block edges are load-bearing: probe one interface at a time
# --------------------------------------------------------------------------

def _boundary_faces_of(S, vert):
    """Which of a lattice block's six faces lie on the domain boundary."""
    out = []
    for ax in (0, 1, 2):
        for side in (0, 1):
            G = tfi.face_grid(vert, ax, side)
            n = 0
            for i in range(G.shape[0] - 1):
                for j in range(G.shape[1] - 1):
                    fk = frozenset((int(G[i, j]), int(G[i + 1, j]),
                                    int(G[i + 1, j + 1]), int(G[i, j + 1])))
                    if len(S.f2h.get(fk, ())) == 1:
                        n += 1
            if n:
                out.append((ax, side))
    return out


def _refill_boundary_error(S, vert, lab, target_h, samples=4):
    """Refill one block at `target_h` and measure how far its refilled
    boundary lands from the input surface.

    This is the quantity a merge has to preserve. Block quality metrics do
    not see it: on cand_002 the merged structure refills to a min scaled
    Jacobian of 0.5668 -- better than the unmerged 0.5470 -- while its
    boundary sits 0.47 off the geometry instead of 0.023. The cells are
    beautifully shaped and in the wrong place."""
    dims = np.array(vert.shape) - 1
    counts = [max(1, int(round(tfi.axis_length(S.P, vert, ax) / target_h)))
              for ax in (0, 1, 2)]
    X = tfi.refill_block(S.P, vert, counts)
    pts = []
    for ax, side in _boundary_faces_of(S, vert):
        sl = [slice(None)] * 3
        sl[ax] = 0 if side == 0 else -1
        F = X[tuple(sl)]                       # (n, m, 3) refilled face
        # sample inside each new quad, where chordal error lives
        for a in np.linspace(0.25, 0.75, samples):
            for b in np.linspace(0.25, 0.75, samples):
                pts.append((1 - a) * (1 - b) * F[:-1, :-1]
                           + a * (1 - b) * F[1:, :-1]
                           + a * b * F[1:, 1:]
                           + (1 - a) * b * F[:-1, 1:])
    if not pts:
        return 0.0
    Q = np.vstack([p.reshape(-1, 3) for p in pts])
    d = np.linalg.norm(cb.project_to_surface(lab, Q) - Q, axis=1)
    return float(np.percentile(d, 95))


def probe_interfaces(S, lab, target_h=0.05, tol=1.5, verbose=True):
    """Test every block interface on its own: is it load-bearing?

    Merging two blocks removes the block edge between them, and the refill
    then has to bridge whatever the geometry does there. Rather than guess
    which edges matter -- an earlier attempt guessed "the ones on feature
    curves" and was wrong -- each interface is merged in isolation, refilled,
    and kept only if the boundary error does not grow by more than `tol`.

    Linear in the number of interfaces, where enumerating block subsets is
    combinatorial: 34 probes here against a subset enumeration that did not
    finish in ten minutes at group size 6."""
    cells_of = S.cells_of()
    root_of = {h: S.root(S.blk0[h]) for h in range(len(S.hexes))}
    adj = block_graph(S, cells_of, root_of)
    lat = {}
    for r in cells_of:
        got = tfi.block_lattice(S.hexes, cells_of[r], S.f2h)
        if got is not None:
            lat[r] = tfi._fix_handedness(S.P, got[1])
    err = {r: _refill_boundary_error(S, lat[r], lab, target_h) for r in lat}

    pairs = sorted({(min(a, b), max(a, b)) for a in adj for b in adj[a]})
    rows = []
    for a, b in pairs:
        if a not in lat or b not in lat:
            continue
        got = tfi.block_lattice(S.hexes, list(cells_of[a]) + list(cells_of[b]),
                                S.f2h)
        if got is None:
            rows.append({"pair": (a, b), "safe": False, "why": "no lattice"})
            continue
        v = tfi._fix_handedness(S.P, got[1])
        e = _refill_boundary_error(S, v, lab, target_h)
        base = max(err[a], err[b])
        safe = e <= max(base * tol, base + 1e-6)
        rows.append({"pair": (a, b), "safe": bool(safe), "sep": base,
                     "merged": e, "ratio": e / base if base else float("inf")})
    if verbose:
        good = [r for r in rows if r.get("safe")]
        print(f"[merge_ilp] probed {len(rows)} interfaces at h={target_h}: "
              f"{len(good)} keep the boundary within {tol:.1f}x, "
              f"{len(rows) - len(good)} are load-bearing")
        for r in sorted(rows, key=lambda r: -r.get("ratio", 0))[:12]:
            if "ratio" in r:
                print(f"[merge_ilp]   {str(r['pair']):>10}  separate "
                      f"{r['sep']:.5f} -> merged {r['merged']:.5f}  "
                      f"({r['ratio']:5.1f}x) {'SAFE' if r['safe'] else ''}")
            else:
                print(f"[merge_ilp]   {str(r['pair']):>10}  {r['why']}")
    return rows

if __name__ == "__main__":
    main()


def class_spread(S, chosen, cells_of):
    """Largest ratio of axis lengths inside one direction class.

    The quantity that decides whether a merged structure can still be
    refilled. A direction class carries ONE cell count for all its axes, so a
    class holding axes of 0.06 and 1.01 has no good answer: count 1 makes the
    long axis 20x too coarse, count 20 makes the short one 16x too fine.

    Measured on cand_002: the 22-block structure keeps every class within
    1.3x, and refills to a boundary error of 0.004. Merged to 7 blocks the
    worst class spreads 16.7x, eight axes come out more than twice as fine as
    they need, and the boundary error is 0.31 -- while every block metric
    still looks excellent (min scaled Jacobian 0.5668, 0 inverted). This is
    the number that sees the damage, and it costs milliseconds where probing
    by refill costs minutes."""
    import numpy as np
    P, H = S.P, S.hexes
    bid = np.zeros(len(H), int)
    for i, g in enumerate(chosen):
        for r in g:
            for h in cells_of[r]:
                bid[h] = i
    f2h, _e = __import__("base_complex").build_topology(H)
    lat, missing = tfi.lattices(P, H, bid, f2h, verbose=False)
    if missing:
        return float("inf"), None
    classes = tfi.direction_classes(lat, f2h, H, bid, verbose=False)
    worst, worst_ci = 0.0, None
    for ci, c in enumerate(classes):
        L = [tfi.axis_length(P, lat[r][1], ax) for r, ax in c]
        sp = max(L) / max(min(L), 1e-12)
        if sp > worst:
            worst, worst_ci = sp, ci
    return worst, (classes[worst_ci] if worst_ci is not None else None)


def solve_bounded_spread(S, max_spread=2.0, max_group=4, min_cells=10,
                         rounds=12, verbose=True):
    """Merge as far as the refill can still follow.

    Solve the cover, measure the resulting class spread, and if it exceeds
    `max_spread` forbid the largest merge feeding the offending class and
    solve again. Greedy only in which merge it drops; every round is an exact
    cover over the remaining candidates, and each costs a solve plus a spread
    evaluation -- no refill anywhere."""
    cells_of = S.cells_of()
    banned, best = set(), None
    for rnd in range(rounds):
        out = solve(S, max_group, min_cells, "lattice", "all", False,
                    banned=banned, verbose=False)
        if out is None:
            break
        sp, cls = class_spread(S, out["chosen"], cells_of)
        if verbose:
            print(f"[merge_ilp] round {rnd + 1}: {out['blocks_after']} blocks, "
                  f"worst class spread {sp:.1f}x, {len(banned)} merges banned")
        best = out
        if sp <= max_spread:
            return out, sp
        if cls is None:
            break
        guilty = {r for r, _ax in cls}
        cand = [g for g in out["merges"] if guilty & set(g)]
        if not cand:
            break
        banned.add(max(cand, key=len))
    return best, None
