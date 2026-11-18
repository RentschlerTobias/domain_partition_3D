"""Global sheet collapse: beam search instead of greedy single steps.

`clean_blocks.collapse_mesh_sheets` takes the best SINGLE sheet per round
and stops as soon as no single step (excess faces, blocks) strictly
improves. Measured on one geometry (machine_0004) whose three runs share
the same singularity graph and the same raw base complex (84 blocks): 12,
22, and 75 blocks -- depending on which sheet falls first and when the
guard kicks in.

Here:
* Beam search over collapse orders (breadth `width`, depth `depth`).
  A step does NOT have to improve immediately; the best state reached
  anywhere in the search tree wins.
* Two guards:
  - "raw": as before, relative to the raw mesh from AlgoHex (no additional
    inverted cell, min sJ not below the raw value).
  - "abs": absolute floor for min sJ (`sj_floor`), independent of how good
    the raw mesh happened to be. The "raw" guard punishes a good raw mesh
    with a stricter guard -- base_b (min sJ 0.18) stayed at 75 blocks,
    base_a (min sJ -0.76) collapsed down to 12.
  Both keep the Hausdorff ceiling against the input surface.
  - "struct": structure only. Inverted cells of the intermediate hex mesh
    do not count at all -- the TFI refills every block freshly at the end.
    Quality serves only as the tie-break. This is the mode the dataset
    relabelling uses; the block faces are projected onto the exact
    geometry later, so the ceiling is a flat absolute Hausdorff tolerance.
* Candidates of a state are evaluated in parallel (fork pool).

Hook-in: `patch(width, depth, guard, sj_floor, workers)` replaces
`clean_blocks.collapse_mesh_sheets`; the rest of clean_blocks runs
unchanged (cavities, labels, merge/split, validator, export).

Result across 17 runs (12 perturbation runs of machine_0004, 1 smoke
sample, 2 base repeats, 2 more started from saved blocks.vtk incl. the
real cluster sample machine_0004_n2000): every run ends at 12 blocks
(12 cuboids, 0 inverted, VALID) with the SAME combinatorial topology,
where the greedy collapse produced 12, 22, or 75.
"""

import time
from multiprocessing import get_context

import numpy as np

from hexdp.core import clean_blocks as cb

_G = {}


def _eval(i):
    P, hexes, labeller, sheets = _G["P"], _G["hexes"], _G["labeller"], _G["sheets"]
    Q, H, _drop = cb.collapse_sheet(P, hexes, sheets[i])
    if not len(H):
        return i, None
    st = cb._structure_stats(Q, H, labeller)
    return i, st


def _children(P, hexes, labeller, workers):
    sheets = cb.mesh_sheets(hexes)
    _G.update(P=P, hexes=hexes, labeller=labeller, sheets=sheets)
    if workers > 1:
        with get_context("fork").Pool(workers) as pool:
            res = pool.map(_eval, range(len(sheets)))
    else:
        res = [_eval(i) for i in range(len(sheets))]
    return sheets, res


def _ok(st, bar, ceiling, guard, sj_floor):
    if st is None:
        return False
    if isinstance(ceiling, tuple):              # local tolerance: distance / allowed distance per input face
        if st["hausdorff_rel"] > ceiling[1] + 1e-9:
            return False
    elif st["hausdorff"] > ceiling + 1e-9:
        return False
    if st.get("non_ball", 0) > bar.get("non_ball", 0):
        # a block that is not a topological ball (wrapped around the blade hole or
        # glued to itself) cannot become a cuboid; with the keep-BL tolerance 0.1
        # the beam otherwise walks into 1- and 4-block states (0036, 0020)
        return False
    if guard == "struct":
        # structure only: inverted cells of the intermediate hex mesh do not
        # count, the TFI refills every block at the end (Gao repairs only then)
        return True
    if st["inverted"] > bar["inverted"]:
        return False
    floor = bar["min_sj"] if guard == "raw" else sj_floor
    return st["min_sj"] >= floor - 1e-9


def topology_issues(P, hexes):
    """Topological defects of the block structure a candidate mesh decomposes into, as kind -> count:
    a face shared by more than two blocks, a sheet with different cell counts, a hanging vertex, a block
    that is no structured grid. Inverted cells are not counted (the struct guard leaves them to the
    untangling). Measured: without this, 8 of 670 keep-BL samples ended with a block face shared by four
    blocks (all 8- and 14-block results), because quality and Hausdorff alone do not see it."""
    import contextlib
    import io
    from collections import Counter
    from hexdp.io.blocks_vtk import from_grids, structured_grids
    with contextlib.redirect_stdout(io.StringIO()):
        f2h, e2h = cb.bc.build_topology(hexes)
        sing = cb.bc.singular_edges(hexes, P, f2h, e2h)
        sheet_of, _n = cb.bfx.label_sheets(hexes, f2h, sing)
        blocks = cb.bc.blocks_from_cut(hexes, f2h, set(sheet_of))
    bid = np.zeros(len(hexes), int)
    for i, b in enumerate(blocks):
        bid[list(b)] = i
    try:
        grids, _ = structured_grids(P, hexes, bid)
    except ValueError:
        return Counter(not_structured=1)
    return Counter(i.kind for i in from_grids(P, grids).check() if i.kind != "inverted_cells")


def _key(st):
    # fewer excess faces, then fewer blocks, then better quality
    return (st["excess"], st["blocks"], -st["min_sj"])


class RawMeshOffSurface(RuntimeError):
    """The raw hex mesh does not match the input surface; nothing to collapse."""


def beam_collapse(P, hexes, labeller, width=6, depth=6, guard="abs",
                  sj_floor=0.05, workers=6, guard_ref=None, verbose=True,
                  hausdorff_tol=0.05, max_base_hausdorff=0.15, topology_guard=True):
    base = cb._structure_stats(P, hexes, labeller)
    base_topo = topology_issues(P, hexes) if topology_guard else None
    # The ceiling below follows a raw Hausdorff above the tolerance. On a raw
    # mesh that is off the surface altogether that frees every collapse: 17
    # of 905 papaya samples (raw Hausdorff 0.26-0.50, the 862 good ones
    # 0.025-0.064) collapsed to ONE 7-faced block, and clean_blocks then died
    # on the median of zero sheets. Refuse such a mesh up front.
    if max_base_hausdorff is not None and base["hausdorff"] > max_base_hausdorff:
        raise RawMeshOffSurface(
            f"raw hex mesh is off the input surface: Hausdorff {base['hausdorff']:.4f} "
            f"> {max_base_hausdorff} ({base['blocks']} blocks, {base['cells']} cells)")
    bar = dict(guard_ref) if guard_ref else {"inverted": base["inverted"],
                                             "min_sj": base["min_sj"]}
    bar.setdefault("non_ball", base.get("non_ball", 0))
    # Hausdorff ceiling: the raw value is itself an accident of the tet
    # meshing (base_a 0.0395, base_b 0.0400) and decided alone over 22
    # versus 75 blocks -- at base_b the 84->22 sheet came in 0.0019 above
    # it. With "abs" an absolute tolerance applies; the block faces are
    # projected onto the exact geometry later anyway.
    ceiling = base["hausdorff"] if guard == "raw" else max(base["hausdorff"], hausdorff_tol)
    if verbose and getattr(labeller, "allow", None) is not None:
        al = labeller.allow
        print(f"[beam] local tolerance per input face: min {al.min():.4f}, median {np.median(al):.4f}, "
              f"max {al.max():.4f}; {100 * np.mean(al < al.max() - 1e-12):.1f}% of the faces below the maximum")
    if guard != "raw" and getattr(labeller, "allow", None) is not None:
        # local tolerance (allowed distance per input face from the local gap): relative ceiling 1, or the raw
        # mesh's own relative deviation if that is larger
        ceiling = ("rel", max(base["hausdorff_rel"], 1.0))
    # with "struct" the collapse may sacrifice quality -> quality only breaks ties
    if verbose:
        fl = bar["min_sj"] if guard == "raw" else sj_floor
        print(f"[beam] start {base}")
        print(f"[beam] guard={guard}: inverted <= {bar['inverted']}, min sJ >= "
              f"{fl:.4f}, Hausdorff <= "
              f"{(f'{ceiling[1]:.3f} x local tolerance' if isinstance(ceiling, tuple) else f'{ceiling:.5f}')}; "
              f"width {width}, depth {depth}")
    beam = [(base, P, hexes, [])]
    best = (base, P, hexes, [])
    seen = set()
    for d in range(depth):
        t0 = time.time()
        cand = []
        for st0, Pb, Hb, path in beam:
            sheets, res = _children(Pb, Hb, labeller, workers)
            for i, st in res:
                if not _ok(st, bar, ceiling, guard, sj_floor):
                    continue
                sig = (st["cells"], st["blocks"], st["excess"], round(st["min_sj"], 6))
                if sig in seen:
                    continue
                seen.add(sig)
                cand.append((st, Pb, Hb, sheets[i], path + [i]))
        if not cand:
            if verbose:
                print(f"[beam] depth {d + 1}: no admissible collapse, stop")
            break
        cand.sort(key=lambda c: _key(c[0]))
        beam, rejected = [], 0
        for st, Pb, Hb, es, path in cand:
            if len(beam) >= width:
                break
            Q, H, _ = cb.collapse_sheet(Pb, Hb, es)
            if topology_guard:
                # no topological defect the raw structure did not have (checked lazily: only for
                # candidates that would enter the beam)
                topo = topology_issues(Q, H)
                if any(n > base_topo.get(k, 0) for k, n in topo.items()):
                    rejected += 1
                    continue
            beam.append((st, Q, H, path))
            if _key(st) < _key(best[0]):
                best = (st, Q, H, path)
        if not beam:
            if verbose:
                print(f"[beam] depth {d + 1}: {rejected} candidates rejected by the topology guard, stop")
            break
        if verbose:
            tops = ", ".join(f"{b[0]['blocks']}bl/ex{b[0]['excess']}/sJ{b[0]['min_sj']:.2f}"
                             for b in beam)
            rej = f", {rejected} rejected by the topology guard" if rejected else ""
            print(f"[beam] depth {d + 1}: {len(cand)} admissible{rej}, beam [{tops}]  "
                  f"best {best[0]['blocks']} blocks ({time.time() - t0:.0f}s)")
    st, P, hexes, path = best
    if verbose:
        print(f"[beam] result: {st} via sheets {path}")
    log = [dict(base, sheet=None), dict(st, sheet=path)]
    return P, hexes, log


def min_patch_gap(points, tris, labels):
    """Smallest distance between two labelled surfaces that share no vertex.

    The absolute Hausdorff tolerance must stay below the narrowest gap of the
    input, or a collapse can close it: canada guide vane, envelope 0.022 from
    the inlet, collapsed 12 -> 5 blocks at Hausdorff 0.045 (< tol 0.05) by
    letting the O-ring open into the inlet. Vertex-to-vertex distance, so it
    overestimates the true gap by at most half an edge length."""
    from scipy.spatial import cKDTree
    labels = np.asarray(labels)
    verts = {int(l): np.unique(tris[labels == l]) for l in np.unique(labels)}
    best, pair = np.inf, None
    ls = sorted(verts)
    for i, a in enumerate(ls):
        tree = cKDTree(points[verts[a]])
        for b in ls[i + 1:]:
            if np.intersect1d(verts[a], verts[b]).size:
                continue
            d, _ = tree.query(points[verts[b]])
            if d.min() < best:
                best, pair = float(d.min()), (a, b)
    return best, pair


def local_allow(points, tris, labels, tol, frac):
    """Allowed boundary deviation per input triangle: min(tol, frac x local gap), the local gap being the distance
    from the triangle's vertices to the nearest vertex of a surface that shares no vertex with its own surface
    (adjacent surfaces meet at gap 0 by construction). Away from narrow gaps the full tolerance applies."""
    from scipy.spatial import cKDTree
    labels = np.asarray(labels)
    verts = {int(l): np.unique(tris[labels == l]) for l in np.unique(labels)}
    trees = {l: cKDTree(points[v]) for l, v in verts.items()}
    gap_v = np.full(len(points), np.inf)
    for a, va in verts.items():
        for b, vb in verts.items():
            if a == b or np.intersect1d(va, vb).size:
                continue
            d, _ = trees[b].query(points[va])
            gap_v[va] = np.minimum(gap_v[va], d)
    gap_t = gap_v[tris].min(1)
    return np.minimum(tol, frac * gap_t)


def patch(width=6, depth=6, guard="abs", sj_floor=0.05, workers=6, hausdorff_tol=0.05,
          max_base_hausdorff=0.15, topology_guard=True):
    def collapse(P, hexes, labeller, max_rounds=4, verbose=True, repair=None,
                 guard_ref=None):
        P, hexes, log = beam_collapse(P, hexes, labeller, width, depth, guard,
                                      sj_floor, workers, guard_ref, verbose,
                                      hausdorff_tol, max_base_hausdorff, topology_guard)
        if repair is not None:
            P = repair(P.copy(), hexes)
        return P, hexes, log
    cb.collapse_mesh_sheets = collapse


if __name__ == "__main__":
    import argparse
    import sys
    from hexdp.strategies.algohex import tet_prep as v5
    ap = argparse.ArgumentParser()
    ap.add_argument("ovm")
    ap.add_argument("--input-vtk", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--width", type=int, default=6)
    ap.add_argument("--depth", type=int, default=6)
    ap.add_argument("--guard", choices=("raw", "abs", "struct"), default="abs")
    ap.add_argument("--sj-floor", type=float, default=0.05)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--hausdorff-tol", type=float, default=0.05)
    ap.add_argument("--untangle-rounds", type=int, default=6)
    ap.add_argument("--max-raw-hausdorff", type=float, default=0.15,
                    help="refuse a raw hex mesh farther than this from the input surface")
    ap.add_argument("--gap-mode", choices=("local", "global"), default="local",
                    help="local: the tolerance per input face is min(tol, gap-frac x local gap); global: one "
                         "tolerance min(tol, gap-frac x narrowest gap) for the whole surface (old behaviour)")
    ap.add_argument("--gap-frac", type=float, default=0.5,
                    help="cap the Hausdorff tolerance at this fraction of the narrowest gap "
                         "between two labelled input surfaces that share no vertex; 0 = off")
    ap.add_argument("--max-drift-frac", type=float, default=0.5,
                    help="untangling may move boundary vertices up to this fraction of the Hausdorff tolerance off "
                         "the input surface to repair inverted cells inherited from the hex mesh (0: the old 1e-3)")
    ap.add_argument("--force-ball", action="store_true",
                    help="cut non-ball base blocks (a block wrapped around a through hole) by extra face sheets")
    ap.add_argument("--boundary-singularities", action="store_true",
                    help="boundary edges with a hex valence other than 2 (and sharp edges) also break blocks: for "
                         "meshes with a regular boundary such as a polycube map (default: interior singular edges)")
    ap.add_argument("--no-topology-guard", action="store_true",
                    help="do not reject collapses that create topological defects (old behaviour)")
    ap.add_argument("--allow-invalid", action="store_true",
                    help="write and exit 0 even if the final block structure fails BlockStructure.check")
    a = ap.parse_args()
    tol = a.hausdorff_tol
    if a.gap_frac > 0:
        import meshio
        _m = meshio.read(a.input_vtk)
        _cd = next(iter(_m.cell_data))
        _T = np.vstack([c.data for c in _m.cells if c.type == "triangle"])
        _L = np.concatenate([np.asarray(d).ravel() for c, d in zip(_m.cells, _m.cell_data[_cd])
                             if c.type == "triangle"])
        gap, pair = min_patch_gap(np.asarray(_m.points, float), _T, _L)
        tol = min(a.hausdorff_tol, a.gap_frac * gap)
        print(f"[beam] narrowest gap {gap:.4f} between labels {pair} -> Hausdorff tol "
              f"{tol:.4f} (flag {a.hausdorff_tol})")
        if a.gap_mode == "local":
            # computed from the labeller's own surface (same file, but no assumption on the triangle order)
            cb.SURFACE_ALLOW = lambda P, T, ids: local_allow(P, T, ids, a.hausdorff_tol, a.gap_frac)
    if a.force_ball:
        cb.FORCE_BALL = True
    if a.boundary_singularities:
        from hexdp.core import base_complex as _bcx
        _bcx.INTERIOR_ONLY = False
    if a.max_drift_frac > 0:
        cb.UNTANGLE_MAX_DRIFT = max(cb.UNTANGLE_MAX_DRIFT, a.max_drift_frac * tol)
        print(f"[beam] untangle max boundary drift {cb.UNTANGLE_MAX_DRIFT:.4f}")
    patch(a.width, a.depth, a.guard, a.sj_floor, a.workers, tol,
          a.max_raw_hausdorff, not a.no_topology_guard)
    if a.ovm.endswith(".vtk"):
        # start from a saved blocks.vtk (fine hex mesh AFTER the greedy
        # collapse): BlockStructure otherwise reads .ovm only. The block IDs
        # are discarded -- the blocks follow from the sheets.
        import meshio
        from hexdp.io import ovm as ovm_io
        _m = meshio.read(a.ovm)
        _H = np.vstack([c.data for c in _m.cells if c.type == "hexahedron"]).astype(np.int64)
        _P = np.asarray(_m.points, float)
        ovm_io.read_ovm = lambda path: (_P, None, None, None)
        ovm_io.ovm_to_cells = lambda P, e, f, poly: (_H, 0)
    try:
        S, before, after, v = cb.postprocess(
            a.ovm, a.input_vtk, v5.NAMES, collapse_rounds=1, untangle_mesh=True,
            untangle_rounds=a.untangle_rounds)
    except RawMeshOffSurface as e:
        print(f"[beam] REJECT: {e}")
        sys.exit(3)
    cb.write_blocks(S, a.out)
    print("\nbefore/after:")
    for k in before:
        print(f"  {k:22s} {before[k]}  ->  {after[k]}")
    # the written structure must pass the same check an agent or the dataset applies to it
    from hexdp.core.blockstructure import invalid
    from hexdp.io.blocks_vtk import read_blocks_vtk
    bad = invalid(read_blocks_vtk(a.out).check())
    if bad:
        for i in bad[:10]:
            print(f"[beam] INVALID {i.kind}: {i.message}")
        if not a.allow_invalid:
            print(f"[beam] REJECT: the block structure is invalid ({len(bad)} issues); see above")
            sys.exit(4)
    sys.exit(0)
