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

import clean_blocks as cb

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
    if st["hausdorff"] > ceiling + 1e-9:
        return False
    if guard == "struct":
        # structure only: inverted cells of the intermediate hex mesh do not
        # count, the TFI refills every block at the end (Gao repairs only then)
        return True
    if st["inverted"] > bar["inverted"]:
        return False
    floor = bar["min_sj"] if guard == "raw" else sj_floor
    return st["min_sj"] >= floor - 1e-9


def _key(st):
    # fewer excess faces, then fewer blocks, then better quality
    return (st["excess"], st["blocks"], -st["min_sj"])


def beam_collapse(P, hexes, labeller, width=6, depth=6, guard="abs",
                  sj_floor=0.05, workers=6, guard_ref=None, verbose=True,
                  hausdorff_tol=0.05):
    base = cb._structure_stats(P, hexes, labeller)
    bar = dict(guard_ref) if guard_ref else {"inverted": base["inverted"],
                                             "min_sj": base["min_sj"]}
    # Hausdorff ceiling: the raw value is itself an accident of the tet
    # meshing (base_a 0.0395, base_b 0.0400) and decided alone over 22
    # versus 75 blocks -- at base_b the 84->22 sheet came in 0.0019 above
    # it. With "abs" an absolute tolerance applies; the block faces are
    # projected onto the exact geometry later anyway.
    ceiling = base["hausdorff"] if guard == "raw" else max(base["hausdorff"], hausdorff_tol)
    # with "struct" the collapse may sacrifice quality -> quality only breaks ties
    if verbose:
        fl = bar["min_sj"] if guard == "raw" else sj_floor
        print(f"[beam] start {base}")
        print(f"[beam] guard={guard}: inverted <= {bar['inverted']}, min sJ >= "
              f"{fl:.4f}, Hausdorff <= {ceiling:.5f}; width {width}, depth {depth}")
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
        beam = []
        for st, Pb, Hb, es, path in cand[:width]:
            Q, H, _ = cb.collapse_sheet(Pb, Hb, es)
            beam.append((st, Q, H, path))
            if _key(st) < _key(best[0]):
                best = (st, Q, H, path)
        if verbose:
            tops = ", ".join(f"{b[0]['blocks']}bl/ex{b[0]['excess']}/sJ{b[0]['min_sj']:.2f}"
                             for b in beam)
            print(f"[beam] depth {d + 1}: {len(cand)} admissible, beam [{tops}]  "
                  f"best {best[0]['blocks']} blocks ({time.time() - t0:.0f}s)")
    st, P, hexes, path = best
    if verbose:
        print(f"[beam] result: {st} via sheets {path}")
    log = [dict(base, sheet=None), dict(st, sheet=path)]
    return P, hexes, log


def patch(width=6, depth=6, guard="abs", sj_floor=0.05, workers=6, hausdorff_tol=0.05):
    def collapse(P, hexes, labeller, max_rounds=4, verbose=True, repair=None,
                 guard_ref=None):
        P, hexes, log = beam_collapse(P, hexes, labeller, width, depth, guard,
                                      sj_floor, workers, guard_ref, verbose,
                                      hausdorff_tol)
        if repair is not None:
            P = repair(P.copy(), hexes)
        return P, hexes, log
    cb.collapse_mesh_sheets = collapse


if __name__ == "__main__":
    import argparse
    import sys
    import tet_prep_v5 as v5
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
    a = ap.parse_args()
    patch(a.width, a.depth, a.guard, a.sj_floor, a.workers, a.hausdorff_tol)
    if a.ovm.endswith(".vtk"):
        # start from a saved blocks.vtk (fine hex mesh AFTER the greedy
        # collapse): BlockStructure otherwise reads .ovm only. The block IDs
        # are discarded -- the blocks follow from the sheets.
        import meshio
        import ovm_io
        _m = meshio.read(a.ovm)
        _H = np.vstack([c.data for c in _m.cells if c.type == "hexahedron"]).astype(np.int64)
        _P = np.asarray(_m.points, float)
        ovm_io.read_ovm = lambda path: (_P, None, None, None)
        ovm_io.ovm_to_cells = lambda P, e, f, poly: (_H, 0)
    S, before, after, v = cb.postprocess(
        a.ovm, a.input_vtk, v5.NAMES, collapse_rounds=1, untangle_mesh=True,
        untangle_rounds=a.untangle_rounds)
    cb.write_blocks(S, a.out)
    print("\nbefore/after:")
    for k in before:
        print(f"  {k:22s} {before[k]}  ->  {after[k]}")
    sys.exit(0)
