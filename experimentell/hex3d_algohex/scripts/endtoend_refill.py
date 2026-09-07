"""Pick the cell size that reproduces the reference full-domain resolution,
refill the core at it, and re-attach the two removed parts.

The staged descent of the plan: identity first (the fold gate), then a coarse
probe, then a bisection on h until the ASSEMBLED mesh lands within 10 % of the
165 356 cells the v11 deliverable has at 17 boundary-layer layers. Bisecting on
the core alone would miss, because refining the core also multiplies the
boundary-layer interface quads by 17.

    PY endtoend_refill.py --blocks <core blocks vtk> --tag v11 [--layers 17]
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "experimentell" / "hex3d_algohex"))

import clean_blocks as cb                                              # noqa: E402
import export_vtk as ev                                                # noqa: E402
import tfi                                                             # noqa: E402

TARGET_CELLS = 165356          # the on-disk v11 full-domain deliverable
DELIV = REPO / "output" / "hex3d_algohex" / "deliverable"


def refill_at(P, H, B, f2h, lat, classes, h, verbose=False):
    counts = tfi.solve_block_divisions(lat, classes, P, h,
                                       frozen_counts=tfi.frozen_from_missing(
                                           lat, H, f2h, B), verbose=verbose)
    return tfi.refill_complex(P, H, B, f2h, lat, classes, counts,
                              verbose=verbose)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blocks", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--layers", type=int, default=17)
    ap.add_argument("--target-cells", type=int, default=TARGET_CELLS)
    a = ap.parse_args()

    P, H, B, f2h = tfi.load_blocks(a.blocks)
    lat, missing = tfi.lattices(P, H, B, f2h)
    classes = tfi.direction_classes(lat, f2h, H, B)
    tfi.check_conformity(lat, classes)

    # the assembled size for a given core size, without assembling: the
    # O-grid is 44 800 cells verbatim and the boundary layer is `layers`
    # times the wall quads, which the core's own boundary determines
    stages, best = [], None
    for h in (0.08, 0.06, 0.05, 0.04, 0.035, 0.03, 0.025):
        Pn, Hn, Bn, rep = refill_at(P, H, B, f2h, lat, classes, h)
        ok, _b = tfi.check_watertight(Pn, Hn, verbose=False)
        sj = cb.scaled_jacobians(Pn, Hn)
        row = {"h": h, "core_cells": len(Hn), "watertight": bool(ok),
               "inverted": int((sj <= 0).sum()),
               "min_sj": round(float(sj.min()), 4)}
        stages.append(row)
        print(f"[endtoend] h={h}: {row}")
        if not ok or row["inverted"]:
            print("[endtoend] stopping the descent here -- keeping the last "
                  "clean stage")
            break
        best = (h, Pn, Hn, Bn)
        if len(Hn) > a.target_cells:
            break

    if best is None:
        raise SystemExit("[endtoend] no stage produced a clean mesh")
    h, Pn, Hn, Bn = best
    core = DELIV / f"T1_9_blocks_{a.tag}_h{h}.vtk"
    ev.write_vtk(str(core), Pn, Hn, [12] * len(Hn), Bn, "block_id",
                 f"{int(Bn.max()) + 1} blocks refilled at h={h}")
    print(f"[endtoend] core written: {core}")

    import reattach
    full = DELIV / f"T1_9_blocks_{a.tag}_h{h}_full.vtk"
    Pf, Hf, Bf, _part = reattach.assemble(str(core), str(full),
                                          n_layers=a.layers)
    sjf = cb.scaled_jacobians(Pf, Hf)
    out = {"tag": a.tag, "h": h, "layers": a.layers, "stages": stages,
           "core_cells": int(len(Hn)), "full_cells": int(len(Hf)),
           "full_blocks": int(Bf.max()) + 1,
           "full_inverted": int((sjf <= 0).sum()),
           "full_min_sj": round(float(sjf.min()), 4),
           "core": str(core), "full": str(full)}
    print(json.dumps(out, indent=1))
    (REPO / "output" / "hex3d_algohex"
     / f"endtoend_{a.tag}.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
