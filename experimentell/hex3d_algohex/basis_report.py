"""One scorecard row per candidate basis, for choosing which run to build on.

The v11-vs-v14 decision was made on numbers scattered across three documents,
and the two that decide whether a basis is usable for TFI were missing from
all of them: whether it has a pinch, and whether any direction class is pinned
to a single division. A block one cell across forces its whole class to one
division everywhere -- v14's 41 tiny blocks are not a cosmetic problem, they
are a resolution constraint on the entire complex.

Reports counts and tail quantiles, never shares: this branch has been fooled
at least six times by a share moving the right way while the structure got
worse (HANDOFF.md "Beware counts that look like progress").

    PY basis_report.py --tag v11 [--tag v15 ...] --out output/.../basis_scorecard.md
"""

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "experimentell" / "hex3d_algohex"))

import base_complex as bc                                              # noqa: E402
import clean_blocks as cb                                              # noqa: E402
import showcase as sc                                                  # noqa: E402
import tet_prep_v5 as v5                                               # noqa: E402
import tfi                                                             # noqa: E402

OUT = REPO / "output" / "hex3d_algohex"
DELIV = OUT / "deliverable"
INPUT_VTK = REPO / "data" / "T1_9" / "T1_9_tet_v5.vtk"


def raw_stats(tag, ovm=None):
    """From the AlgoHex output itself: cells, inverted, singular graph."""
    ovm = Path(ovm) if ovm else OUT / f"T1_9_hex_{tag}.ovm"
    if not ovm.exists():
        return {}
    import ovm_io
    P, e, f, poly = ovm_io.read_ovm(str(ovm))
    H, _skipped = ovm_io.ovm_to_cells(P, e, f, poly)
    sj = cb.scaled_jacobians(P, H)
    f2h, e2h = bc.build_topology(H)
    sing, _val, arcs, _bnd = sc._singular_graph(P, H, topo=(f2h, e2h))
    return {"raw_cells": len(H), "raw_inverted": int((sj <= 0).sum()),
            "raw_min_sj": round(float(sj.min()), 4),
            "raw_mean_sj": round(float(sj.mean()), 4),
            "singular_edges": len(sing), "singular_arcs": len(arcs)}


def feature_edges(input_vtk):
    """Line cells of an AlgoHex input mesh = its feature graph."""
    import meshio
    m = meshio.read(str(input_vtk))
    return int(sum(len(b.data) for b in m.cells if b.type == "line"))


def block_stats(tag, blocks=None):
    """From the postprocessed deliverable: the numbers a basis is judged on."""
    path = Path(blocks) if blocks else DELIV / f"T1_9_blocks_{tag}.vtk"
    if not path.exists():
        return {}
    S, bid, lab = cb.read_blocks_vtk(str(path), str(INPUT_VTK))
    cells = S.cells_of()
    sizes = sorted(len(c) for c in cells.values())
    stat = {r: cb.cuboid_status(S.patches(r)) for r in cells}
    sj = cb.scaled_jacobians(S.P, S.hexes)
    v = cb.HexBlockValidator(S)
    hd = v.boundary_hausdorff(lab)
    pinch = cb.detect_pinch(S, verbose=False)

    # what the block structure is worth to TFI: a class pinned to one division
    # holds the whole complex at that resolution
    f2h, _e2h = bc.build_topology(S.hexes)
    lat, missing = tfi.lattices(S.P, S.hexes, np.asarray(bid, int), f2h,
                                verbose=False)
    classes = tfi.direction_classes(lat, f2h, S.hexes, np.asarray(bid, int),
                                    verbose=False)
    pinned = sum(1 for c in classes
                 if min(int(lat[r][0][ax]) for r, ax in c) <= 1)
    return {"blocks": len(cells),
            "non_cuboid": sum(1 for s, _d in stat.values() if s != "cuboid"),
            "pinches": len([p for p in pinch if p["faces"]]),
            "tiny_lt10": sum(1 for s in sizes if s < 10),
            "smallest_block": sizes[0] if sizes else 0,
            "cells": len(S.hexes),
            "inverted": int((sj <= 0).sum()),
            "min_sj": round(float(sj.min()), 4),
            "sj_p5": round(float(np.percentile(sj, 5)), 4),
            "hausdorff": round(float(hd), 4) if hd is not None else None,
            "validator": "VALID" if v.is_valid() else "INVALID",
            "tfi_classes": len(classes),
            "tfi_classes_at_1": pinned,
            "tfi_non_lattice": len(missing)}


COLUMNS = [("feature_edges", "feat edges"), ("raw_cells", "raw cells"),
           ("raw_inverted", "raw inv"), ("raw_min_sj", "raw min sJ"),
           ("raw_mean_sj", "raw mean sJ"), ("singular_edges", "sing edges"),
           ("singular_arcs", "sing arcs"), ("blocks", "blocks"),
           ("non_cuboid", "non-cuboid"), ("pinches", "pinches"),
           ("tiny_lt10", "tiny <10"), ("smallest_block", "smallest"),
           ("cells", "cells"), ("inverted", "inverted"),
           ("min_sj", "min sJ"), ("sj_p5", "sJ p5"),
           ("hausdorff", "hausdorff"), ("validator", "validator"),
           ("tfi_classes", "TFI classes"), ("tfi_classes_at_1", "classes at 1"),
           ("tfi_non_lattice", "non-lattice")]


def rank(rows):
    """Lexicographic selection rule, as agreed before any run was made:
    0 inverted (hard), then no non-cuboid and no pinch, then no tiny block,
    then fewest blocks, then Hausdorff, then no class pinned to one division."""
    def key(t):
        d = rows[t]
        return (d.get("inverted", 9e9) > 0,
                d.get("non_cuboid", 9e9) + d.get("pinches", 9e9),
                d.get("tiny_lt10", 9e9),
                d.get("blocks", 9e9),
                d.get("hausdorff") or 9e9,
                d.get("tfi_classes_at_1", 9e9))
    return sorted(rows, key=key)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", action="append", default=[])
    ap.add_argument("--map", default=None,
                    help="JSON {tag: {ovm, blocks, input}} for bases whose "
                         "files are not at the conventional paths -- a "
                         "re-postprocessed variant must NOT overwrite the "
                         "deliverable it is being compared against")
    ap.add_argument("--out", default=str(OUT / "basis_scorecard.md"))
    a = ap.parse_args()
    amap = json.loads(Path(a.map).read_text()) if a.map else {}
    tags = a.tag + [t for t in amap if t not in a.tag]

    rows = {}
    for tag in tags:
        print(f"[basis_report] {tag} ...")
        spec = amap.get(tag, {})
        d = {}
        d.update(raw_stats(tag, spec.get("ovm")))
        d.update(block_stats(tag, spec.get("blocks")))
        if spec.get("input"):
            d["feature_edges"] = feature_edges(REPO / spec["input"])
        rows[tag] = d
        print(f"[basis_report]   {d}")

    order = rank(rows)
    lines = ["| metric | " + " | ".join(order) + " |",
             "|---|" + "---|" * len(order)]
    for key, label in COLUMNS:
        vals = [str(rows[t].get(key, "-")) for t in order]
        if all(v == "-" for v in vals):
            continue
        lines.append(f"| {label} | " + " | ".join(vals) + " |")
    lines.append("")
    lines.append(f"Ranked best first by the rule fixed before the runs: "
                 f"0 inverted (hard), then non-cuboid + pinches, then blocks "
                 f"under 10 cells, then block count, then Hausdorff, then "
                 f"direction classes pinned to one division. Winner by that "
                 f"rule: **{order[0]}**.")
    Path(a.out).write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\n[basis_report] wrote {a.out}")


if __name__ == "__main__":
    main()
