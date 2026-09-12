"""Element histogram of a Gmsh 2.2 ASCII MSH file.

Prints, per element type: total count, count without a physical group
(tag 0), and the count per geometric entity id. Intended as a quick
structural check of a dtOO export before tet_prep_v5 reads it.

The dataset machine meshes must look like the T1_9 / candidate meshes:
volume elements et=4 (tets), et=5 (hexes), et=6 (prisms), et=7 (pyramids)
plus tagged 2D elements et=2/3. A mesh with only et=11 (quadratic tets)
is the wrong boundedVolume: `ruWithRounding_mechMesh` is the structural
(FEM) runner solid, not the fluid grid this pipeline reads.

Usage:
    .venv/bin/python experimentell/hex3d_algohex/scripts/msh_histogram.py FILE.msh
"""

import sys
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from dp3d.extraction import parse_msh  # noqa: E402

ETYPE = {
    1: "line",
    2: "tri",
    3: "quad",
    4: "tet",
    5: "hex",
    6: "prism",
    7: "pyramid",
    8: "line2",
    9: "tri2",
    10: "quad2",
    11: "tet2",
}


def main(path: str) -> None:
    _nodes, elements = parse_msh(path)
    print(f"file  {path}")
    print(f"elements {len(elements)}")

    by_type: Counter[int] = Counter()
    no_phys: Counter[int] = Counter()
    by_geom: Counter[tuple[int, int]] = Counter()
    for etype, phys, geom, _tags in elements:
        by_type[etype] += 1
        if phys == 0:
            no_phys[etype] += 1
        by_geom[(etype, geom)] += 1

    print("\nby element type:")
    for etype in sorted(by_type):
        name = ETYPE.get(etype, f"et{etype}")
        print(f"  et={etype:<3} {name:<8} {by_type[etype]:>8}   phys=0: {no_phys[etype]}")

    print("\nby (etype, geom), geom with >100 elements:")
    for (etype, geom), count in sorted(by_geom.items()):
        if count > 100:
            name = ETYPE.get(etype, f"et{etype}")
            print(f"  et={etype:<3} geom={geom:<4} {name:<8} {count:>8}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        raise SystemExit(2)
    main(sys.argv[1])
