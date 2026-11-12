"""B7: the B5/B6 mesh as an OpenFOAM case with the dtOO patches, plus checkMesh.

    export_foam.py <mesh.vtk> --case <dir> [--json report.json]

Input: the hexahedra and the boundary quads of mesh_b5.py (cell data "patch" on the
quads = index into mesh_b5.PATCHES: RU_HUB, RU_SHROUD, RU_BLADE, RU_INLET, RU_OUTLET,
RU_PERIOPS, RU_PERIOSS -- the names dtOO gives them, so the dtOO OpenFOAM setup
(xml/tistos_ru_of.xml: walls, mapped inlet, cyclicAMI pair, function objects) applies
unchanged). Every boundary face of the hex complex must carry exactly one patch.

Patch types as in dtOO: hub / shroud / blade wall, inlet / outlet patch, the periodic
pair cyclicAMI rotating about the machine axis (z).

The polyMesh is built by datagen/check_mesh_openfoam.build_polymesh (face winding
verbatim, so an inverted cell stays visible) and checked with checkMesh
-allTopology -allGeometry against the same gate (no negative volumes,
non-orthogonality <= 70 deg, skewness <= 4).
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "datagen"))

import check_mesh_openfoam as cmo                                     # noqa: E402
sys.path.insert(0, str(HERE))
from mesh_b5 import PATCHES                                           # noqa: E402

PATCH_TYPE = {"RU_HUB": "wall", "RU_SHROUD": "wall", "RU_BLADE": "wall",
              "RU_INLET": "patch", "RU_OUTLET": "patch",
              "RU_PERIOPS": "cyclicAMI", "RU_PERIOSS": "cyclicAMI"}
AMI_PARTNER = {"RU_PERIOPS": "RU_PERIOSS", "RU_PERIOSS": "RU_PERIOPS"}


def read(path):
    import meshio
    m = meshio.read(str(path))
    H = m.cells_dict["hexahedron"].astype(np.int64)
    qi = [i for i, c in enumerate(m.cells) if c.type == "quad"][0]
    Q = m.cells[qi].data.astype(np.int64)
    pid = np.asarray(m.cell_data["patch"][qi]).astype(int)
    names = dict(enumerate(PATCHES.values()))          # legacy VTK drops field_data
    return np.asarray(m.points, float), H, Q, pid, names


def write_boundary(pm, patches, n_internal):
    """patches: [(name, nFaces)] in face order after the internal faces."""
    out = [cmo._POLY_HEADER.format(cls="polyBoundaryMesh", obj="boundary"), f"{len(patches)}\n(\n"]
    start = n_internal
    for name, nf in patches:
        t = PATCH_TYPE[name]
        out.append(f"    {name}\n    {{\n        type            {t};\n")
        if t == "wall":
            out.append("        inGroups        1 ( wall );\n")
        if t == "cyclicAMI":
            out.append("        inGroups        1 ( cyclicAMI );\n"
                       f"        neighbourPatch  {AMI_PARTNER[name]};\n"
                       "        transform       rotational;\n"
                       "        rotationAxis    (0 0 1);\n"
                       "        rotationCentre  (0 0 0);\n")
        out.append(f"        nFaces          {nf};\n        startFace       {start};\n    }}\n")
        start += nf
    out.append(")\n")
    (pm / "boundary").write_text("".join(out))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("mesh")
    ap.add_argument("--case", required=True, help="OpenFOAM case directory to write")
    ap.add_argument("--json", default=None)
    ap.add_argument("--no-check", action="store_true")
    a = ap.parse_args()
    P, H, Q, pid, names = read(a.mesh)
    faces, owner, neigh = cmo.build_polymesh(P, H)
    nb = int((neigh < 0).sum())
    ni = len(faces) - nb
    # internal faces in upper-triangular order (owner, then neighbour), as OpenFOAM wants
    io = np.lexsort((neigh[:ni], owner[:ni]))
    faces[:ni], owner[:ni], neigh[:ni] = faces[:ni][io], owner[:ni][io], neigh[:ni][io]
    # patch of every boundary face, by its vertex set
    qkey = {tuple(sorted(q)): int(p) for q, p in zip(Q.tolist(), pid.tolist())}
    bf = faces[ni:]
    bp = np.array([qkey.get(tuple(sorted(f)), -1) for f in bf.tolist()])
    missing = int((bp < 0).sum())
    if missing or len(qkey) != nb:
        raise SystemExit(f"[b7] {missing} boundary faces without patch, {len(qkey)} patch quads "
                         f"for {nb} boundary faces")
    order = np.argsort(bp, kind="stable")
    faces[ni:], owner[ni:] = bf[order], owner[ni:][order]
    neigh[ni:] = -1
    case = Path(a.case)
    cmo.write_case(case, P, faces, owner, neigh)
    patches = [(names[k], int((bp == k).sum())) for k in sorted(set(bp.tolist()))]
    write_boundary(case / "constant" / "polyMesh", patches, ni)
    rep = {"mesh": str(a.mesh), "case": str(case), "cells": int(len(H)), "points": int(len(P)),
           "faces": int(len(faces)), "patches": dict(patches)}
    if not a.no_check:
        text = cmo.run_checkmesh(case)
        (case / "checkMesh.log").write_text(text)
        r = cmo.parse_report(text)
        v = cmo.evaluate(r)
        rep.update({"checkmesh": r, "passed": v["passed"], "reasons": v["reasons"]})
    if a.json:
        Path(a.json).write_text(json.dumps(rep, indent=1))
    print(json.dumps({k: rep[k] for k in rep if k != "checkmesh"}, indent=1))
    if "checkmesh" in rep:
        print(json.dumps(rep["checkmesh"], indent=1))
    return 0 if rep.get("passed", True) else 1


if __name__ == "__main__":
    sys.exit(main())
