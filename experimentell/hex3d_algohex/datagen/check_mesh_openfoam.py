#!/usr/bin/env python3
"""OpenFOAM checkMesh as an ADDITIONAL mesh-quality gate (hex meshes).

    check_mesh_openfoam.py <mesh.vtk|mesh.vtu> [-o report.json]

This is a *second, independent* gate next to the P1 checker
(``datagen/check_keepbl_sample.py``), which stays untouched.  Where the P1
checker works on the exporter's own block bookkeeping, this one hands the mesh
to OpenFOAM's ``checkMesh -allTopology -allGeometry`` and reports what
OpenFOAM thinks of it.

What it does
------------
1. Reads a hexahedral mesh (``.vtk`` / ``.vtu``) with meshio.
2. Writes it out as a native OpenFOAM ``polyMesh`` (points / faces / owner /
   neighbour / boundary) -- no gmsh round-trip, so the vertex order and the
   face orientation are fully under our control.
3. Runs ``checkMesh -allTopology -allGeometry`` in that case.
4. Parses the report and writes a JSON summary.

Pass criteria (exit code 0 only if ALL hold)
--------------------------------------------
  * zero cells with negative (or zero) volume
  * max non-orthogonality <= 65 deg
  * max skewness <= 4

OpenFOAM is v2606 here; the environment is pulled in with
``source /usr/lib/openfoam/openfoam2606/etc/bashrc``.

Exit codes: 0 = PASS, 1 = FAIL (mesh quality), 2 = error (bad input,
OpenFOAM missing, ...).
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------
# configuration
# --------------------------------------------------------------------------

FOAM_BASHRC = os.environ.get(
    "FOAM_BASHRC", "/usr/lib/openfoam/openfoam2606/etc/bashrc"
)
CHECKMESH_ARGS = "-allTopology -allGeometry"

MAX_NON_ORTHOGONALITY = 65.0      # deg
MAX_SKEWNESS = 4.0
VOLUME_EPS = 0.0                  # cells with volume <= this count as bad

# OpenFOAM's hexahedron face model (src/OpenFOAM/lnInclude/hexCell.C).
# The winding of each face is such that the face normal points OUT of the
# cell, which is what OpenFOAM expects for the owner face.
OF_HEX_FACES = (
    (0, 4, 7, 3),      # x-min
    (1, 2, 6, 5),      # x-max
    (0, 1, 5, 4),      # y-min
    (3, 7, 6, 2),      # y-max
    (0, 3, 2, 1),      # z-min
    (4, 5, 6, 7),      # z-max
)

# VTK/meshio and OpenFOAM use the same hexahedron node ordering:
#   0,1,2,3 = base quad, 4,5,6,7 = top quad, same winding sense.
# Verified against checkMesh: a right-handed VTK hex becomes a positive-volume
# OpenFOAM hex without any re-indexing, and an inverted one is reported as a
# negative-volume cell.
VTK_TO_OF_HEX = (0, 1, 2, 3, 4, 5, 6, 7)


class MeshCheckError(RuntimeError):
    """Raised for unusable input or a broken OpenFOAM run."""


# --------------------------------------------------------------------------
# reading the mesh
# --------------------------------------------------------------------------

def read_hex_mesh(path: Path):
    """Return (points, hexes) with hexes in OpenFOAM node order."""
    try:
        import meshio
    except ImportError as exc:  # pragma: no cover - environment problem
        raise MeshCheckError(f"meshio is not available: {exc}") from exc

    try:
        mesh = meshio.read(str(path))
    except Exception as exc:
        raise MeshCheckError(f"cannot read mesh {path}: {exc}") from exc

    hex_blocks = [cb for cb in mesh.cells if cb.type == "hexahedron"]
    if not hex_blocks:
        types = sorted({cb.type for cb in mesh.cells})
        raise MeshCheckError(
            f"{path} contains no hexahedra (cell types found: {types or 'none'})"
        )

    points = np.asarray(mesh.points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise MeshCheckError(f"{path}: expected 3D points, got {points.shape}")

    hexes = np.concatenate([np.asarray(cb.data, dtype=np.int64)
                             for cb in hex_blocks], axis=0)
    if hexes.ndim != 2 or hexes.shape[1] != 8:
        raise MeshCheckError(f"{path}: expected 8 nodes per hex, got {hexes.shape}")

    # VTK -> OpenFOAM node order
    hexes = hexes[:, list(VTK_TO_OF_HEX)]

    # drop cells with an out-of-range node index (defensive: a corrupt file)
    valid = (hexes >= 0).all(axis=1) & (hexes < len(points)).all(axis=1)
    if not valid.all():
        n_bad = int((~valid).sum())
        raise MeshCheckError(
            f"{path}: {n_bad} hexahedra reference points outside the point list"
        )

    return points, hexes


# --------------------------------------------------------------------------
# building the OpenFOAM polyMesh
# --------------------------------------------------------------------------

def _face_normal(points: np.ndarray, face: np.ndarray) -> np.ndarray:
    """Newell normal of a (planar or not) quad face."""
    n = np.zeros(3)
    for i in range(len(face)):
        a = points[face[i]]
        b = points[face[(i + 1) % len(face)]]
        n[0] += (a[1] - b[1]) * (a[2] + b[2])
        n[1] += (a[2] - b[2]) * (a[0] + b[0])
        n[2] += (a[0] - b[0]) * (a[1] + b[1])
    return n


def build_polymesh(points: np.ndarray, hexes: np.ndarray):
    """Build (faces, owner, neighbour) for a pure-hexahedral mesh.

    The winding of every face is taken *verbatim* from the input hexes.
    Nothing is re-wound, because an inverted cell is exactly what this
    gate has to see: if the faces of a cell were silently turned the right
    way, OpenFOAM would no longer notice the negative volume.

    Faces are deduplicated by their *undirected* vertex set, so two cells
    sharing a face collapse into one OpenFOAM face no matter which way
    round they list it.  The owner is the lower cell index and keeps that
    cell's own winding.  For a consistently oriented complex the owner
    face therefore points out of the owner cell, which is what OpenFOAM
    expects.  A cell that is inverted (or a whole block that is mirrored)
    winds its faces the wrong way, and that is precisely what this gate
    must report - so the wrong winding is deliberately kept and left for
    OpenFOAM to see.
    """
    n_cells = len(hexes)
    n_raw = 6 * n_cells

    raw_faces = np.empty((n_raw, 4), dtype=np.int64)
    for fi, f in enumerate(OF_HEX_FACES):
        raw_faces[fi::6] = hexes[:, list(f)]

    # --- group the raw faces by their undirected vertex set -------------
    keys = np.sort(raw_faces, axis=1)
    uniq_keys, inv, counts = np.unique(
        keys, axis=0, return_inverse=True, return_counts=True
    )
    n_faces = len(uniq_keys)

    order = np.argsort(inv, kind="stable")
    slot_start = np.concatenate([[0], np.cumsum(counts)[:-1]])
    slot_end = slot_start + counts

    faces = np.empty((n_faces, 4), dtype=np.int64)
    owner = np.empty(n_faces, dtype=np.int64)
    neighbour = np.full(n_faces, -1, dtype=np.int64)

    for u in range(n_faces):
        slots = order[slot_start[u]:slot_end[u]]
        cells = slots // 6
        if len(cells) > 2:
            raise MeshCheckError(
                f"non-manifold face shared by {len(cells)} cells "
                f"(face {u}); the mesh is not a valid hex complex"
            )

        # The winding is taken verbatim from the input hexes; nothing is
        # re-wound here.  OpenFOAM's convention is that the owner face is
        # wound outwards from the owner cell.  For a consistently oriented
        # complex that holds automatically.  A cell that is inverted (or a
        # whole block that is mirrored) winds its faces the wrong way, and
        # that is precisely what this gate must report - so the wrong
        # winding is deliberately kept and left for OpenFOAM to see.
        if len(cells) == 2:
            c0, c1 = int(cells[0]), int(cells[1])
            s0, s1 = int(slots[0]), int(slots[1])
            own, oth = (c0, c1) if c0 < c1 else (c1, c0)
            own_slot = s0 if own == c0 else s1
            faces[u] = raw_faces[own_slot]
            owner[u] = own
            neighbour[u] = oth
        else:
            faces[u] = raw_faces[int(slots[0])]
            owner[u] = int(cells[0])
            neighbour[u] = -1

    # --- OpenFOAM wants the boundary faces last and contiguous ------------
    # The neighbour list only covers the internal (leading) faces, and the
    # boundary patch starts at the first boundary face.
    is_bnd = neighbour < 0
    reorder = np.nonzero(~is_bnd)[0].tolist() + np.nonzero(is_bnd)[0].tolist()
    faces = faces[reorder]
    owner = owner[reorder]
    neighbour = neighbour[reorder]

    return faces, owner, neighbour


_POLY_HEADER = (
    "/*--------------------------------*- C++ -*"
    "----------------------------------*\\\n"
    "| =========                 |"
    "                                                 |\n"
    "| \\\\      /  F ield         |"
    " OpenFOAM: The Open Source CFD Toolbox           |\n"
    "|  \\\\    /   O peration     |"
    "                                                 |\n"
    "|   \\\\  /    A nd           |"
    "                                                 |\n"
    "|    \\\\/     M anipulation  |"
    "                                                 |\n"
    "\\*--------------------------------------------------------------"
    "-------------*/\n"
    "FoamFile\n"
    "{{\n"
    "    version     2.0;\n"
    "    format      ascii;\n"
    "    class       {cls};\n"
    '    location    "constant/polyMesh";\n'
    "    object      {obj};\n"
    "}}\n"
    "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * "
    "* * * * * * * //\n"
    "\n"
)


def write_case(case: Path, points: np.ndarray, faces: np.ndarray,
               owner: np.ndarray, neighbour: np.ndarray) -> None:
    """Write a minimal but complete OpenFOAM case around the polyMesh."""
    pm = case / "constant" / "polyMesh"
    pm.mkdir(parents=True, exist_ok=True)
    sysdir = case / "system"
    sysdir.mkdir(parents=True, exist_ok=True)

    (sysdir / "controlDict").write_text(
        "FoamFile { version 2.0; format ascii; class dictionary; "
        "object controlDict; }\n"
        "application checkMesh;\n"
        "startTime 0;\n"
        "stopAt endTime;\n"
        "endTime 1;\n"
        "deltaT 1;\n"
        "writeControl timeStep;\n"
        "writeInterval 1;\n"
        "writeFormat ascii;\n"
        "writePrecision 6;\n"
        "writeCompression off;\n"
        "timeFormat general;\n"
        "timePrecision 6;\n"
        "runTimeModifiable true;\n"
    )
    (sysdir / "fvSchemes").write_text(
        "FoamFile { version 2.0; format ascii; class dictionary; "
        "object fvSchemes; }\n"
        "ddtSchemes { default steady; }\n"
        "gradSchemes { default Gauss linear; }\n"
        "divSchemes { default none; }\n"
        "laplacianSchemes { default Gauss linear corrected; }\n"
        "interpolationSchemes { default linear; }\n"
        "snGradSchemes { default corrected; }\n"
    )
    (sysdir / "fvSolution").write_text(
        "FoamFile { version 2.0; format ascii; class dictionary; "
        "object fvSolution; }\n"
        "solvers { }\n"
    )

    with open(pm / "points", "w") as fh:
        fh.write(_POLY_HEADER.format(cls="vectorField", obj="points"))
        fh.write("%d\n(\n" % len(points))
        for p in points:
            fh.write("(%.17g %.17g %.17g)\n" % (p[0], p[1], p[2]))
        fh.write(")\n")

    with open(pm / "faces", "w") as fh:
        fh.write(_POLY_HEADER.format(cls="faceList", obj="faces"))
        fh.write("%d\n(\n" % len(faces))
        for f in faces:
            fh.write("4(%d %d %d %d)\n" % tuple(f))
        fh.write(")\n")

    with open(pm / "owner", "w") as fh:
        fh.write(_POLY_HEADER.format(cls="labelList", obj="owner"))
        fh.write("%d\n(\n" % len(owner))
        for o in owner:
            fh.write("%d\n" % o)
        fh.write(")\n")

    internal = neighbour[neighbour >= 0]
    with open(pm / "neighbour", "w") as fh:
        fh.write(_POLY_HEADER.format(cls="labelList", obj="neighbour"))
        fh.write("%d\n(\n" % len(internal))
        for n in internal:
            fh.write("%d\n" % n)
        fh.write(")\n")

    n_boundary = int((neighbour < 0).sum())
    n_internal = len(faces) - n_boundary
    with open(pm / "boundary", "w") as fh:
        fh.write(_POLY_HEADER.format(cls="polyBoundaryMesh", obj="boundary"))
        fh.write("%d\n(\n" % 1)
        fh.write(
            "    outer\n"
            "    {\n"
            "        type            wall;\n"
            "        inGroups        1 ( wall );\n"
            "        nFaces          %d;\n"
            "        startFace       %d;\n"
            "    }\n"
            ")\n" % (n_boundary, n_internal)
        )


# --------------------------------------------------------------------------
# running checkMesh
# --------------------------------------------------------------------------

def foam_available() -> bool:
    if not Path(FOAM_BASHRC).exists():
        return False
    return shutil.which("bash") is not None


def run_checkmesh(case: Path, timeout: int = 3600):
    """Run checkMesh in `case`; return the combined output text."""
    cmd = (
        f"source {FOAM_BASHRC} >/dev/null 2>&1; "
        f"cd '{case}'; "
        f"checkMesh {CHECKMESH_ARGS} 2>&1"
    )
    try:
        proc = subprocess.run(
            ["bash", "-c", cmd],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired as exc:
        raise MeshCheckError(f"checkMesh timed out after {timeout}s") from exc
    out = (proc.stdout or "") + (proc.stderr or "")
    if "FOAM FATAL" in out:
        raise MeshCheckError("checkMesh failed:\n" + out[-4000:])
    return out


# --------------------------------------------------------------------------
# parsing
# --------------------------------------------------------------------------

_RE_NEG_VOL = re.compile(
    r"Zero or negative cell volume detected\.\s+"
    r"Minimum negative volume:\s*(-?[0-9.eE+-]+),\s*"
    r"Number of negative volume cells:\s*(\d+)",
    re.S,
)
_RE_MIN_VOL = re.compile(
    r"Min volume\s*=\s*(-?[0-9.eE+-]+?)\s*\.?\s*"
    r"Max volume\s*=\s*(-?[0-9.eE+-]+?)\s*\.?\s*"
    r"Total volume\s*=\s*(-?[0-9.eE+-]+?)\s*\.?\s*"
)
_RE_NORTH = re.compile(
    r"Mesh non-orthogonality Max:\s*([0-9.eE+-]+)\s*average:\s*([0-9.eE+-]+)"
)
_RE_NORTH_SEVERE = re.compile(
    r"Number of severely non-orthogonal\s*\(\s*>\s*70 degrees?\s*\)\s*faces:\s*(\d+)"
)
_RE_NORTH_ERR = re.compile(r"Number of non-orthogonality errors:\s*(\d+)")
_RE_SKEW = re.compile(r"Max skewness\s*=\s*([0-9.eE+-]+)")
_RE_SKEW_FACES = re.compile(r"(\d+)\s+highly skew faces detected")
_RE_ASPECT = re.compile(r"Max aspect ratio:?\s*=?\s*([0-9.eE+-]+)")
_RE_DET = re.compile(
    r"Cell determinant \(wellposedness\)\s*:\s*minimum:\s*(-?[0-9.eE+-]+)"
)
_RE_ORIENT = re.compile(
    r"Error in face pyramids:\s*(\d+)\s*faces are incorrectly oriented"
)
_RE_PYRAMID = re.compile(r"pyramids:\s*(\d+)")


def parse_report(text: str) -> dict:
    """Pull the numbers we gate on out of the checkMesh output."""
    out: dict = {}

    m = _RE_NEG_VOL.search(text)
    if m:
        out["negative_volume_cells"] = int(m.group(2))
        out["min_negative_volume"] = float(m.group(1))
    else:
        out["negative_volume_cells"] = 0
        out["min_negative_volume"] = 0.0

    m = _RE_MIN_VOL.search(text)
    if m:
        out["min_volume"] = float(m.group(1))
        out["max_volume"] = float(m.group(2))
        out["total_volume"] = float(m.group(3))
    else:
        out["min_volume"] = None
        out["max_volume"] = None
        out["total_volume"] = None

    m = _RE_NORTH.search(text)
    if m:
        out["max_non_orthogonality"] = float(m.group(1))
        out["average_non_orthogonality"] = float(m.group(2))
    else:
        out["max_non_orthogonality"] = 0.0
        out["average_non_orthogonality"] = 0.0

    m = _RE_NORTH_SEVERE.search(text)
    out["severe_non_orthogonal_faces"] = int(m.group(1)) if m else 0

    m = _RE_NORTH_ERR.search(text)
    out["non_orthogonality_errors"] = int(m.group(1)) if m else 0

    m = _RE_SKEW.search(text)
    out["max_skewness"] = float(m.group(1)) if m else 0.0

    m = _RE_SKEW_FACES.search(text)
    out["highly_skew_faces"] = int(m.group(1)) if m else 0

    m = _RE_ASPECT.search(text)
    out["max_aspect_ratio"] = float(m.group(1)) if m else None

    m = _RE_DET.search(text)
    out["min_cell_determinant"] = float(m.group(1)) if m else None

    m = _RE_ORIENT.search(text)
    out["incorrectly_oriented_faces"] = int(m.group(1)) if m else 0

    m = _RE_PYRAMID.search(text)
    out["pyramid_cells"] = int(m.group(1)) if m else 0

    out["failed_mesh_checks"] = len(
        re.findall(r"^\s*\*\*\*", text, re.M)
    )
    return out


def evaluate(report: dict) -> dict:
    """Apply the pass criteria."""
    reasons = []

    neg = int(report.get("negative_volume_cells", 0))
    if neg > 0:
        reasons.append(
            f"{neg} cell(s) with non-positive volume "
            f"(minimum {report.get('min_negative_volume')})"
        )

    minvol = report.get("min_volume")
    if minvol is not None and minvol <= VOLUME_EPS:
        reasons.append(f"minimum cell volume {minvol} <= {VOLUME_EPS}")

    nonorth = float(report.get("max_non_orthogonality", 0.0))
    if nonorth > MAX_NON_ORTHOGONALITY:
        reasons.append(
            f"max non-orthogonality {nonorth:g} deg > {MAX_NON_ORTHOGONALITY:g} deg"
        )

    skew = float(report.get("max_skewness", 0.0))
    if skew > MAX_SKEWNESS:
        reasons.append(f"max skewness {skew:g} > {MAX_SKEWNESS:g}")

    if int(report.get("incorrectly_oriented_faces", 0)) > 0:
        reasons.append(
            f"{report['incorrectly_oriented_faces']} face(s) with incorrect "
            "orientation"
        )

    return {
        "passed": not reasons,
        "reasons": reasons,
        "criteria": {
            "negative_volume_cells": 0,
            "max_non_orthogonality_deg": MAX_NON_ORTHOGONALITY,
            "max_skewness": MAX_SKEWNESS,
        },
    }


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------

def check_mesh(mesh_path: Path, keep_case: Path | None = None) -> dict:
    """Run the whole check and return the report dict."""
    mesh_path = Path(mesh_path)
    if not mesh_path.is_file():
        raise MeshCheckError(f"mesh file not found: {mesh_path}")

    if not foam_available():
        raise MeshCheckError(
            f"OpenFOAM environment not found (expected {FOAM_BASHRC}); "
            "cannot run checkMesh"
        )

    points, hexes = read_hex_mesh(mesh_path)

    # next to the mesh, not in $TMPDIR: /tmp is a 64 MB tmpfs in the container and a
    # 100k-cell polyMesh case does not fit there ("No space left on device")
    tmpdir = Path(tempfile.mkdtemp(prefix="foamcheck.", dir=Path(mesh_path).resolve().parent))
    try:
        case = tmpdir / "case"
        case.mkdir()
        faces, owner, neighbour = build_polymesh(points, hexes)
        write_case(case, points, faces, owner, neighbour)

        text = run_checkmesh(case)

        report = parse_report(text)
        report["mesh"] = str(mesh_path)
        report["n_points"] = int(len(points))
        report["n_cells"] = int(len(hexes))
        report["n_faces"] = int(len(faces))
        report["foam_bashrc"] = FOAM_BASHRC
        report["checkmesh_args"] = CHECKMESH_ARGS

        verdict = evaluate(report)
        report["passed"] = verdict["passed"]
        report["reasons"] = verdict["reasons"]
        report["criteria"] = verdict["criteria"]

        if keep_case is not None:
            keep_case = Path(keep_case)
            keep_case.parent.mkdir(parents=True, exist_ok=True)
            if keep_case.exists():
                shutil.rmtree(keep_case)
            shutil.copytree(case, keep_case)
            report["case_dir"] = str(keep_case)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    return report


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="OpenFOAM checkMesh gate for hexahedral meshes."
    )
    ap.add_argument("mesh", help="hex mesh (.vtk / .vtu)")
    ap.add_argument("-o", "--out", default=None,
                    help="write the JSON report here (default: stdout only)")
    ap.add_argument("--keep-case", default=None,
                    help="copy the generated OpenFOAM case here")
    ap.add_argument("--quiet", action="store_true",
                    help="do not print the report to stdout")
    args = ap.parse_args(argv)

    try:
        report = check_mesh(Path(args.mesh), keep_case=args.keep_case)
    except MeshCheckError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    if args.out:
        out_path = Path(args.out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    if not args.quiet:
        print(json.dumps(report, indent=2, sort_keys=True))

    if report["passed"]:
        print("PASS: OpenFOAM checkMesh found no mesh-quality problem")
        return 0

    print("FAIL: " + "; ".join(report["reasons"]), file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
