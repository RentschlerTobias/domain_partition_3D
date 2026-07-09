#!/usr/bin/env python3
"""
Verify the output of scripts/extract_ring_surfaces.py.

Reads the three output VTK files (inner, outer, both), counts faces,
computes radius statistics, and checks surface_type values.

Writes a pass/fail report to .omo/evidence/task-2-extract-hub-shroud-from-vtk.txt.
"""

import sys
from pathlib import Path

import numpy as np

# Reuse the parser from the extraction script (do NOT modify it).
sys.path.insert(0, str(Path(__file__).resolve().parent))
from extract_ring_surfaces import parse_vtk, compute_centroid_radii  # noqa: E402


WORKDIR = Path("/home/t1dde/Duty/projects/domain_partition/domain_partition_3D")
EVIDENCE_PATH = WORKDIR / ".omo" / "evidence" / "task-2-extract-hub-shroud-from-vtk.txt"

# Expected values (from inherited wisdom + learnings.md)
EXPECTED = {
    "inner": {"count": 2800, "max_radius": 0.6, "surface_type": 0},
    "outer": {"count": 2800, "min_radius": 1.8, "surface_type": 1},
    "both":  {"count": 5600, "surface_types": {0, 1}},
}

# Hub/shroud band sits between the two rings; any face with radius in this
# range would indicate a leak from the middle band.
HUB_SHROUD_BAND = (0.6, 1.8)


def verify_file(label, path, expected):
    """Verify a single output VTK file. Returns (passed, report_lines)."""
    lines = []
    lines.append(f"=== {label.upper()} RING: {path.name} ===")

    if not path.exists():
        lines.append(f"  FAIL: file not found: {path}")
        return False, lines

    points, cells, cell_types, cell_data = parse_vtk(path)
    n_faces = len(cells)
    radii = compute_centroid_radii(points, cells)

    lines.append(f"  faces: {n_faces} (expected {expected['count']})")
    lines.append(f"  radius: min={radii.min():.4f}, max={radii.max():.4f}, "
                 f"mean={radii.mean():.4f}")

    # Surface type check
    if "surface_type" not in cell_data:
        lines.append("  FAIL: missing surface_type cell data")
        return False, lines

    st = np.asarray(cell_data["surface_type"])
    unique_st = sorted(set(int(v) for v in st))
    lines.append(f"  surface_type values present: {unique_st}")

    checks = []

    # 1. Face count
    checks.append(("face count", n_faces == expected["count"]))

    # 2. Radius range
    if "max_radius" in expected:
        checks.append((f"all radii < {expected['max_radius']}",
                       bool(radii.max() < expected["max_radius"])))
    if "min_radius" in expected:
        checks.append((f"all radii > {expected['min_radius']}",
                       bool(radii.min() > expected["min_radius"])))

    # 3. Surface type values
    if "surface_type" in expected:
        checks.append((f"surface_type == {expected['surface_type']}",
                       unique_st == [expected["surface_type"]]))
    if "surface_types" in expected:
        checks.append((f"surface_type in {sorted(expected['surface_types'])}",
                       set(unique_st) == expected["surface_types"]))

    # 4. Hub/shroud exclusion: no face should fall in the middle band
    in_band = ((radii >= HUB_SHROUD_BAND[0]) & (radii <= HUB_SHROUD_BAND[1])).sum()
    checks.append((f"no faces in hub/shroud band {HUB_SHROUD_BAND}",
                   in_band == 0))
    lines.append(f"  faces in hub/shroud band: {in_band}")

    # 5. Cell type sanity (all quads = 9)
    unique_ct = sorted(set(cell_types))
    lines.append(f"  cell_types present: {unique_ct}")
    checks.append(("all cells are quads (type 9)", unique_ct == [9]))

    passed = all(ok for _, ok in checks)
    for name, ok in checks:
        lines.append(f"  [{'PASS' if ok else 'FAIL'}] {name}")

    lines.append(f"  RESULT: {'PASS' if passed else 'FAIL'}")
    lines.append("")
    return passed, lines


def main():
    files = [
        ("inner", WORKDIR / "T1_9_inner_ring.vtk", EXPECTED["inner"]),
        ("outer", WORKDIR / "T1_9_outer_ring.vtk", EXPECTED["outer"]),
        ("both",  WORKDIR / "T1_9_both_ring.vtk",  EXPECTED["both"]),
    ]

    all_lines = []
    all_lines.append("Verification of scripts/extract_ring_surfaces.py output")
    all_lines.append(f"Generated: {np.datetime64('now')}")
    all_lines.append(f"Workdir:   {WORKDIR}")
    all_lines.append("")

    overall_pass = True
    for label, path, expected in files:
        ok, report = verify_file(label, path, expected)
        all_lines.extend(report)
        overall_pass = overall_pass and ok

    # Cross-check: both = inner + outer (face counts and surface_type counts)
    inner_path = WORKDIR / "T1_9_inner_ring.vtk"
    outer_path = WORKDIR / "T1_9_outer_ring.vtk"
    both_path  = WORKDIR / "T1_9_both_ring.vtk"
    if inner_path.exists() and outer_path.exists() and both_path.exists():
        _, ic, _, idata = parse_vtk(inner_path)
        _, oc, _, odata = parse_vtk(outer_path)
        _, bc, _, bdata = parse_vtk(both_path)
        all_lines.append("=== CROSS-CHECK: both = inner + outer ===")
        all_lines.append(f"  inner faces: {len(ic)}, outer faces: {len(oc)}, "
                         f"both faces: {len(bc)}")
        sum_ok = len(bc) == len(ic) + len(oc)
        all_lines.append(f"  [{'PASS' if sum_ok else 'FAIL'}] "
                         f"both == inner + outer ({len(bc)} == {len(ic)} + {len(oc)})")

        # surface_type counts in both file
        bst = np.asarray(bdata["surface_type"])
        n_inner = int((bst == 0).sum())
        n_outer = int((bst == 1).sum())
        all_lines.append(f"  both surface_type=0 count: {n_inner} (expected {len(ic)})")
        all_lines.append(f"  both surface_type=1 count: {n_outer} (expected {len(oc)})")
        split_ok = (n_inner == len(ic)) and (n_outer == len(oc))
        all_lines.append(f"  [{'PASS' if split_ok else 'FAIL'}] "
                         f"both splits correctly into inner/outer")
        overall_pass = overall_pass and sum_ok and split_ok
        all_lines.append("")

    all_lines.append("=" * 60)
    all_lines.append(f"OVERALL: {'PASS' if overall_pass else 'FAIL'}")
    all_lines.append("=" * 60)

    report_text = "\n".join(all_lines)
    print(report_text)

    EVIDENCE_PATH.parent.mkdir(parents=True, exist_ok=True)
    EVIDENCE_PATH.write_text(report_text + "\n")
    print(f"\nEvidence written to: {EVIDENCE_PATH}")

    return 0 if overall_pass else 1


if __name__ == "__main__":
    sys.exit(main())
