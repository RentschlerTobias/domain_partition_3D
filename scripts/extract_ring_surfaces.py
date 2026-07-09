#!/usr/bin/env python3
"""
Extract inner/outer ring boundary faces from T1_9_hex_boundaries.vtk.

The inner ring is the cylindrical wall around the blade hole (smallest radial
distance from the Z-axis). The outer ring is the outer cylindrical wall
(largest radial distance). Hub/shroud horizontal surfaces sit in the middle
of the radial distribution and are excluded by the percentile-based filter.

Usage:
    python scripts/extract_ring_surfaces.py            # inner ring (default)
    python scripts/extract_ring_surfaces.py --outer    # outer ring only
    python scripts/extract_ring_surfaces.py --both     # inner + outer
"""

import argparse
import sys
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# VTK I/O
# ---------------------------------------------------------------------------

def parse_vtk(filename):
    """Parse ASCII VTK unstructured grid file.

    Returns
    -------
    points : np.ndarray, shape (n_points, 3)
    cells  : list[list[int]]  (one list of node indices per cell)
    cell_types : list[int]
    cell_data : dict[str, list[int]]
    """
    with open(filename, "r") as f:
        lines = f.readlines()

    points = []
    cells = []
    cell_types = []
    cell_data = {}

    i = 0
    n = len(lines)
    while i < n:
        line = lines[i].strip()

        if line.startswith("POINTS"):
            num_points = int(line.split()[1])
            i += 1
            for _ in range(num_points):
                coords = lines[i].split()
                points.append([float(c) for c in coords[:3]])
                i += 1
            continue

        if line.startswith("CELLS"):
            num_cells = int(line.split()[1])
            i += 1
            for _ in range(num_cells):
                tokens = lines[i].split()
                cells.append([int(t) for t in tokens[1:]])
                i += 1
            continue

        if line.startswith("CELL_TYPES"):
            num_types = int(line.split()[1])
            i += 1
            for _ in range(num_types):
                cell_types.append(int(lines[i].strip()))
                i += 1
            continue

        if line.startswith("CELL_DATA"):
            num_data = int(line.split()[1])
            i += 1
            while i < n and lines[i].strip().startswith("SCALARS"):
                data_name = lines[i].split()[1]
                i += 2  # skip LOOKUP_TABLE line
                values = []
                for _ in range(num_data):
                    values.append(int(lines[i].strip()))
                    i += 1
                cell_data[data_name] = values
            continue

        i += 1

    return np.asarray(points, dtype=float), cells, cell_types, cell_data


def write_vtk(filename, points, cells, cell_types, cell_data):
    """Write ASCII VTK unstructured grid file mirroring the input format."""
    with open(filename, "w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("T1_9 Ring Surface Extraction\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")

        # POINTS
        f.write(f"POINTS {len(points)} float\n")
        for p in points:
            f.write(f"{p[0]} {p[1]} {p[2]}\n")

        # CELLS
        total_size = sum(1 + len(c) for c in cells)
        f.write(f"CELLS {len(cells)} {total_size}\n")
        for c in cells:
            f.write(f"{len(c)} " + " ".join(str(n) for n in c) + "\n")

        # CELL_TYPES
        f.write(f"CELL_TYPES {len(cell_types)}\n")
        for t in cell_types:
            f.write(f"{t}\n")

        # CELL_DATA
        f.write(f"CELL_DATA {len(cells)}\n")
        for name, values in cell_data.items():
            f.write(f"SCALARS {name} int 1\n")
            f.write("LOOKUP_TABLE default\n")
            for v in values:
                f.write(f"{v}\n")


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def compute_centroid_radii(points, cells):
    """Return radial distance sqrt(cx^2 + cy^2) of each cell's centroid."""
    centroids = np.zeros((len(cells), 3), dtype=float)
    for idx, cell in enumerate(cells):
        centroids[idx] = points[cell].mean(axis=0)
    return np.sqrt(centroids[:, 0] ** 2 + centroids[:, 1] ** 2)


def select_by_radius(radii, mode, inner_pct=25.0, outer_pct=75.0):
    """Return boolean mask selecting inner / outer / both rings.

    Inner ring = faces whose centroid radius is below the inner_pct percentile.
    Outer ring = faces whose centroid radius is above the outer_pct percentile.
    Hub/shroud horizontal faces sit in the middle band and are excluded.
    """
    inner_thr = np.percentile(radii, inner_pct)
    outer_thr = np.percentile(radii, outer_pct)

    inner_mask = radii <= inner_thr
    outer_mask = radii >= outer_thr

    if mode == "inner":
        return inner_mask, inner_thr, outer_thr
    if mode == "outer":
        return outer_mask, inner_thr, outer_thr
    if mode == "both":
        return inner_mask | outer_mask, inner_thr, outer_thr
    raise ValueError(f"Unknown mode: {mode}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Extract inner/outer ring boundary faces from "
                    "T1_9_hex_boundaries.vtk based on radial distance from Z-axis."
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--inner", action="store_true",
                       help="Extract inner ring only (default).")
    group.add_argument("--outer", action="store_true",
                       help="Extract outer ring only.")
    group.add_argument("--both", action="store_true",
                       help="Extract both inner and outer rings.")
    parser.add_argument("--input", default="T1_9_hex_boundaries.vtk",
                        help="Input VTK file (default: T1_9_hex_boundaries.vtk).")
    parser.add_argument("--output", default=None,
                        help="Output VTK file (default: T1_9_<mode>_ring.vtk).")
    args = parser.parse_args()

    if args.outer:
        mode = "outer"
    elif args.both:
        mode = "both"
    else:
        mode = "inner"

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"ERROR: input file not found: {input_path}", file=sys.stderr)
        return 1

    output_path = Path(args.output) if args.output else Path(f"T1_9_{mode}_ring.vtk")

    print(f"Reading {input_path} ...")
    points, cells, cell_types, cell_data = parse_vtk(input_path)
    print(f"  points={len(points)}, cells={len(cells)}, cell_types={len(cell_types)}")

    radii = compute_centroid_radii(points, cells)
    mask, inner_thr, outer_thr = select_by_radius(radii, mode)

    # Build filtered arrays
    sel_indices = np.where(mask)[0]
    sel_cells = [cells[i] for i in sel_indices]
    sel_types = [cell_types[i] for i in sel_indices]
    sel_data = {name: [values[i] for i in sel_indices]
                for name, values in cell_data.items()}

    # Add surface_type cell data (0 = inner, 1 = outer)
    if mode == "both":
        sel_data["surface_type"] = [
            0 if radii[i] <= inner_thr else 1 for i in sel_indices
        ]
    elif mode == "inner":
        sel_data["surface_type"] = [0] * len(sel_indices)
    else:  # outer
        sel_data["surface_type"] = [1] * len(sel_indices)

    print(f"Mode: {mode}")
    print(f"  inner threshold (p25): {inner_thr:.4f}")
    print(f"  outer threshold (p75): {outer_thr:.4f}")
    print(f"  selected faces: {len(sel_cells)} / {len(cells)}")
    if len(sel_indices) > 0:
        sel_radii = radii[sel_indices]
        print(f"  selected radius range: [{sel_radii.min():.4f}, {sel_radii.max():.4f}]")
        print(f"  selected radius mean:  {sel_radii.mean():.4f}")

    write_vtk(output_path, points, sel_cells, sel_types, sel_data)
    print(f"Wrote {output_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
