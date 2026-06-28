#!/usr/bin/env python3
"""
Phase 1: 3D Surface Mesh Loading + Visualisierung

Input:  T1_9_hub_raw.stl (raw triangulated surface)
Output: output/T1_9/hub/hub_phase1_mesh.vtk
        output/T1_9/hub/hub_phase1_boundary.vtk

Compute:
- Facet normals (recomputed, consistently oriented)
- Tangential basis per facet (e1, e2)
- Boundary edges (count == 1)
- Boundary nodes (unique endpoints)
"""

import meshio
import numpy as np
from pathlib import Path
import sys

# Paths
STL_FILE = Path("/root/repos/block_structured_meshing/T1_9_hub_raw.stl")
OUT_DIR = Path("/root/repos/block_structured_meshing/output/T1_9/hub")


def load_stl(stl_file):
    """Load STL with meshio."""
    mesh = meshio.read(stl_file)
    points = mesh.points
    cells = mesh.cells[0].data
    print(f"Loaded: {points.shape[0]} points, {cells.shape[0]} triangles")
    return points, cells


def compute_facet_normals(points, cells):
    """
    Recompute normals. Orient consistently via flood-fill from seed facet.
    Returns: normals (n_tris, 3), e1 (n_tris, 3), e2 (n_tris, 3)
    """
    n_tris = cells.shape[0]
    normals = np.zeros((n_tris, 3))
    e1 = np.zeros((n_tris, 3))
    e2 = np.zeros((n_tris, 3))

    # Compute raw normals from vertex order
    for i, tri in enumerate(cells):
        v0 = points[tri[0]]
        v1 = points[tri[1]]
        v2 = points[tri[2]]
        e01 = v1 - v0
        e02 = v2 - v0
        n = np.cross(e01, e02)
        norm = np.linalg.norm(n)
        if norm > 1e-12:
            n = n / norm
        normals[i] = n
        # e1 = first edge, e2 = n x e1
        e1_i = e01 / (np.linalg.norm(e01) + 1e-12)
        e2_i = np.cross(n, e1_i)
        e2_i = e2_i / (np.linalg.norm(e2_i) + 1e-12)
        e1[i] = e1_i
        e2[i] = e2_i

    # Build edge -> facet adjacency
    edge_to_facets = {}
    for fi, tri in enumerate(cells):
        edges = [
            (int(min(tri[0], tri[1])), int(max(tri[0], tri[1]))),
            (int(min(tri[1], tri[2])), int(max(tri[1], tri[2]))),
            (int(min(tri[2], tri[0])), int(max(tri[2], tri[0]))),
        ]
        for e in edges:
            if e not in edge_to_facets:
                edge_to_facets[e] = []
            edge_to_facets[e].append(fi)

    # Consistent orientation via BFS
    oriented = np.zeros(n_tris, dtype=bool)
    queue = [0]
    oriented[0] = True
    head = 0

    while head < len(queue):
        fi = queue[head]
        head += 1
        tri = cells[fi]
        edges = [
            (int(min(tri[0], tri[1])), int(max(tri[0], tri[1]))),
            (int(min(tri[1], tri[2])), int(max(tri[1], tri[2]))),
            (int(min(tri[2], tri[0])), int(max(tri[2], tri[0]))),
        ]
        for e in edges:
            for fj in edge_to_facets[e]:
                if fj == fi or oriented[fj]:
                    continue
                # Check if shared edge has same orientation
                tri_j = cells[fj]
                # Find edge in tri_j
                e_set = set(e)
                for k in range(3):
                    a = tri_j[k]
                    b = tri_j[(k + 1) % 3]
                    if set([a, b]) == e_set:
                        # Same vertex order -> flip needed
                        normals[fj] = -normals[fj]
                        e1[fj] = -e1[fj]  # e2 = n x e1 stays consistent
                        break
                oriented[fj] = True
                queue.append(fj)

    print(f"Oriented {oriented.sum()}/{n_tris} facets")
    return normals, e1, e2


def find_boundary_edges(cells):
    """Find edges that belong to exactly one facet."""
    edge_counts = {}
    for tri in cells:
        edges = [
            (int(min(tri[0], tri[1])), int(max(tri[0], tri[1]))),
            (int(min(tri[1], tri[2])), int(max(tri[1], tri[2]))),
            (int(min(tri[2], tri[0])), int(max(tri[2], tri[0]))),
        ]
        for e in edges:
            edge_counts[e] = edge_counts.get(e, 0) + 1

    boundary_edges = [e for e, c in edge_counts.items() if c == 1]
    print(f"Boundary edges: {len(boundary_edges)} (interior: {sum(1 for c in edge_counts.values() if c == 2)})")
    return np.array(boundary_edges)


def export_mesh_vtk(points, cells, normals, e1, e2, boundary_nodes, out_file):
    """Export mesh with cell vectors + point boundary flag."""
    boundary_flag = np.zeros(len(points), dtype=int)
    boundary_flag[boundary_nodes] = 1

    mesh = meshio.Mesh(
        points=points,
        cells=[("triangle", cells)],
        cell_data={
            "normal": [normals],
            "e1": [e1],
            "e2": [e2],
        },
        point_data={
            "boundary": boundary_flag,
        },
    )
    mesh.write(out_file)
    print(f"Saved mesh: {out_file}")


def export_boundary_vtk(points, boundary_edges, out_file):
    """Export boundary edges as lines."""
    lines = boundary_edges.reshape(-1, 2)
    mesh = meshio.Mesh(
        points=points,
        cells=[("line", lines)],
    )
    mesh.write(out_file)
    print(f"Saved boundary: {out_file}")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Phase 1: Surface Mesh Loading")
    print("=" * 60)

    # 1. Load
    points, cells = load_stl(STL_FILE)

    # 2. Normals + tangential basis
    normals, e1, e2 = compute_facet_normals(points, cells)

    # 3. Boundary edges
    boundary_edges = find_boundary_edges(cells)
    boundary_nodes = np.unique(boundary_edges.flatten())
    print(f"Boundary nodes: {len(boundary_nodes)}")

    # 4. Export
    export_mesh_vtk(
        points, cells, normals, e1, e2, boundary_nodes,
        OUT_DIR / "hub_phase1_mesh.vtk"
    )
    export_boundary_vtk(
        points, boundary_edges,
        OUT_DIR / "hub_phase1_boundary.vtk"
    )

    print("\n" + "=" * 60)
    print("Phase 1 complete")
    print("=" * 60)


if __name__ == "__main__":
    main()
