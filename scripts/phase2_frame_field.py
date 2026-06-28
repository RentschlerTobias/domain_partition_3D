#!/usr/bin/env python3
"""
Phase 2: 3D Surface Frame Field

Input:  output/T1_9/hub/hub_phase1_mesh.vtk
Output: output/T1_9/hub/hub_phase2_field.vtk

Algorithm:
1. Load Phase 1 mesh with normals/e1/e2 per facet and boundary flags
2. Compute vertex-local tangent basis (n_v, e1_v, e2_v)
3. Compute boundary conditions: project boundary edges into vertex-local basis
4. Assemble stiffness matrix in vertex-local coordinates (rotate per facet)
5. Apply Dirichlet BC for boundary vertices
6. Solve sparse linear system
7. Normalize per vertex
8. Map field back to 3D for visualization
9. Export VTK with point vectors
"""

import meshio
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
from pathlib import Path
import time

# Paths
IN_VTK = Path("/root/repos/block_structured_meshing/output/T1_9/hub/hub_phase1_mesh.vtk")
OUT_DIR = Path("/root/repos/block_structured_meshing/output/T1_9/hub")
OUT_VTK = OUT_DIR / "hub_phase2_field.vtk"


def load_phase1(vtk_file):
    """Load Phase 1 VTK with points, cells, normals, e1, e2, boundary."""
    mesh = meshio.read(vtk_file)
    points = mesh.points
    cells = mesh.cells[0].data
    normals = mesh.cell_data["normal"][0]
    e1 = mesh.cell_data["e1"][0]
    e2 = mesh.cell_data["e2"][0]
    boundary = mesh.point_data["boundary"]
    
    print(f"Loaded: {points.shape[0]} points, {cells.shape[0]} triangles")
    print(f"Cell data: normal {normals.shape}, e1 {e1.shape}, e2 {e2.shape}")
    print(f"Boundary points: {boundary.sum()}")
    return points, cells, normals, e1, e2, boundary


def compute_vertex_basis(points, cells, normals, e1, e2):
    """
    Compute vertex-local tangent basis.
    n_v = average of adjacent facet normals
    e1_v = projection of global x-axis onto tangent plane
    e2_v = n_v x e1_v
    """
    n_points = len(points)
    n_v = np.zeros((n_points, 3))
    count = np.zeros(n_points)
    
    # Average facet normals at vertices
    for fi, tri in enumerate(cells):
        for vi in tri:
            n_v[vi] += normals[fi]
            count[vi] += 1
    
    for i in range(n_points):
        if count[i] > 0:
            n_v[i] /= count[i]
        norm = np.linalg.norm(n_v[i])
        if norm > 1e-12:
            n_v[i] /= norm
    
    # e1_v = projection of x-axis onto tangent plane
    e1_v = np.zeros((n_points, 3))
    x_axis = np.array([1.0, 0.0, 0.0])
    for i in range(n_points):
        # Project x_axis onto tangent plane
        proj = x_axis - np.dot(x_axis, n_v[i]) * n_v[i]
        norm = np.linalg.norm(proj)
        if norm > 1e-12:
            e1_v[i] = proj / norm
        else:
            # x_axis parallel to normal, use y_axis
            y_axis = np.array([0.0, 1.0, 0.0])
            proj = y_axis - np.dot(y_axis, n_v[i]) * n_v[i]
            norm = np.linalg.norm(proj)
            if norm > 1e-12:
                e1_v[i] = proj / norm
            else:
                e1_v[i] = np.array([0.0, 0.0, 1.0])  # fallback
    
    # e2_v = n_v x e1_v
    e2_v = np.cross(n_v, e1_v)
    for i in range(n_points):
        norm = np.linalg.norm(e2_v[i])
        if norm > 1e-12:
            e2_v[i] /= norm
    
    print(f"Vertex basis computed for {n_points} vertices")
    return n_v, e1_v, e2_v


def compute_boundary_conditions(points, cells, boundary, e1_v, e2_v):
    """
    Compute boundary frame directions.
    For each boundary vertex, find adjacent boundary edges, project into
    vertex-local tangent plane, normalize, average, compute angle.
    Returns: boundary_nodes array, boundary_angles array
    """
    n_points = len(points)
    
    # Find boundary edges
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
    
    # Build edge adjacency per vertex
    vertex_edges = {i: [] for i in range(n_points)}
    for e in boundary_edges:
        vertex_edges[e[0]].append(e)
        vertex_edges[e[1]].append(e)
    
    boundary_nodes = []
    boundary_angles = []
    
    for i in range(n_points):
        if not boundary[i]:
            continue
        
        edges = vertex_edges[i]
        if len(edges) < 2:
            print(f"Warning: boundary vertex {i} has {len(edges)} edges")
            continue
        
        # Get edge vectors (from vertex i)
        edge_vecs = []
        for e in edges:
            if e[0] == i:
                vec = points[e[1]] - points[i]
            else:
                vec = points[e[0]] - points[i]
            edge_vecs.append(vec)
        
        # Project into tangent plane
        proj_vecs = []
        for vec in edge_vecs:
            proj = vec - np.dot(vec, n_v[i]) * n_v[i]
            norm = np.linalg.norm(proj)
            if norm > 1e-12:
                proj_vecs.append(proj / norm)
        
        if len(proj_vecs) < 2:
            # Use just one edge if projection failed
            angle = np.arctan2(np.dot(edge_vecs[0], e2_v[i]), 
                              np.dot(edge_vecs[0], e1_v[i])) % (2 * np.pi)
        else:
            # Average two projected edges
            avg = (proj_vecs[0] + proj_vecs[1]) / 2
            norm = np.linalg.norm(avg)
            if norm > 1e-12:
                avg = avg / norm
            angle = np.arctan2(np.dot(avg, e2_v[i]), 
                              np.dot(avg, e1_v[i])) % (2 * np.pi)
        
        boundary_nodes.append(i)
        boundary_angles.append(angle)
    
    boundary_nodes = np.array(boundary_nodes, dtype=int)
    boundary_angles = np.array(boundary_angles)
    
    print(f"Boundary conditions: {len(boundary_nodes)} vertices, angles in [0, 2π)")
    return boundary_nodes, boundary_angles


def compute_local_stiffness_matrix_2d(coords_2d):
    """
    Compute 2D Laplacian stiffness matrix for triangle.
    coords_2d: (3, 2) array of vertex coordinates in 2D
    Returns: (6, 6) matrix
    """
    x = coords_2d[:, 0]
    y = coords_2d[:, 1]
    
    area = 0.5 * np.abs((x[1] - x[0]) * (y[2] - y[0]) - 
                        (x[2] - x[0]) * (y[1] - y[0]))
    
    if area < 1e-12:
        return np.zeros((6, 6))
    
    b = np.zeros(3)
    c = np.zeros(3)
    for i in range(3):
        j = (i + 1) % 3
        k = (i + 2) % 3
        b[i] = y[j] - y[k]
        c[i] = x[k] - x[j]
    
    grad_N = np.array([[b[i], c[i]] for i in range(3)]) / (2 * area)
    
    A_e = np.zeros((6, 6))
    for i in range(3):
        for j in range(3):
            A_e[2*i, 2*j] = area * np.dot(grad_N[i], grad_N[j])
            A_e[2*i+1, 2*j+1] = area * np.dot(grad_N[i], grad_N[j])
    
    return A_e


def assemble_stiffness_matrix(points, cells, normals, e1, e2, e1_v, e2_v, 
                              boundary_nodes, boundary_angles):
    """
    Assemble stiffness matrix with rotation to vertex-local basis.
    """
    n_points = len(points)
    n_elements = len(cells)
    num_dofs = n_points * 2
    
    rows, cols, vals = [], [], []
    b = np.zeros(num_dofs)
    
    for e in range(n_elements):
        tri = cells[e]
        
        # Project vertices to facet-local 2D coordinates
        coords_2d = np.zeros((3, 2))
        for i in range(3):
            vi = tri[i]
            v = points[vi]
            coords_2d[i, 0] = np.dot(v, e1[e])
            coords_2d[i, 1] = np.dot(v, e2[e])
        
        # Compute local stiffness matrix in facet-local basis
        A_e_local = compute_local_stiffness_matrix_2d(coords_2d)
        
        # Compute rotation matrices from vertex-local to facet-local
        R = np.zeros((6, 6))
        for i in range(3):
            vi = tri[i]
            # R_i: 2x2 rotation from vertex-local to facet-local
            R_i = np.array([
                [np.dot(e1_v[vi], e1[e]), np.dot(e2_v[vi], e1[e])],
                [np.dot(e1_v[vi], e2[e]), np.dot(e2_v[vi], e2[e])]
            ])
            R[2*i:2*i+2, 2*i:2*i+2] = R_i
        
        # Transform to vertex-local basis: A_e = R^T @ A_e_local @ R
        A_e = R.T @ A_e_local @ R
        
        # Assemble
        for i in range(3):
            vi = tri[i]
            for j in range(3):
                vj = tri[j]
                for di in range(2):
                    for dj in range(2):
                        rows.append(2*vi + di)
                        cols.append(2*vj + dj)
                        vals.append(A_e[2*i+di, 2*j+dj])
    
    # Build sparse matrix
    A = coo_matrix((vals, (rows, cols)), shape=(num_dofs, num_dofs)).tocsr().tolil()
    
    # Apply Dirichlet BC
    for idx, bi in enumerate(boundary_nodes):
        angle = boundary_angles[idx]
        # Representation field: (cos 4θ, sin 4θ)
        ref_angle = 4 * angle
        bx = np.cos(ref_angle)
        by = np.sin(ref_angle)
        
        dof_x = 2 * bi
        dof_y = 2 * bi + 1
        
        A.rows[dof_x] = [dof_x]
        A.data[dof_x] = [1.0]
        A.rows[dof_y] = [dof_y]
        A.data[dof_y] = [1.0]
        
        b[dof_x] = bx
        b[dof_y] = by
    
    A_sparse = A.tocsr()
    return A_sparse, b


def solve_frame_field(A, b):
    """Solve sparse linear system."""
    print("Solving sparse system...")
    t0 = time.time()
    u = spsolve(A, b)
    t1 = time.time()
    print(f"Solved in {t1-t0:.2f}s")
    return u


def normalize_field(u, n_points):
    """Normalize field per vertex."""
    u_x = u[0::2]
    u_y = u[1::2]
    
    norms = np.sqrt(u_x**2 + u_y**2)
    norms = np.maximum(norms, 1e-8)
    
    u_x_norm = u_x / norms
    u_y_norm = u_y / norms
    
    return u_x_norm, u_y_norm


def map_to_3d(u_x, u_y, e1_v, e2_v):
    """Map 2D field vectors to 3D."""
    n_points = len(u_x)
    field_3d = np.zeros((n_points, 3))
    
    for i in range(n_points):
        field_3d[i] = u_x[i] * e1_v[i] + u_y[i] * e2_v[i]
    
    return field_3d


def compute_field_angle(u_x, u_y):
    """Compute angle of field in local basis."""
    angles = np.arctan2(u_y, u_x) % (2 * np.pi)
    return angles


def export_field_vtk(points, cells, field_3d, angles, boundary, out_file):
    """Export mesh with frame field vectors."""
    mesh = meshio.Mesh(
        points=points,
        cells=[("triangle", cells)],
        point_data={
            "frame": field_3d,
            "frame_angle": angles,
            "boundary": boundary,
        },
    )
    mesh.write(out_file)
    print(f"Saved field: {out_file}")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("Phase 2: 3D Surface Frame Field")
    print("=" * 60)
    
    # 1. Load Phase 1
    points, cells, normals, e1, e2, boundary = load_phase1(IN_VTK)
    
    # 2. Compute vertex basis
    global n_v
    n_v, e1_v, e2_v = compute_vertex_basis(points, cells, normals, e1, e2)
    
    # 3. Boundary conditions
    boundary_nodes, boundary_angles = compute_boundary_conditions(
        points, cells, boundary, e1_v, e2_v
    )
    
    # 4. Assemble stiffness matrix
    A, b = assemble_stiffness_matrix(
        points, cells, normals, e1, e2, e1_v, e2_v,
        boundary_nodes, boundary_angles
    )
    
    # 5. Solve
    u = solve_frame_field(A, b)
    
    # 6. Normalize
    u_x, u_y = normalize_field(u, len(points))
    
    # 7. Map to 3D
    field_3d = map_to_3d(u_x, u_y, e1_v, e2_v)
    
    # 8. Compute angle
    angles = compute_field_angle(u_x, u_y)
    
    # 9. Export
    export_field_vtk(points, cells, field_3d, angles, boundary, OUT_VTK)
    
    print("\n" + "=" * 60)
    print("Phase 2 complete")
    print("=" * 60)


if __name__ == "__main__":
    main()
