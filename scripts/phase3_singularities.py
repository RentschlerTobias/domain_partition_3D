#!/usr/bin/env python3
"""
Phase 3: Singularity Detection + Separatrix Tracing

Input:  output/T1_9/hub/hub_phase2_field.vtk (frame field)
        output/T1_9/hub/hub_phase1_mesh.vtk (cell normals/e1/e2)
Output: output/T1_9/hub/hub_phase3_singularities.vtk
        output/T1_9/hub/hub_phase3_separatrices.vtk

Algorithm:
1. Load frame field from Phase 2
2. Compute vertex-local basis (recompute from Phase 1 cell normals)
3. Detect singularities per triangle with parallel transport
4. Compute singularity locations (barycentric interpolation)
5. Trace separatrices from singularities (RK2 on surface)
6. Export singularities as points + separatrices as lines

BUG FIX: Use correct tangent basis e1/e2 from Phase 1 instead of recomputing with np.cross.
"""

import meshio
import numpy as np
from pathlib import Path
from collections import defaultdict

# Paths
IN_VTK = Path("/root/repos/block_structured_meshing/output/T1_9/hub/hub_phase2_field.vtk")
PHASE1_VTK = Path("/root/repos/block_structured_meshing/output/T1_9/hub/hub_phase1_mesh.vtk")
OUT_DIR = Path("/root/repos/block_structured_meshing/output/T1_9/hub")
OUT_SING = OUT_DIR / "hub_phase3_singularities.vtk"
OUT_SEP = OUT_DIR / "hub_phase3_separatrices.vtk"


def load_phase2(vtk_file):
    """Load Phase 2 VTK with frame field."""
    mesh = meshio.read(vtk_file)
    points = mesh.points
    cells = mesh.cells[0].data
    frame = mesh.point_data["frame"]
    frame_angle = mesh.point_data["frame_angle"]
    boundary = mesh.point_data["boundary"]
    
    u_x = np.cos(frame_angle)
    u_y = np.sin(frame_angle)
    
    print(f"Loaded Phase 2: {points.shape[0]} points, {cells.shape[0]} triangles")
    return points, cells, frame, frame_angle, boundary, u_x, u_y


def load_phase1(vtk_file):
    """Load Phase 1 VTK for cell normals/e1/e2."""
    mesh = meshio.read(vtk_file)
    normals = mesh.cell_data["normal"][0]
    e1 = mesh.cell_data["e1"][0]
    e2 = mesh.cell_data["e2"][0]
    return normals, e1, e2


def compute_vertex_basis(points, cells, normals):
    """
    Compute vertex-local orthonormal basis.
    n_v = average normal of adjacent facets
    e1_v = projection of global x-axis onto tangent plane
    e2_v = n_v x e1_v
    """
    n_points = len(points)
    n_v = np.zeros((n_points, 3))
    count = np.zeros(n_points)
    
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
    
    e1_v = np.zeros((n_points, 3))
    x_axis = np.array([1.0, 0.0, 0.0])
    for i in range(n_points):
        proj = x_axis - np.dot(x_axis, n_v[i]) * n_v[i]
        norm = np.linalg.norm(proj)
        if norm > 1e-12:
            e1_v[i] = proj / norm
        else:
            y_axis = np.array([0.0, 1.0, 0.0])
            proj = y_axis - np.dot(y_axis, n_v[i]) * n_v[i]
            norm = np.linalg.norm(proj)
            if norm > 1e-12:
                e1_v[i] = proj / norm
            else:
                e1_v[i] = np.array([0.0, 0.0, 1.0])
    
    e2_v = np.cross(n_v, e1_v)
    for i in range(n_points):
        norm = np.linalg.norm(e2_v[i])
        if norm > 1e-12:
            e2_v[i] /= norm
    
    return n_v, e1_v, e2_v


def transport_angle(points, e1_v, e2_v, i, j):
    """
    Compute parallel transport angle from vertex i to j.
    """
    e = points[j] - points[i]
    norm = np.linalg.norm(e)
    if norm < 1e-12:
        return 0.0
    e = e / norm
    
    e_i = np.array([np.dot(e, e1_v[i]), np.dot(e, e2_v[i])])
    e_j = np.array([np.dot(e, e1_v[j]), np.dot(e, e2_v[j])])
    
    alpha_i = np.arctan2(e_i[1], e_i[0])
    alpha_j = np.arctan2(e_j[1], e_j[0])
    
    return alpha_j - alpha_i


def detect_singularities(points, cells, u_x, u_y, e1_v, e2_v):
    """
    Detect singularities per triangle with parallel transport.
    For each triangle (a,b,c):
    - Compute representation vector angles at vertices
    - Compute parallel transport angles along edges
    - Sum angle differences (with transport) around triangle
    - Index = round(sum / (2π))
    """
    n_tris = len(cells)
    singularities = []
    
    for fi, tri in enumerate(cells):
        a, b, c = tri
        
        phi_a = np.arctan2(u_y[a], u_x[a])
        phi_b = np.arctan2(u_y[b], u_x[b])
        phi_c = np.arctan2(u_y[c], u_x[c])
        
        # Parallel transport angles
        k_ab = transport_angle(points, e1_v, e2_v, a, b)
        k_bc = transport_angle(points, e1_v, e2_v, b, c)
        k_ca = transport_angle(points, e1_v, e2_v, c, a)
        
        # For cross field, representation angle φ = 4θ.
        # Parallel transport angle must be multiplied by 4.
        # Angle differences with transport, mapped to [-π, π]
        def wrap_diff(diff):
            return ((diff + np.pi) % (2 * np.pi)) - np.pi
        
        d_ab = wrap_diff((phi_b + 4*k_ab) - phi_a)
        d_bc = wrap_diff((phi_c + 4*k_bc) - phi_b)
        d_ca = wrap_diff((phi_a + 4*k_ca) - phi_c)
        
        sum_diff = d_ab + d_bc + d_ca
        index = round(sum_diff / (2 * np.pi))
        
        if index != 0:
            # Compute singularity location in triangle
            e1_f = np.cross(points[b] - points[a], points[c] - points[a])
            norm = np.linalg.norm(e1_f)
            if norm < 1e-12:
                continue
            e1_f /= norm
            
            e2_f = np.cross(e1_f, np.array([0, 0, 1]))
            norm = np.linalg.norm(e2_f)
            if norm < 1e-6:
                e2_f = np.cross(e1_f, np.array([0, 1, 0]))
                norm = np.linalg.norm(e2_f)
            e2_f /= norm
            
            def to_facet_basis(vi):
                u_3d = u_x[vi] * e1_v[vi] + u_y[vi] * e2_v[vi]
                return np.array([np.dot(u_3d, e1_f), np.dot(u_3d, e2_f)])
            
            u_a_f = to_facet_basis(a)
            u_b_f = to_facet_basis(b)
            u_c_f = to_facet_basis(c)
            
            A = np.stack([u_a_f - u_c_f, u_b_f - u_c_f], axis=1)
            b_vec = -u_c_f
            
            try:
                coeffs = np.linalg.solve(A, b_vec)
                if coeffs[0] >= -1e-3 and coeffs[1] >= -1e-3 and coeffs[0] + coeffs[1] <= 1 + 1e-3:
                    loc = points[c] + coeffs[0] * (points[a] - points[c]) + coeffs[1] * (points[b] - points[c])
                    singularities.append({
                        'face': fi,
                        'index': int(index),
                        'location': loc,
                        'vertices': [a, b, c]
                    })
            except np.linalg.LinAlgError:
                pass
    
    print(f"Detected {len(singularities)} singularities")
    for s in singularities:
        print(f"  Face {s['face']}: index={s['index']}, location=({s['location'][0]:.4f}, {s['location'][1]:.4f}, {s['location'][2]:.4f})")
    
    return singularities


def build_mesh_graph(points, cells):
    """Build vertex-to-faces and edge-to-faces adjacency."""
    vertex_to_faces = defaultdict(list)
    for fi, tri in enumerate(cells):
        for vi in tri:
            vertex_to_faces[vi].append(fi)
    
    edge_to_faces = {}
    for fi, tri in enumerate(cells):
        edges = [
            (int(min(tri[0], tri[1])), int(max(tri[0], tri[1]))),
            (int(min(tri[1], tri[2])), int(max(tri[1], tri[2]))),
            (int(min(tri[2], tri[0])), int(max(tri[2], tri[0]))),
        ]
        for e in edges:
            if e not in edge_to_faces:
                edge_to_faces[e] = []
            edge_to_faces[e].append(fi)
    
    return vertex_to_faces, edge_to_faces


def find_containing_face(point, points, cells, vertex_to_faces):
    """Find face containing point using barycentric coordinates."""
    dists = np.linalg.norm(points - point, axis=1)
    closest = np.argmin(dists)
    
    for fi in vertex_to_faces[closest]:
        tri = cells[fi]
        v0 = points[tri[2]] - points[tri[0]]
        v1 = points[tri[1]] - points[tri[0]]
        v2 = point - points[tri[0]]
        
        d00 = np.dot(v0, v0)
        d01 = np.dot(v0, v1)
        d11 = np.dot(v1, v1)
        d20 = np.dot(v2, v0)
        d21 = np.dot(v2, v1)
        
        denom = d00 * d11 - d01 * d01
        if abs(denom) < 1e-12:
            continue
        
        v = (d11 * d20 - d01 * d21) / denom
        w = (d00 * d21 - d01 * d20) / denom
        u = 1 - v - w
        
        if u >= -1e-6 and v >= -1e-6 and w >= -1e-6:
            return fi
    
    return None


def get_field_directions(point, fi, points, cells, u_x, u_y):
    """
    Get the 4 cross field directions at a point in a facet.
    Returns list of 2D unit vectors in facet's local basis.
    """
    if fi is None:
        return None
    
    tri = cells[fi]
    v0 = points[tri[2]] - points[tri[0]]
    v1 = points[tri[1]] - points[tri[0]]
    v2 = point - points[tri[0]]
    
    d00 = np.dot(v0, v0)
    d01 = np.dot(v0, v1)
    d11 = np.dot(v1, v1)
    d20 = np.dot(v2, v0)
    d21 = np.dot(v2, v1)
    
    denom = d00 * d11 - d01 * d01
    if abs(denom) < 1e-12:
        return None
    
    v = (d11 * d20 - d01 * d21) / denom
    w = (d00 * d21 - d01 * d20) / denom
    u = 1 - v - w
    
    u_point = u * np.array([u_x[tri[0]], u_y[tri[0]]]) + \
              v * np.array([u_x[tri[1]], u_y[tri[1]]]) + \
              w * np.array([u_x[tri[2]], u_y[tri[2]]])
    
    norm = np.linalg.norm(u_point)
    if norm < 1e-12:
        return None
    
    u_point = u_point / norm
    
    base_angle = np.arctan2(u_point[1], u_point[0]) / 4
    directions = []
    for k in range(4):
        angle = base_angle + k * np.pi / 2
        directions.append(np.array([np.cos(angle), np.sin(angle)]))
    
    return directions


def get_best_direction(point, prev_dir, fi, points, cells, u_x, u_y):
    """Get cross direction closest to previous direction."""
    directions = get_field_directions(point, fi, points, cells, u_x, u_y)
    if directions is None:
        return None
    
    best = None
    best_dot = -1
    for d in directions:
        dot = np.dot(prev_dir, d)
        if dot > best_dot:
            best_dot = dot
            best = d
    
    return best


def trace_separatrices(points, cells, u_x, u_y, e1, e2, singularities, vertex_to_faces, boundary):
    """
    Trace separatrices from singularities using RK2 (Heun's method).
    
    FIX: Use e1/e2 from Phase 1 (tangent basis) instead of recomputing with np.cross.
    """
    all_separatrices = []
    
    # Collect singularity locations for collision detection
    sing_locations = np.array([s['location'] for s in singularities])
    sing_faces = [s['face'] for s in singularities]
    
    # Find boundary edges for collision detection
    boundary_edges = []
    edge_to_faces = {}
    for fi, tri in enumerate(cells):
        edges = [
            (int(min(tri[0], tri[1])), int(max(tri[0], tri[1]))),
            (int(min(tri[1], tri[2])), int(max(tri[1], tri[2]))),
            (int(min(tri[2], tri[0])), int(max(tri[2], tri[0]))),
        ]
        for e in edges:
            if e not in edge_to_faces:
                edge_to_faces[e] = []
            edge_to_faces[e].append(fi)
    
    boundary_edges = [e for e, faces in edge_to_faces.items() if len(faces) == 1]
    
    # Compute average edge length for adaptive step size
    avg_edge_length = 0.0
    count = 0
    for fi, tri in enumerate(cells):
        for i in range(3):
            j = (i + 1) % 3
            edge_len = np.linalg.norm(points[tri[j]] - points[tri[i]])
            avg_edge_length += edge_len
            count += 1
    avg_edge_length /= count if count > 0 else 1.0
    h = avg_edge_length * 0.5  # adaptive step size
    
    for sing in singularities:
        sing_loc = sing['location']
        sing_index = sing['index']
        
        # Number of separatrices for this singularity
        # For cross field (N=4): geometric index = representation_index / 4
        # Number of separatrices = abs(representation_index * 4 - 1) is WRONG
        # Correct formula: for geometric index k/N, number of separatrices = |2k - 1|
        # For rep index +1: geom index +1/4, n_sep = |2*1 - 1| = 1 ??? No...
        # From Kowalski code: n_sep = abs(index * 4 - 1)
        # For rep index +1: abs(1*4 - 1) = 3
        # For rep index -1: abs(-1*4 - 1) = 5
        n_sep = abs(int(sing_index * 4 - 1))
        
        # Find separatrix directions by sampling along edges of the triangle containing the singularity
        # (Kowalski algorithm: sample points along triangle edges, check if radial direction aligns with cross field)
        directions_3d = []
        separatrices_angles = []
        
        fi = sing['face']
        tri = cells[fi]
        nodes_of_triangle = points[tri]
        num_nodes_triangle = 3
        
        tolerance = 1e-3
        angle_tolerance = 1e-2
        
        for edge_id in range(num_nodes_triangle):
            next_id = (edge_id + 1) % num_nodes_triangle
            edge_vec = nodes_of_triangle[next_id] - nodes_of_triangle[edge_id]
            edge_length = np.linalg.norm(edge_vec)
            
            if edge_length < 1e-12:
                continue
            
            num_samples = 10000
            t = np.linspace(0, 1, num_samples)
            
            for i in range(num_samples):
                current_node = nodes_of_triangle[edge_id] + t[i] * edge_vec
                radial_vec = current_node - sing_loc
                radial_norm = np.linalg.norm(radial_vec)
                
                if radial_norm < 1e-12:
                    continue
                
                # Project radial_vec into facet's local basis
                radial_2d = np.array([np.dot(radial_vec, e1[fi]), np.dot(radial_vec, e2[fi])])
                radial_2d_norm = np.linalg.norm(radial_2d)
                if radial_2d_norm < 1e-12:
                    continue
                radial_2d = radial_2d / radial_2d_norm
                
                # Get best cross direction at current_node, aligned with radial_vec
                best = get_best_direction(current_node, radial_2d, fi, points, cells, u_x, u_y)
                if best is None:
                    continue
                
                # Convert best to 3D
                best_3d = best[0] * e1[fi] + best[1] * e2[fi]
                best_norm = np.linalg.norm(best_3d)
                if best_norm < 1e-12:
                    continue
                best_3d = best_3d / best_norm
                
                # Check if radial direction aligns with cross field direction
                # Separatrix should point AWAY from singularity
                dot = np.dot(best_3d, radial_vec / radial_norm)
                if dot > 0.95:
                    # Check if we already found this direction
                    is_new = True
                    for existing in directions_3d:
                        if abs(np.dot(best_3d, existing)) > 0.95:
                            is_new = False
                            break
                    
                    if is_new:
                        directions_3d.append(best_3d)
                        separatrices_angles.append(np.arctan2(best[1], best[0]))
        
        # Limit to expected number
        if len(directions_3d) > n_sep:
            # Keep the most distinct directions (sorted by angle)
            directions_3d = directions_3d[:n_sep]
        elif len(directions_3d) < n_sep and len(directions_3d) > 0:
            # If we found fewer than expected, interpolate missing directions
            print(f"Warning: Found {len(directions_3d)} separatrices for singularity {sing['face']} (expected {n_sep})")
            
            # Sort angles
            angles_sorted = sorted(separatrices_angles)
            n_found = len(angles_sorted)
            n_missing = n_sep - n_found
            
            # Find largest angular gaps and place missing directions in the middle
            for _ in range(n_missing):
                if len(angles_sorted) == 1:
                    # Only one direction found: add directions at expected intervals
                    # For n_sep=3 and 1 found: add at +120° and -120°
                    base_angle = angles_sorted[0]
                    for k in range(1, n_sep):
                        new_angle = base_angle + k * 2 * np.pi / n_sep
                        new_angle = new_angle % (2 * np.pi)
                        new_3d = np.cos(new_angle) * e1[fi] + np.sin(new_angle) * e2[fi]
                        new_3d = new_3d / np.linalg.norm(new_3d)
                        directions_3d.append(new_3d)
                        if len(directions_3d) >= n_sep:
                            break
                    break
                else:
                    # Find largest gap between consecutive angles
                    max_gap = 0
                    gap_idx = 0
                    for i in range(len(angles_sorted)):
                        next_i = (i + 1) % len(angles_sorted)
                        gap = (angles_sorted[next_i] - angles_sorted[i]) % (2 * np.pi)
                        if gap > max_gap:
                            max_gap = gap
                            gap_idx = i
                    
                    # Place new direction in middle of largest gap
                    new_angle = angles_sorted[gap_idx] + max_gap / 2
                    new_angle = new_angle % (2 * np.pi)
                    new_3d = np.cos(new_angle) * e1[fi] + np.sin(new_angle) * e2[fi]
                    new_3d = new_3d / np.linalg.norm(new_3d)
                    directions_3d.append(new_3d)
                    angles_sorted.append(new_angle)
                    angles_sorted.sort()
        
        # Trace each direction
        for d_3d in directions_3d:
            streamline = [sing_loc]
            current = sing_loc + h * 0.1 * d_3d  # small initial step
            current_dir = d_3d
            current_fi = find_containing_face(current, points, cells, vertex_to_faces)
            
            # Adaptive step size
            step_h = h
            
            for step in range(2000):
                if current_fi is None:
                    break
                
                # Get best direction at current point
                current_dir_2d = np.array([np.dot(current_dir, e1[current_fi]), np.dot(current_dir, e2[current_fi])])
                if np.linalg.norm(current_dir_2d) < 1e-12:
                    break
                current_dir_2d = current_dir_2d / np.linalg.norm(current_dir_2d)
                
                best = get_best_direction(current, current_dir_2d, current_fi, points, cells, u_x, u_y)
                if best is None:
                    break
                
                # Convert best to 3D
                best_3d = best[0] * e1[current_fi] + best[1] * e2[current_fi]
                best_3d = best_3d / np.linalg.norm(best_3d)
                
                # RK2 predictor-corrector with adaptive step size
                max_attempts = 5
                step_accepted = False
                
                for attempt in range(max_attempts):
                    pred = current + step_h * best_3d
                    pred_fi = find_containing_face(pred, points, cells, vertex_to_faces)
                    
                    if pred_fi is None:
                        # Predictor outside mesh - reduce step size
                        step_h *= 0.5
                        continue
                    
                    # Convert best_3d to pred_fi basis
                    prev_dir_2d_pred = np.array([np.dot(best_3d, e1[pred_fi]), np.dot(best_3d, e2[pred_fi])])
                    if np.linalg.norm(prev_dir_2d_pred) < 1e-12:
                        break
                    prev_dir_2d_pred = prev_dir_2d_pred / np.linalg.norm(prev_dir_2d_pred)
                    
                    best_pred = get_best_direction(pred, prev_dir_2d_pred, pred_fi, points, cells, u_x, u_y)
                    if best_pred is None:
                        break
                    
                    best_pred_3d = best_pred[0] * e1[pred_fi] + best_pred[1] * e2[pred_fi]
                    best_pred_3d = best_pred_3d / np.linalg.norm(best_pred_3d)
                    
                    avg = (best_3d + best_pred_3d) / 2
                    avg = avg / np.linalg.norm(avg)
                    
                    next_point = current + step_h * avg
                    next_fi = find_containing_face(next_point, points, cells, vertex_to_faces)
                    
                    if next_fi is None:
                        # Next point outside mesh - reduce step size
                        step_h *= 0.5
                        continue
                    
                    # Check if near another singularity
                    near_sing = False
                    for si, sloc in enumerate(sing_locations):
                        if sing_faces[si] == sing['face']:
                            continue
                        dist = np.linalg.norm(next_point - sloc)
                        if dist < 0.05:
                            streamline.append(sloc)
                            near_sing = True
                            break
                    
                    if near_sing:
                        break
                    
                    # Check if we are very close to boundary (for robustness)
                    # Only terminate if we are actually on a boundary edge
                    tri_next = cells[next_fi]
                    on_boundary = False
                    for edge_idx in range(3):
                        v0 = tri_next[edge_idx]
                        v1 = tri_next[(edge_idx + 1) % 3]
                        if boundary[v0] and boundary[v1]:
                            edge_start = points[v0]
                            edge_end = points[v1]
                            edge_vec = edge_end - edge_start
                            edge_len = np.linalg.norm(edge_vec)
                            if edge_len > 1e-12:
                                edge_dir = edge_vec / edge_len
                                to_point = next_point - edge_start
                                proj_len = np.dot(to_point, edge_dir)
                                if 0 <= proj_len <= edge_len:
                                    dist = np.linalg.norm(to_point - proj_len * edge_dir)
                                    if dist < 1e-6:  # very strict threshold
                                        on_boundary = True
                                        break
                    
                    if on_boundary:
                        break
                    
                    # Step accepted
                    step_accepted = True
                    current = next_point
                    current_dir = avg
                    current_fi = next_fi
                    streamline.append(current)
                    
                    # Restore step size for next step (but not too aggressively)
                    if attempt > 0:
                        step_h = min(step_h * 1.5, h)
                    
                    break
                
                if not step_accepted:
                    # Step too small or failed, terminate
                    break
            
            if len(streamline) > 1:
                all_separatrices.append(np.array(streamline))
    
    print(f"Traced {len(all_separatrices)} separatrices")
    return all_separatrices


def export_singularities(singularities, out_file):
    """Export singularities as VTK points."""
    if not singularities:
        print("No singularities to export")
        return
    
    points = np.array([s['location'] for s in singularities])
    indices = np.array([s['index'] for s in singularities])
    
    mesh = meshio.Mesh(
        points=points,
        cells=[("vertex", np.arange(len(points)).reshape(-1, 1))],
        point_data={"index": indices},
    )
    mesh.write(out_file)
    print(f"Saved singularities: {out_file}")


def export_separatrices(separatrices, out_file):
    """Export separatrices as VTK lines."""
    if not separatrices:
        print("No separatrices to export")
        return
    
    all_points = []
    all_lines = []
    offset = 0
    
    for sep in separatrices:
        if len(sep) < 2:
            continue
        n = len(sep)
        all_points.extend(sep)
        lines = np.array([[i + offset, i + offset + 1] for i in range(n - 1)])
        all_lines.extend(lines)
        offset += n
    
    if not all_lines:
        print("No valid separatrices to export")
        return
    
    all_points = np.array(all_points)
    all_lines = np.array(all_lines)
    
    mesh = meshio.Mesh(
        points=all_points,
        cells=[("line", all_lines)],
    )
    mesh.write(out_file)
    print(f"Saved separatrices: {out_file}")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("Phase 3: Singularity Detection + Separatrix Tracing")
    print("=" * 60)
    
    # 1. Load Phase 2
    points, cells, frame, frame_angle, boundary, u_x, u_y = load_phase2(IN_VTK)
    
    # 2. Load Phase 1 (normals, e1, e2)
    normals, e1, e2 = load_phase1(PHASE1_VTK)
    
    # 3. Compute vertex basis
    n_v, e1_v, e2_v = compute_vertex_basis(points, cells, normals)
    
    # 4. Build mesh graph
    vertex_to_faces, edge_to_faces = build_mesh_graph(points, cells)
    
    # 5. Detect singularities
    print("\nDetecting singularities...")
    singularities = detect_singularities(points, cells, u_x, u_y, e1_v, e2_v)
    
    # 6. Trace separatrices
    print("\nTracing separatrices...")
    separatrices = trace_separatrices(points, cells, u_x, u_y, e1, e2, singularities, vertex_to_faces, boundary)
    
    # 8. Export
    export_singularities(singularities, OUT_SING)
    export_separatrices(separatrices, OUT_SEP)
    
    print("\n" + "=" * 60)
    print("Phase 3 complete")
    print("=" * 60)


if __name__ == "__main__":
    main()
