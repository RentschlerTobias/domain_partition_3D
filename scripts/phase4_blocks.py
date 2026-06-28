#!/usr/bin/env python3
"""
Phase 4: Block Partitioning (4-sided faces)

Input:  output/T1_9/hub/hub_phase3_separatrices.vtk
        output/T1_9/hub/hub_phase3_singularities.vtk
        output/T1_9/hub/hub_phase1_boundary.vtk
Output: output/T1_9/hub/hub_phase4_blocks.vtk

Algorithm:
1. Load separatrices + singularities + boundary
2. Build simplified graph:
   - Nodes = singularities + boundary vertices where separatrices hit
   - Boundary edges are merged into segments between hit points
   - Separatrices are edges
3. Find 4-sided faces (cycles with 4 edges)
4. Export as mesh blocks
"""

import meshio
import numpy as np
from pathlib import Path
from collections import defaultdict
import networkx as nx

# Paths
IN_SEP = Path("/root/repos/block_structured_meshing/output/T1_9/hub/hub_phase3_separatrices.vtk")
IN_SING = Path("/root/repos/block_structured_meshing/output/T1_9/hub/hub_phase3_singularities.vtk")
IN_BOUNDARY = Path("/root/repos/block_structured_meshing/output/T1_9/hub/hub_phase1_boundary.vtk")
OUT_DIR = Path("/root/repos/block_structured_meshing/output/T1_9/hub")
OUT_BLOCKS = OUT_DIR / "hub_phase4_blocks.vtk"


def load_separatrices(vtk_file):
    """Load separatrices as list of polylines."""
    mesh = meshio.read(vtk_file)
    points = mesh.points
    lines = mesh.cells[0].data
    
    # Build polylines from connected lines
    polylines = []
    used = set()
    
    for start_line in range(len(lines)):
        if start_line in used:
            continue
        
        poly = [points[lines[start_line][0]], points[lines[start_line][1]]]
        used.add(start_line)
        
        current_end = lines[start_line][1]
        while True:
            found = False
            for li in range(len(lines)):
                if li in used:
                    continue
                if lines[li][0] == current_end:
                    poly.append(points[lines[li][1]])
                    current_end = lines[li][1]
                    used.add(li)
                    found = True
                    break
                elif lines[li][1] == current_end:
                    poly.append(points[lines[li][0]])
                    current_end = lines[li][0]
                    used.add(li)
                    found = True
                    break
            if not found:
                break
        
        current_start = lines[start_line][0]
        while True:
            found = False
            for li in range(len(lines)):
                if li in used:
                    continue
                if lines[li][0] == current_start:
                    poly.insert(0, points[lines[li][1]])
                    current_start = lines[li][1]
                    used.add(li)
                    found = True
                    break
                elif lines[li][1] == current_start:
                    poly.insert(0, points[lines[li][0]])
                    current_start = lines[li][0]
                    used.add(li)
                    found = True
                    break
            if not found:
                break
        
        if len(poly) > 1:
            polylines.append(np.array(poly))
    
    print(f"Loaded {len(polylines)} separatrices")
    for i, poly in enumerate(polylines):
        print(f"  Sep {i}: {len(poly)} points, start={poly[0]}, end={poly[-1]}")
    
    return polylines


def load_singularities(vtk_file):
    """Load singularity points."""
    mesh = meshio.read(vtk_file)
    points = mesh.points
    indices = mesh.point_data["index"]
    
    print(f"Loaded {len(points)} singularities")
    return points, indices


def load_boundary(vtk_file):
    """Load boundary edges and build boundary ring."""
    mesh = meshio.read(vtk_file)
    points = mesh.points
    lines = mesh.cells[0].data
    
    # Build boundary ring
    # Find unique boundary vertices
    boundary_vertices = np.unique(lines.flatten())
    
    # Build adjacency
    adj = defaultdict(list)
    for line in lines:
        a, b = line
        adj[a].append(b)
        adj[b].append(a)
    
    # Find boundary ring (ordered)
    ring = []
    start = boundary_vertices[0]
    current = start
    prev = None
    
    while True:
        ring.append(current)
        neighbors = [n for n in adj[current] if n != prev]
        if not neighbors:
            break
        prev = current
        current = neighbors[0]
        if current == start:
            break
    
    ring_coords = points[ring]
    
    print(f"Loaded boundary: {len(ring)} vertices, {len(lines)} edges")
    print(f"Ring length: {len(ring)}")
    
    return points, lines, ring, ring_coords


def snap_to_boundary(endpoints, boundary_ring, boundary_coords):
    """Snap endpoints to nearest boundary vertex."""
    snapped = []
    for pt in endpoints:
        dists = np.linalg.norm(boundary_coords - pt, axis=1)
        idx = np.argmin(dists)
        snapped.append((idx, boundary_ring[idx], dists[idx]))
    
    return snapped


def build_graph(singularities, separatrices, boundary_ring, boundary_coords):
    """Build simplified graph with merged boundary segments."""
    
    # Get endpoints of separatrices
    endpoints = []
    for sep in separatrices:
        if len(sep) < 2:
            continue
        endpoints.append(sep[0])  # start
        endpoints.append(sep[-1])  # end
    
    # Snap to boundary (find closest boundary vertex for each endpoint)
    snapped = snap_to_boundary(endpoints, boundary_ring, boundary_coords)
    
    # Identify which separatrices connect to which boundary vertices
    hit_vertices = set()
    sep_connections = []
    
    for i, sep in enumerate(separatrices):
        if len(sep) < 2:
            continue
        
        start_idx = 2 * i
        end_idx = 2 * i + 1
        
        start_snap = snapped[start_idx]
        end_snap = snapped[end_idx]
        
        # Check which endpoints are on boundary
        start_on_boundary = start_snap[2] < 1.0  # tolerance
        end_on_boundary = end_snap[2] < 1.0
        
        if start_on_boundary:
            hit_vertices.add(start_snap[1])
        if end_on_boundary:
            hit_vertices.add(end_snap[1])
        
        sep_connections.append({
            'start_snap': start_snap,
            'end_snap': end_snap,
            'start_on_boundary': start_on_boundary,
            'end_on_boundary': end_on_boundary,
        })
    
    print(f"Boundary hit vertices: {len(hit_vertices)}")
    
    # Build graph
    G = nx.Graph()
    
    # Add singularity nodes
    for i, sing in enumerate(singularities):
        G.add_node(f"sing_{i}", pos=sing, type="singularity")
    
    # Add boundary hit nodes
    hit_list = sorted(list(hit_vertices))
    for i, hv in enumerate(hit_list):
        G.add_node(f"bound_{i}", pos=boundary_coords[boundary_ring.index(hv)] if hv in boundary_ring else None, type="boundary")
    
    # Map boundary vertex to node name
    bv_to_node = {}
    for i, hv in enumerate(hit_list):
        bv_to_node[hv] = f"bound_{i}"
    
    # Add boundary segments
    n_ring = len(boundary_ring)
    for i in range(len(hit_list)):
        v1 = hit_list[i]
        v2 = hit_list[(i + 1) % len(hit_list)]
        
        # Find arc length along boundary
        idx1 = boundary_ring.index(v1)
        idx2 = boundary_ring.index(v2)
        
        # Compute arc length
        if idx2 > idx1:
            arc_len = idx2 - idx1
        else:
            arc_len = n_ring - idx1 + idx2
        
        # Add edge
        node1 = bv_to_node[v1]
        node2 = bv_to_node[v2]
        if node1 != node2:
            G.add_edge(node1, node2, type="boundary", weight=arc_len)
    
    # Add separatrix edges
    for i, conn in enumerate(sep_connections):
        sep = separatrices[i]
        
        # Start node
        if conn['start_on_boundary']:
            start_node = bv_to_node[conn['start_snap'][1]]
        else:
            # Add as singularity
            start_node = None
            for j, sing in enumerate(singularities):
                if np.linalg.norm(sep[0] - sing) < 1e-3:
                    start_node = f"sing_{j}"
                    break
        
        # End node
        if conn['end_on_boundary']:
            end_node = bv_to_node[conn['end_snap'][1]]
        else:
            end_node = None
            for j, sing in enumerate(singularities):
                if np.linalg.norm(sep[-1] - sing) < 1e-3:
                    end_node = f"sing_{j}"
                    break
        
        if start_node and end_node and start_node != end_node:
            G.add_edge(start_node, end_node, type="separatrix", weight=len(sep))
    
    print(f"Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    
    return G


def find_4sided_faces(G):
    """Find 4-sided faces (cycles with 4 edges)."""
    faces = []
    
    # Enumerate all simple cycles of length 4
    for cycle in nx.simple_cycles(G.to_directed()):
        if len(cycle) == 4:
            # Check if it's a valid cycle
            is_valid = True
            for i in range(4):
                a = cycle[i]
                b = cycle[(i+1) % 4]
                if not G.has_edge(a, b):
                    is_valid = False
                    break
            
            if is_valid:
                sorted_cycle = tuple(sorted(cycle))
                if sorted_cycle not in [tuple(sorted(f)) for f in faces]:
                    faces.append(tuple(cycle))
    
    print(f"Found {len(faces)} 4-sided faces")
    return faces


def export_blocks(G, faces, out_file):
    """Export blocks as VTK mesh."""
    if not faces:
        print("No faces to export")
        return
    
    vertex_map = {}
    all_points = []
    
    for face in faces:
        for node_id in face:
            if node_id not in vertex_map:
                vertex_map[node_id] = len(all_points)
                all_points.append(G.nodes[node_id]["pos"])
    
    all_points = np.array(all_points)
    
    quads = []
    for face in faces:
        quad = [vertex_map[node_id] for node_id in face]
        quads.append(quad)
    
    quads = np.array(quads)
    
    mesh = meshio.Mesh(
        points=all_points,
        cells=[("quad", quads)],
    )
    mesh.write(out_file)
    print(f"Saved blocks: {out_file}")
    print(f"  Points: {len(all_points)}, Quads: {len(quads)}")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("Phase 4: Block Partitioning")
    print("=" * 60)
    
    # 1. Load
    separatrices = load_separatrices(IN_SEP)
    sing_points, sing_indices = load_singularities(IN_SING)
    bound_points, bound_lines, bound_ring, bound_ring_coords = load_boundary(IN_BOUNDARY)
    
    # 2. Build graph
    G = build_graph(sing_points, separatrices, bound_ring, bound_ring_coords)
    
    # 3. Find faces
    faces = find_4sided_faces(G)
    
    # 4. Export
    export_blocks(G, faces, OUT_BLOCKS)
    
    print("\n" + "=" * 60)
    print("Phase 4 complete")
    print("=" * 60)


if __name__ == "__main__":
    main()
