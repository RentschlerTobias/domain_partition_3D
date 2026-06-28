#!/usr/bin/env python3
"""
Build block structure from streamlines and merge to coarse blocks.
"""
import sys
sys.path.insert(0, '/root/environments/blocking/lib/python3.12/site-packages')
import gmsh
import numpy as np
import networkx as nx
from collections import defaultdict, Counter

STREAMLINE_PATH = '/root/repos/block_structured_meshing/streamlines_s2.msh'
MESH_PATH = '/root/repos/block_structured_meshing/T1_9_hub_remeshed.msh'

def load_streamlines(path):
    """Load streamlines from MSH file."""
    streamlines = []
    with open(path, 'r') as f:
        lines = f.readlines()
        in_nodes = False
        in_elements = False
        nodes = {}
        elements = []
        
        for line in lines:
            if '$Nodes' in line:
                in_nodes = True
                continue
            if '$EndNodes' in line:
                in_nodes = False
                continue
            if '$Elements' in line:
                in_elements = True
                continue
            if '$EndElements' in line:
                in_elements = False
                continue
            
            if in_nodes:
                parts = line.strip().split()
                if len(parts) == 4:
                    tag = int(parts[0])
                    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                    nodes[tag] = (x, y, z)
            
            if in_elements:
                parts = line.strip().split()
                if len(parts) >= 6:
                    n1 = int(parts[-2])
                    n2 = int(parts[-1])
                    elements.append((n1, n2))
    
    # Build streamlines from connected elements
    G = nx.Graph()
    G.add_edges_from(elements)
    
    # Find connected components (each is a streamline)
    for comp in nx.connected_components(G):
        if len(comp) < 2:
            continue
        subG = G.subgraph(comp)
        # Find endpoints
        endpoints = [n for n in subG.nodes() if subG.degree(n) == 1]
        if len(endpoints) == 2:
            # Path
            path = nx.shortest_path(subG, endpoints[0], endpoints[1])
            streamlines.append([nodes[n] for n in path])
        elif len(endpoints) == 0:
            # Loop
            cycle = nx.find_cycle(subG)
            path = [cycle[0][0]] + [v for _, v in cycle]
            streamlines.append([nodes[n] for n in path])
    
    return streamlines

def build_graph_from_streamlines(streamlines):
    """Build graph from streamlines with rounded vertices."""
    # Round to create unique vertices
    point_to_vertex = {}
    vertices_list = []
    vertex_counter = 0
    
    edges = []
    edge_to_streamline = []
    
    for sl_idx, streamline in enumerate(streamlines):
        start_point = np.array(streamline[0])
        end_point = np.array(streamline[-1])
        
        start_key = tuple(np.round(start_point, decimals=6))
        end_key = tuple(np.round(end_point, decimals=6))
        
        if start_key not in point_to_vertex:
            point_to_vertex[start_key] = vertex_counter
            vertices_list.append(start_point)
            vertex_counter += 1
        
        if end_key not in point_to_vertex:
            point_to_vertex[end_key] = vertex_counter
            vertices_list.append(end_point)
            vertex_counter += 1
        
        start_idx = point_to_vertex[start_key]
        end_idx = point_to_vertex[end_key]
        
        edges.append((start_idx, end_idx))
        edge_to_streamline.append(sl_idx)
    
    vertices = np.array(vertices_list)
    return vertices, edges, edge_to_streamline, streamlines

def detect_faces_in_graph(vertices, edges):
    """Detect faces in the streamline graph using planar embedding."""
    G = nx.Graph()
    G.add_edges_from(edges)
    
    is_planar, embedding = nx.check_planarity(G)
    if not is_planar:
        print("WARNING: Graph not planar")
        return []
    
    faces_raw = []
    seen_half_edges = set()
    for u, v in embedding.edges():
        if (u, v) in seen_half_edges:
            continue
        face = embedding.traverse_face(u, v, mark_half_edges=seen_half_edges)
        faces_raw.append(face)
    
    # Filter to quad faces (4 vertices)
    quad_faces = [face for face in faces_raw if len(face) == 4]
    print(f"Total faces: {len(faces_raw)}, Quad faces: {len(quad_faces)}")
    
    return quad_faces

def validate_block_structure(vertices, edges, faces):
    """Basic validation."""
    G = nx.Graph()
    G.add_edges_from(edges)
    
    # Check Euler characteristic
    V = len(vertices)
    E = len(edges)
    F = len(faces)
    euler = V - E + F
    
    # Find boundary loops
    boundary_edges = []
    edge_count = Counter()
    for u, v in edges:
        edge_count[tuple(sorted([u, v]))] += 1
    
    for e, c in edge_count.items():
        if c == 1:
            boundary_edges.append(e)
    
    boundary_loops = 0
    if boundary_edges:
        bG = nx.Graph()
        bG.add_edges_from(boundary_edges)
        for comp in nx.connected_components(bG):
            boundary_loops += 1
    
    expected_euler = 2 - boundary_loops
    print(f"V={V}, E={E}, F={F}, Euler={euler}, Expected={expected_euler}, Boundary loops={boundary_loops}")
    
    if abs(euler - expected_euler) > 0.1:
        print("WARNING: Euler characteristic mismatch")
        return False
    
    return True

def main():
    print("Loading streamlines...")
    streamlines = load_streamlines(STREAMLINE_PATH)
    print(f"Loaded {len(streamlines)} streamlines")
    
    print("\nBuilding graph...")
    vertices, edges, edge_to_streamline, streamlines = build_graph_from_streamlines(streamlines)
    print(f"Vertices: {len(vertices)}, Edges: {len(edges)}")
    
    print("\nDetecting faces...")
    faces = detect_faces_in_graph(vertices, edges)
    
    print("\nValidating...")
    is_valid = validate_block_structure(vertices, edges, faces)
    
    if is_valid:
        print(f"\nBlock structure: {len(faces)} blocks")
        print(f"Target: ~20 blocks")
        
        if len(faces) > 50:
            print("\nWARNING: Too many blocks. Need merging.")
        else:
            print("\nBlock count is acceptable.")
    
    # Export coarse blocks
    export_coarse_blocks(vertices, faces, edges, edge_to_streamline, streamlines, 
                        '/root/repos/block_structured_meshing/T1_9_hub_coarse_blocks.msh')

def export_coarse_blocks(vertices, faces, edges, edge_to_streamline, streamlines, path):
    """Export coarse block structure as MSH file."""
    all_points = []
    all_elements = []
    point_offset = 1
    element_id = 1
    
    # Export vertices as points
    for i, v in enumerate(vertices):
        all_points.append((point_offset + i, v[0], v[1], v[2]))
    
    # Export edges as lines
    for i, (u, v) in enumerate(edges):
        all_elements.append({
            'id': element_id,
            'type': 1,
            'nodes': [point_offset + u, point_offset + v]
        })
        element_id += 1
    
    # Export faces as quads
    for face in faces:
        all_elements.append({
            'id': element_id,
            'type': 3,
            'nodes': [point_offset + v for v in face]
        })
        element_id += 1
    
    with open(path, 'w') as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        f.write("$Nodes\n")
        f.write(f"{len(all_points)}\n")
        for nid, x, y, z in all_points:
            f.write(f"{nid} {x:.10e} {y:.10e} {z:.10e}\n")
        f.write("$EndNodes\n")
        f.write("$Elements\n")
        f.write(f"{len(all_elements)}\n")
        for e in all_elements:
            nodes_str = " ".join(str(n) for n in e['nodes'])
            f.write(f"{e['id']} {e['type']} 2 1 1 {nodes_str}\n")
        f.write("$EndElements\n")
    
    print(f"\nExported coarse blocks to {path}")

if __name__ == '__main__':
    main()
