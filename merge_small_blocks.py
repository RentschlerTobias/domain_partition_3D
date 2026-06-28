#!/usr/bin/env python3
"""
Merge small blocks by removing short separatrices.
"""
import sys
sys.path.insert(0, '/root/environments/blocking/lib/python3.12/site-packages')
import gmsh
import numpy as np
import networkx as nx
from collections import defaultdict

MESH_PATH = '/root/repos/block_structured_meshing/T1_9_hub_remeshed.msh'
OUTPUT_PATH = '/root/repos/block_structured_meshing/T1_9_hub_merged_blocks.msh'
TARGET_BLOCKS = 20

def load_mesh(path):
    gmsh.initialize()
    gmsh.open(path)
    
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    vertices = coords.reshape(-1, 3)
    tag_to_idx = np.full(int(node_tags.max()) + 1, -1)
    tag_to_idx[node_tags.astype(int)] = np.arange(len(node_tags))
    
    surfaces = {}
    for _, s_tag in gmsh.model.getEntities(2):
        _, q_tags = gmsh.model.mesh.getElementsByType(3, s_tag)
        if q_tags.size:
            faces = tag_to_idx[q_tags.astype(int)].reshape(-1, 4)
            surfaces[s_tag] = faces
    
    gmsh.finalize()
    return vertices, surfaces

def extract_separatrix_graph(vertices, faces):
    edges = []
    for face_idx in range(len(faces)):
        face = faces[face_idx]
        for i in range(4):
            edges.append((face[i], face[(i+1)%4]))
    
    adj = defaultdict(list)
    for u, v in edges:
        adj[u].append(v)
    
    edge_count = {}
    for u, v in edges:
        e = tuple(sorted([u, v]))
        edge_count[e] = edge_count.get(e, 0) + 1
    
    boundary_nodes = set()
    for e, c in edge_count.items():
        if c == 1:
            boundary_nodes.update(e)
    
    singularities = []
    for node, neighbors in adj.items():
        val = len(neighbors)
        if node not in boundary_nodes and val != 4:
            singularities.append(node)
    
    graph_edges = set()
    for e, c in edge_count.items():
        if c == 1:
            graph_edges.add(e)
    
    for sing in singularities:
        for neighbor in adj[sing]:
            path = [sing, neighbor]
            prev = sing
            curr = neighbor
            
            while True:
                if curr in boundary_nodes or curr in singularities:
                    break
                
                neighbors_curr = adj[curr]
                if len(neighbors_curr) != 4:
                    break
                
                prev_idx = neighbors_curr.index(prev)
                next_idx = (prev_idx + 2) % 4
                next_node = neighbors_curr[next_idx]
                
                edge = tuple(sorted([curr, next_node]))
                graph_edges.add(edge)
                
                path.append(next_node)
                prev = curr
                curr = next_node
                
                if curr == sing:
                    break
            
            edge = tuple(sorted([sing, neighbor]))
            graph_edges.add(edge)
    
    G = nx.Graph()
    G.add_edges_from(graph_edges)
    
    return G, boundary_nodes

def find_faces(G):
    is_planar, embedding = nx.check_planarity(G)
    if not is_planar:
        return []
    
    faces_raw = []
    seen_half_edges = set()
    for u, v in embedding.edges():
        if (u, v) in seen_half_edges:
            continue
        face = embedding.traverse_face(u, v, mark_half_edges=seen_half_edges)
        faces_raw.append(face)
    
    return faces_raw

def count_corners(G, face):
    """Count corners in a face boundary."""
    # Build boundary edges
    boundary_edges = []
    n = len(face)
    for i in range(n):
        u = face[i]
        v = face[(i+1)%n]
        if G.has_edge(u, v):
            boundary_edges.append((u, v))
    
    # Count corners (vertices with degree != 2 in boundary)
    bG = nx.Graph()
    bG.add_edges_from(boundary_edges)
    corners = [n for n in bG.nodes() if bG.degree(n) != 2]
    return len(corners)

def merge_small_faces(G, target_blocks):
    """Merge small faces by removing short edges."""
    faces = find_faces(G)
    
    # Identify removable edges
    # An edge is removable if it is shared by two faces and removing it
    # creates a new face with 4 corners
    
    # Build edge-to-faces mapping
    edge_to_faces = defaultdict(list)
    for face_idx, face in enumerate(faces):
        n = len(face)
        for i in range(n):
            u = face[i]
            v = face[(i+1)%n]
            edge = tuple(sorted([u, v]))
            edge_to_faces[edge].append(face_idx)
    
    # Internal edges (shared by two faces)
    internal_edges = [e for e, fl in edge_to_faces.items() if len(fl) == 2]
    
    # Score edges by length
    edge_scores = {}
    for edge in internal_edges:
        u, v = edge
        length = np.linalg.norm(vertices[u] - vertices[v])
        edge_scores[edge] = length
    
    # Sort by length
    sorted_edges = sorted(edge_scores.items(), key=lambda x: x[1])
    
    # Try removing edges
    while len(faces) > target_blocks and sorted_edges:
        edge, length = sorted_edges.pop(0)
        
        # Check if edge is still in graph
        if not G.has_edge(edge[0], edge[1]):
            continue
        
        # Get the two faces
        face_indices = edge_to_faces[edge]
        if len(face_indices) != 2:
            continue
        
        f1, f2 = face_indices
        
        # Check if both faces are quads
        if len(faces[f1]) != 4 or len(faces[f2]) != 4:
            continue
        
        # Temporarily remove edge
        G_temp = G.copy()
        G_temp.remove_edge(edge[0], edge[1])
        
        # Check if new faces are valid
        try:
            new_faces = find_faces(G_temp)
            # Find the merged face
            for new_face in new_faces:
                if len(new_face) > 4:
                    corners = count_corners(G_temp, new_face)
                    if corners == 4:
                        # Accept merge
                        G = G_temp
                        faces = new_faces
                        
                        # Update edge_to_faces
                        edge_to_faces = defaultdict(list)
                        for face_idx, face in enumerate(faces):
                            n = len(face)
                            for i in range(n):
                                u = face[i]
                                v = face[(i+1)%n]
                                e = tuple(sorted([u, v]))
                                edge_to_faces[e].append(face_idx)
                        
                        # Recompute internal edges
                        internal_edges = [e for e, fl in edge_to_faces.items() if len(fl) == 2]
                        edge_scores = {}
                        for e in internal_edges:
                            u, v = e
                            length = np.linalg.norm(vertices[u] - vertices[v])
                            edge_scores[e] = length
                        sorted_edges = sorted(edge_scores.items(), key=lambda x: x[1])
                        break
        except:
            pass
    
    return G, faces

def main():
    global vertices
    
    print("Loading mesh...")
    vertices, surfaces = load_mesh(MESH_PATH)
    
    for s_tag, faces in surfaces.items():
        print(f"\n=== Surface {s_tag} ===")
        print(f"Faces: {len(faces)}")
        
        G, boundary_nodes = extract_separatrix_graph(vertices, faces)
        
        print(f"Initial graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
        
        initial_faces = find_faces(G)
        print(f"Initial faces: {len(initial_faces)}")
        
        # Merge small faces
        G_merged, merged_faces = merge_small_faces(G, TARGET_BLOCKS)
        
        print(f"Merged faces: {len(merged_faces)}")
        
        # Count quad faces
        quad_faces = [f for f in merged_faces if len(f) == 4]
        print(f"Quad faces: {len(quad_faces)}")
        
        # Export
        export_blocks(vertices, merged_faces, list(G_merged.edges()), 
                     f'/root/repos/block_structured_meshing/T1_9_hub_merged_s{s_tag}.msh')
    
    print("\nDone.")

def export_blocks(vertices, faces, edges, path):
    all_points = []
    all_elements = []
    point_offset = 1
    element_id = 1
    
    used_vertices = set()
    for face in faces:
        used_vertices.update(face)
    for u, v in edges:
        used_vertices.add(u)
        used_vertices.add(v)
    
    vertex_map = {}
    for i, v_idx in enumerate(sorted(used_vertices)):
        vertex_map[v_idx] = point_offset + i
        v = vertices[v_idx]
        all_points.append((point_offset + i, v[0], v[1], v[2]))
    
    for u, v in edges:
        all_elements.append({
            'id': element_id,
            'type': 1,
            'nodes': [vertex_map[u], vertex_map[v]]
        })
        element_id += 1
    
    for face in faces:
        all_elements.append({
            'id': element_id,
            'type': 3,
            'nodes': [vertex_map[v] for v in face]
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
    
    print(f"\nExported blocks to {path}")

if __name__ == '__main__':
    main()
