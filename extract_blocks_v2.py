#!/usr/bin/env python3
"""
Extract block structure from mesh by building separatrix graph with all vertices.
"""
import sys
sys.path.insert(0, '/root/environments/blocking/lib/python3.12/site-packages')
import gmsh
import numpy as np
import networkx as nx
from collections import defaultdict

MESH_PATH = '/root/repos/block_structured_meshing/T1_9_hub_remeshed.msh'

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
    """Extract the full separatrix graph including all intermediate vertices."""
    # Build edges
    edges = []
    for face_idx in range(len(faces)):
        face = faces[face_idx]
        for i in range(4):
            edges.append((face[i], face[(i+1)%4]))
    
    # Build adjacency
    adj = defaultdict(list)
    for u, v in edges:
        adj[u].append(v)
    
    # Find boundary nodes
    edge_count = {}
    for u, v in edges:
        e = tuple(sorted([u, v]))
        edge_count[e] = edge_count.get(e, 0) + 1
    
    boundary_nodes = set()
    for e, c in edge_count.items():
        if c == 1:
            boundary_nodes.update(e)
    
    # Find singularities
    singularities = []
    for node, neighbors in adj.items():
        val = len(neighbors)
        if node not in boundary_nodes and val != 4:
            singularities.append(node)
    
    print(f"Singularities: {len(singularities)}")
    for node in singularities:
        print(f"  Node {node}: valence {len(adj[node])}")
    
    # Build the separatrix graph
    # Start with boundary edges
    graph_edges = set()
    for e, c in edge_count.items():
        if c == 1:
            graph_edges.add(e)
    
    # Add separatrices (edges from singularities going straight)
    for sing in singularities:
        for neighbor in adj[sing]:
            # Trace path
            path = [sing, neighbor]
            prev = sing
            curr = neighbor
            
            while True:
                if curr in boundary_nodes or curr in singularities:
                    break
                
                neighbors_curr = adj[curr]
                if len(neighbors_curr) != 4:
                    break
                
                # Find opposite neighbor
                prev_idx = neighbors_curr.index(prev)
                next_idx = (prev_idx + 2) % 4
                next_node = neighbors_curr[next_idx]
                
                # Add edge to graph
                edge = tuple(sorted([curr, next_node]))
                graph_edges.add(edge)
                
                path.append(next_node)
                prev = curr
                curr = next_node
                
                if curr == sing:
                    break
            
            # Add the first edge from singularity
            edge = tuple(sorted([sing, neighbor]))
            graph_edges.add(edge)
    
    print(f"Graph edges: {len(graph_edges)}")
    
    # Build graph
    G = nx.Graph()
    G.add_edges_from(graph_edges)
    
    print(f"Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    print(f"Connected components: {nx.number_connected_components(G)}")
    
    return G, boundary_nodes

def find_faces_in_graph(G):
    """Find faces using planar embedding."""
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
    
    print(f"Total faces: {len(faces_raw)}")
    
    # Filter to 4-sided faces
    quad_faces = [face for face in faces_raw if len(face) == 4]
    print(f"Quad faces: {len(quad_faces)}")
    
    return quad_faces

def main():
    print("Loading mesh...")
    vertices, surfaces = load_mesh(MESH_PATH)
    
    for s_tag, faces in surfaces.items():
        print(f"\n=== Surface {s_tag} ===")
        print(f"Faces: {len(faces)}")
        
        G, boundary_nodes = extract_separatrix_graph(vertices, faces)
        
        print("\nFinding faces...")
        quad_faces = find_faces_in_graph(G)
        
        print(f"\nBlock structure: {len(quad_faces)} blocks")
        
        # Export
        if quad_faces:
            export_blocks(vertices, quad_faces, list(G.edges()), 
                         f'/root/repos/block_structured_meshing/T1_9_hub_blocks_s{s_tag}.msh')
    
    print("\nDone.")

def export_blocks(vertices, faces, edges, path):
    """Export block structure as MSH file."""
    all_points = []
    all_elements = []
    point_offset = 1
    element_id = 1
    
    # Export vertices
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
    
    # Export edges
    for u, v in edges:
        all_elements.append({
            'id': element_id,
            'type': 1,
            'nodes': [vertex_map[u], vertex_map[v]]
        })
        element_id += 1
    
    # Export faces
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
