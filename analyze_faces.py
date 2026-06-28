#!/usr/bin/env python3
"""
Analyze the faces found by planar embedding.
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

def main():
    vertices, surfaces = load_mesh(MESH_PATH)
    
    for s_tag, faces in surfaces.items():
        G, boundary_nodes = extract_separatrix_graph(vertices, faces)
        
        is_planar, embedding = nx.check_planarity(G)
        if not is_planar:
            print("Not planar")
            return
        
        faces_raw = []
        seen_half_edges = set()
        for u, v in embedding.edges():
            if (u, v) in seen_half_edges:
                continue
            face = embedding.traverse_face(u, v, mark_half_edges=seen_half_edges)
            faces_raw.append(face)
        
        print(f"Total faces: {len(faces_raw)}")
        
        # Analyze each face
        for i, face in enumerate(faces_raw):
            # Check if outer face (contains boundary nodes)
            has_boundary = any(n in boundary_nodes for n in face)
            print(f"  Face {i}: {len(face)} vertices, boundary={has_boundary}")
            
            if len(face) <= 6:
                # Print vertices
                coords = [vertices[n] for n in face]
                print(f"    Coordinates: {[(round(c[0], 3), round(c[1], 3)) for c in coords]}")

if __name__ == '__main__':
    main()
