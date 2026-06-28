#!/usr/bin/env python3
"""
Complete pipeline for T1_9 hub surface domain partition.
Steps:
1. Extract hub surface from STL
2. Remesh with Gmsh quasi-structured quad algorithm
3. Extract separatrices and build block structure
4. Export coarse blocks
"""
import sys
sys.path.insert(0, '/root/environments/blocking/lib/python3.12/site-packages')
import gmsh
import numpy as np
import networkx as nx
from collections import defaultdict

# Paths
STL_PATH = '/root/repos/block_structured_meshing/T1_9/T1_9_ru_gridGmsh.stl'
HUB_STL_PATH = '/root/repos/block_structured_meshing/T1_9_hub.stl'
REMESHED_PATH = '/root/repos/block_structured_meshing/T1_9_hub_remeshed.msh'
OUTPUT_PATH = '/root/repos/block_structured_meshing/T1_9_hub_blocks.msh'

# Parameters
ELEMENT_SIZE = 1.0
ANGLE = 60

def extract_hub_stl():
    """Extract hub surface (z=0) from STL."""
    print("Extracting hub surface from STL...")
    
    facets = []
    with open(STL_PATH, 'r') as f:
        lines = f.readlines()
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if line.startswith('facet normal'):
                parts = line.split()
                normal = [float(parts[2]), float(parts[3]), float(parts[4])]
                
                i += 1
                if 'outer loop' not in lines[i].strip():
                    i += 1
                    continue
                
                vertices = []
                for j in range(3):
                    i += 1
                    vline = lines[i].strip()
                    parts = vline.split()
                    if parts[0] == 'vertex':
                        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                        vertices.append((x, y, z))
                
                all_z0 = all(abs(v[2]) < 1e-6 for v in vertices)
                if all_z0:
                    facets.append((normal, vertices))
            i += 1
    
    with open(HUB_STL_PATH, 'w') as f:
        f.write("solid hub\n")
        for normal, vertices in facets:
            f.write(f"  facet normal {normal[0]} {normal[1]} {normal[2]}\n")
            f.write("    outer loop\n")
            for v in vertices:
                f.write(f"      vertex {v[0]} {v[1]} {v[2]}\n")
            f.write("    endloop\n")
            f.write("  endfacet\n")
        f.write("endsolid hub\n")
    
    print(f"Extracted {len(facets)} hub facets to {HUB_STL_PATH}")

def remesh_hub():
    """Remesh hub surface with quasi-structured quad algorithm."""
    print("\nRemeshing hub surface...")
    
    gmsh.initialize()
    gmsh.merge(HUB_STL_PATH)
    
    gmsh.model.mesh.classifySurfaces(ANGLE * np.pi / 180, True, False, 180 * np.pi / 180)
    gmsh.model.mesh.createGeometry()
    
    gmsh.option.setNumber("Mesh.Algorithm", 11)
    gmsh.option.setNumber("Mesh.MeshSizeMax", ELEMENT_SIZE)
    gmsh.option.setNumber("Mesh.MeshSizeMin", ELEMENT_SIZE * 0.5)
    
    gmsh.model.mesh.generate(2)
    gmsh.write(REMESHED_PATH)
    
    # Count elements
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    total_nodes = len(node_tags)
    
    total_quads = 0
    for _, s_tag in gmsh.model.getEntities(2):
        elem_types, elem_tags, _ = gmsh.model.mesh.getElements(2, s_tag)
        for et, etags in zip(elem_types, elem_tags):
            if et == 3:
                total_quads += len(etags)
    
    gmsh.finalize()
    
    print(f"Remeshed mesh: {total_nodes} nodes, {total_quads} quads")
    print(f"Saved to {REMESHED_PATH}")
    
    return total_quads

def extract_blocks():
    """Extract block structure from remeshed mesh."""
    print("\nExtracting block structure...")
    
    gmsh.initialize()
    gmsh.open(REMESHED_PATH)
    
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
    
    for s_tag, faces in surfaces.items():
        print(f"\nProcessing surface {s_tag}...")
        
        # Build mesh graph
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
        
        print(f"  Singularities: {len(singularities)}")
        
        # Build separatrix graph
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
        
        print(f"  Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
        
        # Find faces using planar embedding
        is_planar, embedding = nx.check_planarity(G)
        if not is_planar:
            print("  ERROR: Graph not planar")
            continue
        
        faces_raw = []
        seen_half_edges = set()
        for u, v in embedding.edges():
            if (u, v) in seen_half_edges:
                continue
            face = embedding.traverse_face(u, v, mark_half_edges=seen_half_edges)
            faces_raw.append(face)
        
        quad_faces = [face for face in faces_raw if len(face) == 4]
        print(f"  Total faces: {len(faces_raw)}, Quad faces: {len(quad_faces)}")
        
        # Export
        export_blocks(vertices, quad_faces, list(G.edges()), OUTPUT_PATH)
        print(f"  Exported blocks to {OUTPUT_PATH}")
    
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

def main():
    extract_hub_stl()
    remesh_hub()
    extract_blocks()

if __name__ == '__main__':
    main()
