#!/usr/bin/env python3
"""
Robust 3D surface domain partition pipeline for T1_9 hub.
Steps:
1. Load remeshed quad mesh
2. Extract faces, vertices, edges
3. Build adjacency with geometric edge sorting
4. Detect quad faces using planar embedding
5. Merge blocks greedily
6. Validate
7. Export coarse block structure
"""
import sys
sys.path.insert(0, '/root/environments/blocking/lib/python3.12/site-packages')
import gmsh
import numpy as np
import networkx as nx
from collections import defaultdict
import meshio

MESH_PATH = '/root/repos/block_structured_meshing/T1_9_hub_remeshed.msh'
OUTPUT_PATH = '/root/repos/block_structured_meshing/T1_9_hub_blocks.msh'

def load_mesh(path):
    """Load mesh and extract vertices, quad faces per surface."""
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

def build_adjacency_and_sort(vertices, faces):
    """Build adjacency and sort edges geometrically around each node."""
    # Build edges
    edges = []
    for face_idx in range(len(faces)):
        face = faces[face_idx]
        for i in range(4):
            edges.append((face[i], face[(i+1)%4]))
    
    # Build adjacency list
    adj = defaultdict(list)
    for u, v in edges:
        adj[u].append(v)
    
    # Sort edges geometrically around each node
    # Compute local normal (average of face normals)
    node_normals = np.zeros_like(vertices)
    node_face_count = np.zeros(len(vertices))
    
    for face_idx in range(len(faces)):
        face = faces[face_idx]
        p0, p1, p2, p3 = vertices[face]
        # Two triangle normals
        n1 = np.cross(p1 - p0, p2 - p0)
        n2 = np.cross(p2 - p0, p3 - p0)
        n = n1 + n2
        n_norm = np.linalg.norm(n)
        if n_norm > 1e-10:
            n /= n_norm
        for v in face:
            node_normals[v] += n
            node_face_count[v] += 1
    
    for i in range(len(vertices)):
        if node_face_count[i] > 0:
            node_normals[i] /= node_face_count[i]
            nn = np.linalg.norm(node_normals[i])
            if nn > 1e-10:
                node_normals[i] /= nn
    
    # Sort neighbors by angle in tangent plane
    sorted_adj = {}
    for node, neighbors in adj.items():
        if len(neighbors) <= 2:
            sorted_adj[node] = neighbors
            continue
        
        p = vertices[node]
        n = node_normals[node]
        
        # Create basis in tangent plane
        if abs(n[2]) > 0.9:
            t1 = np.array([1.0, 0.0, 0.0])
        else:
            t1 = np.array([0.0, 0.0, 1.0])
        t1 = t1 - np.dot(t1, n) * n
        t1_norm = np.linalg.norm(t1)
        if t1_norm > 1e-10:
            t1 /= t1_norm
        else:
            t1 = np.array([1.0, 0.0, 0.0])
        t2 = np.cross(n, t1)
        
        angles = []
        for nb in neighbors:
            v = vertices[nb] - p
            v = v - np.dot(v, n) * n  # Project to tangent plane
            vn = np.linalg.norm(v)
            if vn > 1e-10:
                v /= vn
            a = np.arctan2(np.dot(v, t2), np.dot(v, t1))
            angles.append(a)
        
        sorted_neighbors = [n for _, n in sorted(zip(angles, neighbors))]
        sorted_adj[node] = sorted_neighbors
    
    return sorted_adj, edges

def detect_quad_faces(vertices, faces):
    """Detect quad faces using planar embedding."""
    # Build graph with sorted edges
    sorted_adj, edge_list = build_adjacency_and_sort(vertices, faces)
    
    # Create networkx graph
    G = nx.Graph()
    G.add_edges_from(edge_list)
    
    # Check planarity
    is_planar, embedding = nx.check_planarity(G)
    if not is_planar:
        print("WARNING: Graph not planar, using 4-cycle detection fallback")
        return detect_quad_faces_fallback(vertices, faces)
    
    # Traverse faces
    faces_raw = []
    seen_half_edges = set()
    for u, v in embedding.edges():
        if (u, v) in seen_half_edges:
            continue
        face = embedding.traverse_face(u, v, mark_half_edges=seen_half_edges)
        faces_raw.append(face)
    
    # Keep only quad faces
    quad_faces = [face for face in faces_raw if len(face) == 4]
    print(f"Total faces: {len(faces_raw)}, Quad faces: {len(quad_faces)}")
    
    if len(quad_faces) == 0:
        return np.empty((4, 0), dtype=np.int64)
    
    return np.array(quad_faces).T

def detect_quad_faces_fallback(vertices, faces):
    """Fallback: detect 4-cycles."""
    # Build edges from faces
    edges = []
    for face_idx in range(len(faces)):
        face = faces[face_idx]
        for i in range(4):
            edges.append((face[i], face[(i+1)%4]))
    
    G = nx.Graph()
    G.add_edges_from(edges)
    
    quad_faces = []
    seen = set()
    for cycle in nx.simple_cycles(G, length_bound=4):
        if len(cycle) == 4:
            face = tuple(sorted(cycle))
            if face not in seen:
                seen.add(face)
                quad_faces.append(list(cycle))
    
    if len(quad_faces) == 0:
        return np.empty((4, 0), dtype=np.int64)
    
    return np.array(quad_faces).T

def trace_streamlines(vertices, faces):
    """Trace streamlines from irregular vertices."""
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
    
    # Trace streamlines
    streamlines = []
    visited_edges = set()
    
    for sing in singularities:
        for neighbor in adj[sing]:
            # Trace from sing to neighbor
            path = [vertices[sing], vertices[neighbor]]
            prev = sing
            curr = neighbor
            
            while True:
                # Check if boundary or singularity
                if curr in boundary_nodes or curr in singularities:
                    break
                
                # Find next node (go straight)
                neighbors_curr = adj[curr]
                if len(neighbors_curr) != 4:
                    break
                
                # Find opposite neighbor
                prev_idx = neighbors_curr.index(prev)
                next_idx = (prev_idx + 2) % 4
                next_node = neighbors_curr[next_idx]
                
                path.append(vertices[next_node])
                prev = curr
                curr = next_node
                
                if curr == sing:
                    break
            
            # Check if edge already visited
            edge_key = tuple(sorted([sing, curr]))
            if edge_key in visited_edges:
                continue
            visited_edges.add(edge_key)
            
            if len(path) > 1:
                streamlines.append(np.array(path))
    
    # Also trace boundaries
    boundary_edges_list = [e for e, c in edge_count.items() if c == 1]
    if boundary_edges_list:
        bG = nx.Graph()
        bG.add_edges_from(boundary_edges_list)
        for comp in nx.connected_components(bG):
            subG = bG.subgraph(comp)
            try:
                cycle = nx.find_cycle(subG)
                path = [cycle[0][0]] + [v for _, v in cycle]
                path_coords = [vertices[n] for n in path]
                streamlines.append(np.array(path_coords))
            except nx.NetworkXNoCycle:
                # Open chain
                endpoints = [n for n in subG.nodes() if subG.degree(n) == 1]
                if len(endpoints) == 2:
                    path = nx.shortest_path(subG, endpoints[0], endpoints[1])
                    path_coords = [vertices[n] for n in path]
                    streamlines.append(np.array(path_coords))
    
    return streamlines

def main():
    print("Loading mesh...")
    vertices, surfaces = load_mesh(MESH_PATH)
    
    for s_tag, faces in surfaces.items():
        print(f"\n=== Surface {s_tag} ===")
        print(f"Vertices: {len(vertices)}, Faces: {len(faces)}")
        
        # Build edges
        edges = []
        for face_idx in range(len(faces)):
            face = faces[face_idx]
            for i in range(4):
                edges.append((face[i], face[(i+1)%4]))
        
        print(f"Edges: {len(edges)}")
        
        # Detect quad faces
        quad_faces = detect_quad_faces(vertices, faces)
        print(f"Detected quad faces: {quad_faces.shape[1] if quad_faces.size > 0 else 0}")
        
        # Trace streamlines
        streamlines = trace_streamlines(vertices, faces)
        print(f"Streamlines: {len(streamlines)}")
        
        # Export streamlines
        if streamlines:
            export_streamlines(streamlines, f'/root/repos/block_structured_meshing/streamlines_s{s_tag}.msh')
    
    print("\nDone.")

def export_streamlines(streamlines, path):
    """Export streamlines to MSH file."""
    all_points = []
    all_elements = []
    point_offset = 1
    element_id = 1
    
    for points in streamlines:
        start_idx = point_offset
        for pt in points:
            all_points.append((point_offset, pt[0], pt[1], pt[2]))
            point_offset += 1
        
        for i in range(len(points) - 1):
            all_elements.append({
                'id': element_id,
                'nodes': [start_idx + i, start_idx + i + 1]
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
            f.write(f"{e['id']} 1 2 1 1 {e['nodes'][0]} {e['nodes'][1]}\n")
        f.write("$EndElements\n")
    
    print(f"Exported streamlines to {path}")

if __name__ == '__main__':
    main()
