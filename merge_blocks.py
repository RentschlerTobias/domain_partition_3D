#!/usr/bin/env python3
"""
Merge quad mesh faces into coarse blocks.
"""
import sys
sys.path.insert(0, '/root/environments/blocking/lib/python3.12/site-packages')
import gmsh
import numpy as np
import networkx as nx
from collections import defaultdict

MESH_PATH = '/root/repos/block_structured_meshing/T1_9_hub_remeshed.msh'
OUTPUT_PATH = '/root/repos/block_structured_meshing/T1_9_hub_coarse_blocks.msh'
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

def build_face_graph(faces):
    """Build dual graph where each face is a node and adjacent faces are connected."""
    # Build edge-to-faces mapping
    edge_to_faces = defaultdict(list)
    for face_idx in range(len(faces)):
        face = faces[face_idx]
        for i in range(4):
            edge = tuple(sorted([face[i], face[(i+1)%4]]))
            edge_to_faces[edge].append(face_idx)
    
    # Build adjacency
    adj = defaultdict(set)
    for edge, face_list in edge_to_faces.items():
        if len(face_list) == 2:
            f1, f2 = face_list
            adj[f1].add(f2)
            adj[f2].add(f1)
    
    return adj, edge_to_faces

def get_block_boundary_edges(faces, face_indices, edge_to_faces):
    """Get boundary edges of a block (faces that are only in this block)."""
    block_edges = defaultdict(int)
    for fi in face_indices:
        face = faces[fi]
        for i in range(4):
            edge = tuple(sorted([face[i], face[(i+1)%4]]))
            block_edges[edge] += 1
    
    # Boundary edges are those that appear only once
    boundary = [edge for edge, count in block_edges.items() if count == 1]
    return boundary

def can_merge(faces, block1, block2, edge_to_faces):
    """Check if two blocks can be merged and result is a valid quad."""
    combined = set(block1) | set(block2)
    
    # Get boundary of combined block
    boundary = get_block_boundary_edges(faces, combined, edge_to_faces)
    
    # A valid quad block should have exactly 4 boundary edges
    # But actually, a block can be an n-gon. For transfinite interpolation,
    # we need exactly 4 sides.
    
    # Count boundary edges connected to each boundary vertex
    boundary_vertices = set()
    for u, v in boundary:
        boundary_vertices.add(u)
        boundary_vertices.add(v)
    
    # Build boundary graph
    bG = nx.Graph()
    bG.add_edges_from(boundary)
    
    # A valid quad block should have 4 corners (vertices with degree != 2 in boundary)
    corners = [n for n in bG.nodes() if bG.degree(n) != 2]
    
    # For a quad block, we need exactly 4 corners
    return len(corners) == 4, len(corners)

def merge_faces_greedy(faces, target_blocks):
    """Greedy merge of faces to reduce to target number of blocks."""
    adj, edge_to_faces = build_face_graph(faces)
    
    # Start with each face as its own block
    blocks = [{i} for i in range(len(faces))]
    block_id = {i: i for i in range(len(faces))}
    
    # Merge loop
    while len(blocks) > target_blocks:
        # Find best merge (smallest resulting block with 4 corners)
        best_merge = None
        best_score = float('inf')
        
        # Check all adjacent pairs
        for block_idx in range(len(blocks)):
            for face_idx in blocks[block_idx]:
                for neighbor in adj[face_idx]:
                    neighbor_block = block_id[neighbor]
                    if neighbor_block == block_idx:
                        continue
                    
                    can_merge_result, n_corners = can_merge(faces, blocks[block_idx], blocks[neighbor_block], edge_to_faces)
                    if can_merge_result:
                        # Score: prefer merges that create smaller blocks
                        score = len(blocks[block_idx]) + len(blocks[neighbor_block])
                        if score < best_score:
                            best_score = score
                            best_merge = (block_idx, neighbor_block)
        
        if best_merge is None:
            print(f"No valid merge found. Stopping at {len(blocks)} blocks.")
            break
        
        # Merge blocks
        b1, b2 = best_merge
        merged = blocks[b1] | blocks[b2]
        blocks[b1] = merged
        blocks.pop(b2)
        
        # Update block_id
        for face_idx in merged:
            block_id[face_idx] = b1
        
        # Renumber
        new_blocks = []
        new_id = {}
        for i, block in enumerate(blocks):
            if block:
                new_blocks.append(block)
                for face_idx in block:
                    new_id[face_idx] = i
        
        blocks = new_blocks
        block_id = new_id
    
    return blocks

def get_block_vertices_and_edges(faces, block_faces, edge_to_faces):
    """Get vertices and edges for a block."""
    # Get boundary edges
    boundary = get_block_boundary_edges(faces, block_faces, edge_to_faces)
    
    # Order boundary edges into a cycle
    bG = nx.Graph()
    bG.add_edges_from(boundary)
    
    # Find corners
    corners = [n for n in bG.nodes() if bG.degree(n) != 2]
    
    # If not 4 corners, this is not a quad block
    if len(corners) != 4:
        return None, None
    
    # Order boundary cycle
    try:
        cycle = nx.find_cycle(bG)
        ordered_vertices = [cycle[0][0]] + [v for _, v in cycle]
        
        # Extract edges (the boundary edges themselves)
        edges = boundary
        
        return ordered_vertices, edges
    except:
        return None, None

def main():
    print("Loading mesh...")
    vertices, surfaces = load_mesh(MESH_PATH)
    
    for s_tag, faces in surfaces.items():
        print(f"\n=== Surface {s_tag} ===")
        print(f"Faces: {len(faces)}")
        
        # Build dual graph
        adj, edge_to_faces = build_face_graph(faces)
        
        # Merge faces
        print(f"Merging to target {TARGET_BLOCKS} blocks...")
        blocks = merge_faces_greedy(faces, TARGET_BLOCKS)
        print(f"Result: {len(blocks)} blocks")
        
        # Check each block
        valid_blocks = 0
        for i, block in enumerate(blocks):
            ordered_vertices, edges = get_block_vertices_and_edges(faces, block, edge_to_faces)
            if ordered_vertices is not None:
                valid_blocks += 1
                if i < 5:
                    print(f"  Block {i}: {len(block)} faces, {len(ordered_vertices)} vertices")
        
        print(f"Valid quad blocks: {valid_blocks}/{len(blocks)}")
        
        # Export
        export_blocks(vertices, faces, blocks, edge_to_faces, OUTPUT_PATH)
    
    print("\nDone.")

def export_blocks(vertices, faces, blocks, edge_to_faces, path):
    """Export coarse blocks as MSH file."""
    all_points = []
    all_elements = []
    point_offset = 1
    element_id = 1
    
    # Collect unique vertices
    used_vertices = set()
    for block in blocks:
        boundary = get_block_boundary_edges(faces, block, edge_to_faces)
        for u, v in boundary:
            used_vertices.add(u)
            used_vertices.add(v)
    
    # Map vertex index to point ID
    vertex_to_id = {}
    for v_idx in sorted(used_vertices):
        vertex_to_id[v_idx] = point_offset
        all_points.append((point_offset, vertices[v_idx][0], vertices[v_idx][1], vertices[v_idx][2]))
        point_offset += 1
    
    # Export block edges
    for block in blocks:
        boundary = get_block_boundary_edges(faces, block, edge_to_faces)
        bG = nx.Graph()
        bG.add_edges_from(boundary)
        try:
            cycle = nx.find_cycle(bG)
            for u, v in cycle:
                all_elements.append({
                    'id': element_id,
                    'type': 1,
                    'nodes': [vertex_to_id[u], vertex_to_id[v]]
                })
                element_id += 1
        except:
            pass
    
    # Export block faces as quads
    for block in blocks:
        boundary = get_block_boundary_edges(faces, block, edge_to_faces)
        bG = nx.Graph()
        bG.add_edges_from(boundary)
        try:
            corners = [n for n in bG.nodes() if bG.degree(n) != 2]
            if len(corners) == 4:
                # Find corner ordering
                # Start from first corner, trace boundary
                ordered = [corners[0]]
                visited = {corners[0]}
                curr = corners[0]
                while len(ordered) < 4:
                    neighbors = [n for n in bG.neighbors(curr) if n not in visited]
                    if neighbors:
                        curr = neighbors[0]
                        ordered.append(curr)
                        visited.add(curr)
                    else:
                        break
                
                if len(ordered) == 4:
                    all_elements.append({
                        'id': element_id,
                        'type': 3,
                        'nodes': [vertex_to_id[v] for v in ordered]
                    })
                    element_id += 1
        except:
            pass
    
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
    
    print(f"\nExported to {path}")

if __name__ == '__main__':
    main()
