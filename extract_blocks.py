import gmsh
import numpy as np
from collections import defaultdict

OUT_DIR = "/root/repos/block_structured_meshing"

def extract_coarse_blocks(msh_file, name, out_msh):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    
    gmsh.open(msh_file)
    
    # get nodes
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    node_coords = np.array(node_coords).reshape(-1, 3)
    node_map = {int(tag): i for i, tag in enumerate(node_tags)}
    
    # get quads
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(2)
    quads = []
    for etype, etags, enodes in zip(elem_types, elem_tags, elem_node_tags):
        if etype == 3:  # 4-node quad
            n_node_per_elem = 4
            for i in range(len(etags)):
                n = enodes[i*n_node_per_elem:(i+1)*n_node_per_elem]
                quads.append([int(x) for x in n])
    
    print(f"[{name}] {len(quads)} quads, {len(node_tags)} nodes")
    
    # build adjacency
    vertex_edges = defaultdict(set)
    edge_quads = defaultdict(list)
    quad_edges = []
    
    for qi, quad in enumerate(quads):
        q_edges = []
        for i in range(4):
            v1 = quad[i]
            v2 = quad[(i+1)%4]
            edge = tuple(sorted([v1, v2]))
            edge_quads[edge].append(qi)
            vertex_edges[v1].add(v2)
            vertex_edges[v2].add(v1)
            q_edges.append(edge)
        quad_edges.append(q_edges)
    
    # compute valence
    valence = {}
    for v in vertex_edges:
        valence[v] = len(vertex_edges[v])
    
    # identify singular vertices
    singular = set()
    boundary_vertices = set()
    for v in vertex_edges:
        is_boundary = False
        for nb in vertex_edges[v]:
            edge = tuple(sorted([v, nb]))
            if len(edge_quads[edge]) == 1:
                is_boundary = True
                break
        
        if is_boundary:
            boundary_vertices.add(v)
        else:
            if valence[v] != 4:
                singular.add(v)
    
    print(f"[{name}] singular vertices: {len(singular)}")
    print(f"[{name}] boundary vertices: {len(boundary_vertices)}")
    
    # trace separatrices
    # compute cyclic order of edges around each vertex
    vertex_order = {}
    for v in vertex_edges:
        neighbors = list(set(vertex_edges[v]))
        if len(neighbors) < 2:
            continue
        
        # get normal vector
        normals = []
        for nb in neighbors:
            edge = tuple(sorted([v, nb]))
            for qi in edge_quads[edge]:
                quad = quads[qi]
                idx = [node_map[x] for x in quad]
                pts = node_coords[idx]
                v1 = pts[1] - pts[0]
                v2 = pts[2] - pts[1]
                n = np.cross(v1, v2)
                if np.linalg.norm(n) > 1e-10:
                    n = n / np.linalg.norm(n)
                    normals.append(n)
        
        if len(normals) == 0:
            continue
        
        normal = np.mean(normals, axis=0)
        normal = normal / (np.linalg.norm(normal) + 1e-10)
        
        v_pos = node_coords[node_map[v]]
        projected = []
        for nb in neighbors:
            nb_pos = node_coords[node_map[nb]]
            vec = nb_pos - v_pos
            vec = vec - np.dot(vec, normal) * normal
            if np.linalg.norm(vec) > 1e-10:
                projected.append((nb, vec))
        
        if len(projected) < 2:
            continue
        
        ref = projected[0][1]
        ref = ref / np.linalg.norm(ref)
        
        angles = []
        for nb, vec in projected:
            vec = vec / np.linalg.norm(vec)
            angle = np.arctan2(np.dot(np.cross(ref, vec), normal), np.dot(ref, vec))
            if angle < 0:
                angle += 2 * np.pi
            angles.append((angle, nb))
        
        angles.sort()
        vertex_order[v] = [nb for _, nb in angles]
    
    # identify separatrix edges
    separatrix_edges = set()
    visited_paths = set()
    
    for v0 in singular:
        if v0 not in vertex_order:
            continue
        
        order = vertex_order[v0]
        n_edges = len(order)
        
        for i in range(n_edges):
            v1 = order[i]
            
            # trace path
            path = [v0, v1]
            current = v1
            prev = v0
            
            max_iter = 1000
            iter_count = 0
            while iter_count < max_iter:
                iter_count += 1
                
                if current in singular:
                    break
                
                if current not in vertex_order:
                    break
                
                order_curr = vertex_order[current]
                n_curr = len(order_curr)
                
                if n_curr < 2:
                    break
                
                if prev not in order_curr:
                    break
                
                idx = order_curr.index(prev)
                opposite_idx = (idx + n_curr // 2) % n_curr
                next_v = order_curr[opposite_idx]
                
                if next_v == path[-1]:
                    break
                
                path.append(next_v)
                prev = current
                current = next_v
                
                if current in path[:-1]:
                    break
            
            # store edges
            path_tuple = tuple(sorted([tuple(path), tuple(reversed(path))])[0])
            if path_tuple not in visited_paths:
                visited_paths.add(path_tuple)
                for j in range(len(path)-1):
                    edge = tuple(sorted([path[j], path[j+1]]))
                    separatrix_edges.add(edge)
    
    print(f"[{name}] separatrix edges: {len(separatrix_edges)}")
    
    # flood-fill to find blocks
    # a block is a set of quads bounded by separatrix edges and boundary edges
    visited_quads = [False] * len(quads)
    blocks = []
    
    for start_qi in range(len(quads)):
        if visited_quads[start_qi]:
            continue
        
        # flood fill
        block_quads = []
        stack = [start_qi]
        visited_quads[start_qi] = True
        
        while stack:
            qi = stack.pop()
            block_quads.append(qi)
            
            # check neighbors
            for edge in quad_edges[qi]:
                if len(edge_quads[edge]) == 2:
                    other_qi = edge_quads[edge][0] if edge_quads[edge][1] == qi else edge_quads[edge][1]
                    
                    if not visited_quads[other_qi]:
                        # check if edge is a separatrix edge
                        if edge not in separatrix_edges:
                            visited_quads[other_qi] = True
                            stack.append(other_qi)
        
        blocks.append(block_quads)
    
    print(f"[{name}] blocks found: {len(blocks)}")
    
    # analyze block boundaries
    blocks_4 = []
    for block_quads in blocks:
        # find boundary edges of this block
        block_boundary_edges = set()
        for qi in block_quads:
            for edge in quad_edges[qi]:
                # edge is on block boundary if it's a separatrix or boundary or shared with another block
                if edge in separatrix_edges:
                    block_boundary_edges.add(edge)
                elif len(edge_quads[edge]) == 1:
                    block_boundary_edges.add(edge)
                else:
                    # shared with another block
                    other_qi = edge_quads[edge][0] if edge_quads[edge][1] == qi else edge_quads[edge][1]
                    if other_qi not in block_quads:
                        block_boundary_edges.add(edge)
        
        # count boundary segments
        # a boundary segment is a chain of edges
        # count vertices on boundary
        boundary_vertices = set()
        for edge in block_boundary_edges:
            boundary_vertices.add(edge[0])
            boundary_vertices.add(edge[1])
        
        # count corners
        # corners are boundary vertices where direction changes
        n_corners = 0
        for v in boundary_vertices:
            # count boundary edges at this vertex
            n_bnd_edges = 0
            for nb in vertex_edges[v]:
                edge = tuple(sorted([v, nb]))
                if edge in block_boundary_edges:
                    n_bnd_edges += 1
            
            if n_bnd_edges == 2:
                n_corners += 1
        
        # heuristic: if block has 4 corners, it's a 4-sided block
        if n_corners == 4:
            blocks_4.append(block_quads)
    
    print(f"[{name}] 4-sided blocks: {len(blocks_4)}")
    
    # create coarse block mesh
    gmsh.clear()
    
    # collect all vertices on block boundaries
    all_boundary_vertices = set()
    for block_quads in blocks_4:
        for qi in block_quads:
            for v in quads[qi]:
                all_boundary_vertices.add(v)
    
    all_boundary_vertices = sorted(list(all_boundary_vertices))
    node_tag_map = {}
    for i, v in enumerate(all_boundary_vertices):
        tag = i + 1
        node_tag_map[v] = tag
        idx = node_map[v]
        coord = node_coords[idx]
        gmsh.model.addDiscreteEntity(0, tag)
        gmsh.model.mesh.addNodes(0, tag, [tag], list(coord))
    
    # add quads for blocks
    surf = gmsh.model.addDiscreteEntity(2, 1)
    quad_nodes = []
    for block_quads in blocks_4:
        # compute block center
        center = np.zeros(3)
        for qi in block_quads:
            for v in quads[qi]:
                idx = node_map[v]
                center += node_coords[idx]
        center /= (len(block_quads) * 4)
        
        # find 4 corner vertices
        # corners are boundary vertices with 2 boundary edges
        block_boundary_edges = set()
        for qi in block_quads:
            for edge in quad_edges[qi]:
                if edge in separatrix_edges or len(edge_quads[edge]) == 1:
                    block_boundary_edges.add(edge)
                else:
                    other_qi = edge_quads[edge][0] if edge_quads[edge][1] == qi else edge_quads[edge][1]
                    if other_qi not in block_quads:
                        block_boundary_edges.add(edge)
        
        corners = []
        for v in all_boundary_vertices:
            if v not in [q for qi in block_quads for q in quads[qi]]:
                continue
            
            n_bnd_edges = 0
            for nb in vertex_edges[v]:
                edge = tuple(sorted([v, nb]))
                if edge in block_boundary_edges:
                    n_bnd_edges += 1
            
            if n_bnd_edges == 2:
                corners.append(v)
        
        if len(corners) != 4:
            continue
        
        # order corners
        # compute angles from center
        corner_angles = []
        for v in corners:
            idx = node_map[v]
            vec = node_coords[idx] - center
            angle = np.arctan2(vec[1], vec[0])
            if angle < 0:
                angle += 2 * np.pi
            corner_angles.append((angle, v))
        
        corner_angles.sort()
        ordered_corners = [v for _, v in corner_angles]
        
        quad_nodes.extend([node_tag_map[v] for v in ordered_corners])
    
    if quad_nodes:
        gmsh.model.mesh.addElementsByType(1, 3, [], quad_nodes)
    
    gmsh.write(out_msh)
    print(f"[{name}] saved coarse blocks to {out_msh}")
    
    gmsh.finalize()

if __name__ == "__main__":
    extract_coarse_blocks(f"{OUT_DIR}/T1_9_hub_quad.msh", "hub", f"{OUT_DIR}/T1_9_hub_blocks.msh")
    extract_coarse_blocks(f"{OUT_DIR}/T1_9_shroud_quad.msh", "shroud", f"{OUT_DIR}/T1_9_shroud_blocks.msh")
