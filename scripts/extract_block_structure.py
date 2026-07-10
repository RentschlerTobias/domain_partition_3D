#!/usr/bin/env python3
"""
Extract actual block structure from hex mesh by finding structured sub-grids.
Use connectivity to identify i,j,k directions and dimensions.
"""

import numpy as np
import json
from collections import defaultdict, deque

def parse_msh(filename):
    nodes = {}
    hex_elements = []
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    in_nodes = False
    for line in lines:
        line = line.strip()
        if line == '$Nodes':
            in_nodes = True
            continue
        if line == '$EndNodes':
            break
        if in_nodes:
            if line.isdigit():
                continue
            parts = line.split()
            tag = int(parts[0])
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            nodes[tag] = np.array([x, y, z])
    
    in_elements = False
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line == '$Elements':
            in_elements = True
            continue
        if line == '$EndElements':
            break
        if in_elements:
            if line.isdigit():
                continue
            parts = line.split()
            elem_type = int(parts[1])
            if elem_type == 5:
                tag = int(parts[0])
                num_tags = int(parts[2])
                geom_tag = int(parts[4]) if num_tags >= 2 else 0
                nodes_idx = [int(p) for p in parts[3 + num_tags:]]
                hex_elements.append({
                    'tag': tag,
                    'geom_tag': geom_tag,
                    'nodes': nodes_idx
                })
    
    return nodes, hex_elements

def build_face_neighbors(hex_elements):
    """Build face neighbor connectivity."""
    # Map face (4 nodes) to elements
    face_to_elems = defaultdict(list)
    
    for elem in hex_elements:
        n = elem['nodes']
        # 6 faces of hexahedron
        faces = [
            frozenset([n[0], n[1], n[2], n[3]]),
            frozenset([n[4], n[5], n[6], n[7]]),
            frozenset([n[0], n[1], n[5], n[4]]),
            frozenset([n[2], n[3], n[7], n[6]]),
            frozenset([n[0], n[3], n[7], n[4]]),
            frozenset([n[1], n[2], n[6], n[5]]),
        ]
        for face in faces:
            face_to_elems[face].append(elem['tag'])
    
    # Build neighbor graph
    neighbors = defaultdict(set)
    for face, elems in face_to_elems.items():
        if len(elems) == 2:
            neighbors[elems[0]].add(elems[1])
            neighbors[elems[1]].add(elems[0])
    
    return neighbors

def find_connected_components(hex_elements, neighbors):
    """Find connected components in the hex mesh."""
    elem_tags = {e['tag'] for e in hex_elements}
    visited = set()
    components = []
    
    for start_tag in elem_tags:
        if start_tag in visited:
            continue
        
        component = []
        queue = deque([start_tag])
        visited.add(start_tag)
        
        while queue:
            tag = queue.popleft()
            component.append(tag)
            
            for neighbor in neighbors[tag]:
                if neighbor not in visited and neighbor in elem_tags:
                    visited.add(neighbor)
                    queue.append(neighbor)
        
        components.append(component)
    
    return components

def estimate_dimensions(elems, nodes, neighbors):
    """Estimate structured grid dimensions by analyzing element chains."""
    elem_map = {e['tag']: e for e in elems}
    
    # Pick an interior element with 6 neighbors
    interior_elem = None
    for e in elems:
        if len(neighbors[e['tag']]) == 6:
            interior_elem = e
            break
    
    if interior_elem is None:
        return None
    
    # Get the 6 neighbors
    neighbor_tags = list(neighbors[interior_elem['tag']])
    
    # Compute vectors from center of interior element to centers of neighbors
    def elem_center(elem):
        return np.mean([nodes[n] for n in elem['nodes']], axis=0)
    
    center = elem_center(interior_elem)
    
    # Vectors to neighbors
    vectors = []
    for nt in neighbor_tags:
        nc = elem_center(elem_map[nt])
        v = nc - center
        vectors.append((nt, v, np.linalg.norm(v)))
    
    # Sort by direction
    vectors.sort(key=lambda x: x[2])
    
    # In a structured grid, opposite neighbors should have opposite vectors
    # Group into 3 pairs (i, j, k directions)
    
    # Simple approach: find 3 orthogonal-ish directions
    used = set()
    directions = []
    
    for i, (nt1, v1, d1) in enumerate(vectors):
        if i in used:
            continue
        
        # Find opposite direction
        best_match = None
        best_dot = -1
        for j, (nt2, v2, d2) in enumerate(vectors):
            if j in used or i == j:
                continue
            dot = np.dot(v1, v2) / (d1 * d2)
            if dot < -0.5 and dot < best_dot:
                best_dot = dot
                best_match = j
        
        if best_match is not None:
            used.add(i)
            used.add(best_match)
            # Average direction
            dir_vec = (v1 - vectors[best_match][1]) / 2
            directions.append((dir_vec, d1))
    
    return {
        'interior_elem': interior_elem['tag'],
        'num_neighbors': len(neighbor_tags),
        'directions': len(directions),
        'vectors': [(v.tolist(), float(d)) for v, d in directions]
    }

def trace_direction(start_elem, direction_vec, elem_map, nodes, neighbors, tol=0.5):
    """Trace along a direction to find number of elements in that direction."""
    def elem_center(elem):
        return np.mean([nodes[n] for n in elem['nodes']], axis=0)
    
    chain = [start_elem['tag']]
    current = start_elem
    
    while True:
        center = elem_center(current)
        best_next = None
        best_dot = -1
        
        for nt in neighbors[current['tag']]:
            if nt in chain:
                continue
            nc = elem_center(elem_map[nt])
            v = nc - center
            v_norm = np.linalg.norm(v)
            if v_norm < 0.001:
                continue
            
            v_unit = v / v_norm
            dir_unit = direction_vec / np.linalg.norm(direction_vec)
            dot = np.dot(v_unit, dir_unit)
            
            if dot > tol and dot > best_dot:
                best_dot = dot
                best_next = nt
        
        if best_next is None:
            break
        
        chain.append(best_next)
        current = elem_map[best_next]
    
    return len(chain)

def main():
    msh_file = 'T1_9/T1_9_ru_gridGmsh.msh'
    
    print("Parsing mesh...")
    nodes, hex_elements = parse_msh(msh_file)
    print(f"Loaded {len(nodes)} nodes, {len(hex_elements)} hex elements")
    
    # Group by geom_tag
    regions = defaultdict(list)
    for elem in hex_elements:
        regions[elem['geom_tag']].append(elem)
    
    # Analyze each region
    all_block_info = {}
    
    for geom_tag in sorted(regions.keys()):
        elems = regions[geom_tag]
        print(f"\n{'='*60}")
        print(f"Region {geom_tag}: {len(elems)} elements")
        print(f"{'='*60}")
        
        # Build neighbors
        neighbors = build_face_neighbors(elems)
        
        # Find connected components
        components = find_connected_components(elems, neighbors)
        print(f"Connected components: {len(components)}")
        for i, comp in enumerate(components):
            print(f"  Component {i+1}: {len(comp)} elements")
        
        # Estimate dimensions
        dims = estimate_dimensions(elems, nodes, neighbors)
        if dims:
            print(f"\nDimension estimation:")
            print(f"  Interior element: {dims['interior_elem']}")
            print(f"  Found {dims['directions']} direction pairs")
            for i, (vec, dist) in enumerate(dims['vectors']):
                print(f"  Direction {i+1}: vec={vec}, dist={dist:.4f}")
        
        # Try to count elements in each direction from a corner element
        # Find element with minimum neighbors (corner)
        corner_elem = min(elems, key=lambda e: len(neighbors[e['tag']]))
        print(f"\nCorner element (min neighbors): {corner_elem['tag']}")
        print(f"  Neighbors: {len(neighbors[corner_elem['tag']])}")
        
        # Try to trace in 3 directions from corner
        if dims and dims['directions'] >= 1:
            elem_map = {e['tag']: e for e in elems}
            
            for i, (vec, _) in enumerate(dims['vectors']):
                count = trace_direction(corner_elem, np.array(vec), elem_map, nodes, neighbors)
                print(f"  Direction {i+1} chain length: {count} elements")
        
        all_block_info[geom_tag] = {
            'num_elements': len(elems),
            'num_components': len(components),
            'component_sizes': [len(c) for c in components],
            'dimension_estimate': dims
        }
    
    # Save info
    with open('hex_block_structure.json', 'w') as f:
        json.dump(all_block_info, f, indent=2, default=lambda x: x.tolist() if hasattr(x, 'tolist') else x)
    print("\nStructure saved to hex_block_structure.json")

if __name__ == '__main__':
    main()
