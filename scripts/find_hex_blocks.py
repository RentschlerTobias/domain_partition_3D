#!/usr/bin/env python3
"""
Analyze hex mesh to find implicit block structure.
Look for structured sub-grids within each region.
"""

import numpy as np
import json
from collections import defaultdict

def parse_msh(filename):
    """Parse Gmsh 2.2 ASCII file, return nodes and hex elements."""
    nodes = {}
    hex_elements = []
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Parse nodes
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
            nodes[tag] = (x, y, z)
    
    # Parse hex elements
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
            if elem_type == 5:  # Hexahedron
                tag = int(parts[0])
                num_tags = int(parts[2])
                phys_tag = int(parts[3]) if num_tags >= 1 else 0
                geom_tag = int(parts[4]) if num_tags >= 2 else 0
                nodes_idx = [int(p) for p in parts[3 + num_tags:]]
                hex_elements.append({
                    'tag': tag,
                    'phys_tag': phys_tag,
                    'geom_tag': geom_tag,
                    'nodes': nodes_idx
                })
    
    return nodes, hex_elements

def find_structured_dimensions(elems, nodes):
    """Try to find structured grid dimensions by analyzing regular connectivity."""
    # Build node-to-elements
    node_to_elems = defaultdict(list)
    for elem in elems:
        for nid in elem['nodes']:
            node_to_elems[nid].append(elem)
    
    # For each element, count face neighbors (elements sharing exactly 4 nodes = a face)
    elem_map = {e['tag']: e for e in elems}
    face_neighbors = {}
    
    for elem in elems:
        # Get all faces (each face is 4 nodes)
        n = elem['nodes']
        faces = [
            frozenset([n[0], n[1], n[2], n[3]]),  # bottom
            frozenset([n[4], n[5], n[6], n[7]]),  # top
            frozenset([n[0], n[1], n[5], n[4]]),  # front
            frozenset([n[2], n[3], n[7], n[6]]),  # back
            frozenset([n[0], n[3], n[7], n[4]]),  # left
            frozenset([n[1], n[2], n[6], n[5]]),  # right
        ]
        
        neighbors = []
        for face in faces:
            # Find other elements sharing this face
            # Get one node from face
            sample_node = list(face)[0]
            for other in node_to_elems[sample_node]:
                if other['tag'] != elem['tag']:
                    if len(face.intersection(other['nodes'])) == 4:
                        neighbors.append(other['tag'])
        
        face_neighbors[elem['tag']] = neighbors
    
    # Count how many face neighbors each element has
    neighbor_counts = [len(face_neighbors[e['tag']]) for e in elems]
    
    return face_neighbors, neighbor_counts

def analyze_region_blocks(nodes, hex_elements, geom_tag):
    """Analyze a single region to find block structure."""
    elems = [e for e in hex_elements if e['geom_tag'] == geom_tag]
    
    print(f"\n=== Region {geom_tag} ({len(elems)} elements) ===")
    
    face_neighbors, neighbor_counts = find_structured_dimensions(elems, nodes)
    
    nc = np.array(neighbor_counts)
    print(f"Face neighbor count: min={nc.min()}, max={nc.max()}, mean={nc.mean():.1f}")
    print(f"Elements with 6 face neighbors: {np.sum(nc == 6)} ({100*np.sum(nc == 6)/len(elems):.1f}%)")
    print(f"Elements with 5 face neighbors: {np.sum(nc == 5)}")
    print(f"Elements with 4 face neighbors: {np.sum(nc == 4)}")
    print(f"Elements with <4 face neighbors: {np.sum(nc < 4)}")
    
    # Elements with 6 face neighbors are interior to a structured block
    # Elements with <6 are on boundaries between blocks
    
    # Try to find connected components of elements with 6 neighbors
    interior_elems = [e for e in elems if len(face_neighbors[e['tag']]) == 6]
    print(f"Interior elements: {len(interior_elems)}")
    
    if len(interior_elems) == 0:
        print("No perfect interior elements found - grid may be unstructured or O-grid type")
        
        # Try with 5 neighbors (one face on boundary)
        interior_elems_5 = [e for e in elems if len(face_neighbors[e['tag']]) >= 5]
        print(f"Elements with >=5 face neighbors: {len(interior_elems_5)}")
    
    # Try to estimate dimensions by looking at element centers
    centers = []
    for elem in elems:
        coords = np.array([nodes[n] for n in elem['nodes']])
        centers.append(coords.mean(axis=0))
    centers = np.array(centers)
    
    # Cluster centers by z-layer
    z_values = np.unique(np.round(centers[:,2], 4))
    print(f"Unique Z-layers (rounded): {len(z_values)}")
    print(f"Z-range: [{z_values.min():.4f}, {z_values.max():.4f}]")
    
    # Count elements per z-layer
    elems_per_layer = []
    for z in z_values:
        count = np.sum(np.abs(centers[:,2] - z) < 0.001)
        elems_per_layer.append(count)
    
    print(f"Elements per layer: min={min(elems_per_layer)}, max={max(elems_per_layer)}")
    print(f"Unique layer counts: {sorted(set(elems_per_layer))}")
    
    # If all layers have same count, we can estimate 2D dimensions
    if len(set(elems_per_layer)) == 1:
        layer_count = elems_per_layer[0]
        print(f"Consistent {layer_count} elements per layer")
        
        # Find factors
        factors = []
        for i in range(1, int(np.sqrt(layer_count)) + 1):
            if layer_count % i == 0:
                factors.append((i, layer_count // i))
        print(f"Possible 2D dimensions: {factors}")
    
    return {
        'geom_tag': geom_tag,
        'num_elements': len(elems),
        'num_layers': len(z_values),
        'elements_per_layer': elems_per_layer,
        'neighbor_stats': {
            'min': int(nc.min()),
            'max': int(nc.max()),
            'mean': float(nc.mean()),
            'six_neighbors': int(np.sum(nc == 6))
        }
    }

def write_block_vtk(nodes, hex_elements, filename, block_geom_tags=None):
    """Write selected regions as VTK with block IDs."""
    if block_geom_tags:
        elems = [e for e in hex_elements if e['geom_tag'] in block_geom_tags]
    else:
        elems = hex_elements
    
    # Collect unique nodes
    used_node_tags = set()
    for elem in elems:
        used_node_tags.update(elem['nodes'])
    
    used_node_tags = sorted(used_node_tags)
    node_index_map = {tag: idx for idx, tag in enumerate(used_node_tags)}
    
    with open(filename, 'w') as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("T1_9 Hex Blocks\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        
        f.write(f"POINTS {len(used_node_tags)} float\n")
        for tag in used_node_tags:
            x, y, z = nodes[tag]
            f.write(f"{x} {y} {z}\n")
        
        f.write(f"\nCELLS {len(elems)} {len(elems) * 9}\n")
        for elem in elems:
            f.write("8 ")
            for nid in elem['nodes']:
                f.write(f"{node_index_map[nid]} ")
            f.write("\n")
        
        f.write(f"\nCELL_TYPES {len(elems)}\n")
        for _ in elems:
            f.write("12\n")
        
        # Block ID
        f.write(f"\nCELL_DATA {len(elems)}\n")
        f.write("SCALARS block_id int 1\n")
        f.write("LOOKUP_TABLE default\n")
        for elem in elems:
            f.write(f"{elem['geom_tag']}\n")
    
    print(f"\nBlock VTK written: {filename}")
    print(f"  Regions included: {sorted(set(e['geom_tag'] for e in elems))}")
    print(f"  Total elements: {len(elems)}")

def main():
    msh_file = 'T1_9/T1_9_ru_gridGmsh.msh'
    
    print("Parsing mesh...")
    nodes, hex_elements = parse_msh(msh_file)
    print(f"Loaded {len(nodes)} nodes, {len(hex_elements)} hex elements")
    
    # Analyze each region
    all_stats = {}
    for geom_tag in sorted(set(e['geom_tag'] for e in hex_elements)):
        stats = analyze_region_blocks(nodes, hex_elements, geom_tag)
        all_stats[geom_tag] = stats
    
    # Write combined VTK with all regions
    write_block_vtk(nodes, hex_elements, 'T1_9_hex_blocks.vtk')
    
    # Write individual block VTKs
    for geom_tag in sorted(all_stats.keys()):
        write_block_vtk(nodes, hex_elements, f'T1_9_hex_block_{geom_tag}.vtk', [geom_tag])
    
    # Save stats
    with open('hex_block_stats.json', 'w') as f:
        json.dump(all_stats, f, indent=2)
    print("\nStats saved to hex_block_stats.json")

if __name__ == '__main__':
    main()
