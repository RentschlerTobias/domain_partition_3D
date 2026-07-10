#!/usr/bin/env python3
"""
Analyze hex mesh structure around blade and create VTK visualization.
Understand how the 44,800 hex elements are organized in 5 regions.
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
                # 8 node indices
                nodes_idx = [int(p) for p in parts[3 + num_tags:]]
                hex_elements.append({
                    'tag': tag,
                    'phys_tag': phys_tag,
                    'geom_tag': geom_tag,
                    'nodes': nodes_idx
                })
    
    return nodes, hex_elements

def analyze_hex_regions(nodes, hex_elements):
    """Analyze each hex region's geometry."""
    regions = defaultdict(list)
    
    for hex_elem in hex_elements:
        regions[hex_elem['geom_tag']].append(hex_elem)
    
    print("\n=== HEX REGION ANALYSIS ===")
    print(f"Total hex elements: {len(hex_elements)}\n")
    
    region_stats = {}
    for geom_tag, elems in sorted(regions.items()):
        # Compute bounding box and center
        all_coords = []
        for elem in elems:
            for nid in elem['nodes']:
                all_coords.append(nodes[nid])
        all_coords = np.array(all_coords)
        
        bb_min = all_coords.min(axis=0)
        bb_max = all_coords.max(axis=0)
        center = all_coords.mean(axis=0)
        
        # Count unique nodes
        unique_nodes = set()
        for elem in elems:
            unique_nodes.update(elem['nodes'])
        
        print(f"Region geom_tag={geom_tag}:")
        print(f"  Elements: {len(elems)}")
        print(f"  Unique nodes: {len(unique_nodes)}")
        print(f"  Bounding box:")
        print(f"    X: [{bb_min[0]:.4f}, {bb_max[0]:.4f}]")
        print(f"    Y: [{bb_min[1]:.4f}, {bb_max[1]:.4f}]")
        print(f"    Z: [{bb_min[2]:.4f}, {bb_max[2]:.4f}]")
        print(f"  Center: ({center[0]:.4f}, {center[1]:.4f}, {center[2]:.4f})")
        
        # Compute radial extent
        r = np.sqrt(all_coords[:,0]**2 + all_coords[:,1]**2)
        print(f"  Radius: [{r.min():.4f}, {r.max():.4f}]")
        
        # Angular extent
        theta = np.degrees(np.arctan2(all_coords[:,1], all_coords[:,0]))
        print(f"  Theta: [{theta.min():.1f}°, {theta.max():.1f}°]")
        print()
        
        region_stats[geom_tag] = {
            'num_elements': len(elems),
            'num_nodes': len(unique_nodes),
            'bbox': {
                'min': bb_min.tolist(),
                'max': bb_max.tolist()
            },
            'center': center.tolist(),
            'radius': [float(r.min()), float(r.max())],
            'theta': [float(theta.min()), float(theta.max())]
        }
    
    return region_stats

def write_vtk_unstructured(nodes, hex_elements, filename):
    """Write hex mesh as VTK unstructured grid."""
    # Collect unique nodes used by hex elements
    used_node_tags = set()
    for elem in hex_elements:
        used_node_tags.update(elem['nodes'])
    
    used_node_tags = sorted(used_node_tags)
    node_index_map = {tag: idx for idx, tag in enumerate(used_node_tags)}
    
    with open(filename, 'w') as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("T1_9 Hex Mesh Regions\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        
        # Points
        f.write(f"POINTS {len(used_node_tags)} float\n")
        for tag in used_node_tags:
            x, y, z = nodes[tag]
            f.write(f"{x} {y} {z}\n")
        
        # Cells
        f.write(f"\nCELLS {len(hex_elements)} {len(hex_elements) * 9}\n")
        for elem in hex_elements:
            f.write("8 ")
            for nid in elem['nodes']:
                f.write(f"{node_index_map[nid]} ")
            f.write("\n")
        
        # Cell types
        f.write(f"\nCELL_TYPES {len(hex_elements)}\n")
        for _ in hex_elements:
            f.write("12\n")  # VTK_HEXAHEDRON = 12
        
        # Cell data: geom_tag
        f.write(f"\nCELL_DATA {len(hex_elements)}\n")
        f.write("SCALARS geom_tag int 1\n")
        f.write("LOOKUP_TABLE default\n")
        for elem in hex_elements:
            f.write(f"{elem['geom_tag']}\n")
        
        # Also write element ID
        f.write("\nSCALARS element_id int 1\n")
        f.write("LOOKUP_TABLE default\n")
        for elem in hex_elements:
            f.write(f"{elem['tag']}\n")
    
    print(f"VTK file written: {filename}")
    print(f"  Nodes: {len(used_node_tags)}")
    print(f"  Hex elements: {len(hex_elements)}")

def extract_region_vtk(nodes, hex_elements, geom_tag, filename):
    """Extract single region as VTK."""
    region_elems = [e for e in hex_elements if e['geom_tag'] == geom_tag]
    write_vtk_unstructured(nodes, region_elems, filename)

def analyze_block_structure(nodes, hex_elements):
    """Try to understand the implicit block structure."""
    # Look at element connectivity to find structured grid dimensions
    regions = defaultdict(list)
    for elem in hex_elements:
        regions[elem['geom_tag']].append(elem)
    
    print("\n=== BLOCK STRUCTURE ANALYSIS ===\n")
    
    for geom_tag, elems in sorted(regions.items()):
        print(f"Region {geom_tag} ({len(elems)} elements):")
        
        # Find structured dimensions by analyzing connectivity
        # Build node-to-elements mapping
        node_to_elems = defaultdict(set)
        for elem in elems:
            for nid in elem['nodes']:
                node_to_elems[nid].add(elem['tag'])
        
        # Count element neighbors for each element
        elem_map = {e['tag']: e for e in elems}
        neighbor_counts = []
        
        for elem in elems:
            neighbors = set()
            for nid in elem['nodes']:
                for other_tag in node_to_elems[nid]:
                    if other_tag != elem['tag']:
                        neighbors.add(other_tag)
            neighbor_counts.append(len(neighbors))
        
        neighbor_counts = np.array(neighbor_counts)
        print(f"  Neighbor count: min={neighbor_counts.min()}, max={neighbor_counts.max()}, mean={neighbor_counts.mean():.1f}")
        
        # In a perfect structured grid, interior elements have 6 neighbors
        # (one for each face)
        perfect = np.sum(neighbor_counts == 6)
        print(f"  Elements with 6 neighbors (interior): {perfect} ({100*perfect/len(elems):.1f}%)")
        
        # Try to estimate grid dimensions
        # For a structured grid with N elements, dimensions might be nx*ny*nz = N
        # We need to find factors
        N = len(elems)
        factors = []
        for i in range(1, int(np.sqrt(N)) + 1):
            if N % i == 0:
                factors.append((i, N // i))
        
        print(f"  Possible 2D factorizations of {N}: {factors[:10]}")
        print()

def main():
    msh_file = 'T1_9/T1_9_ru_gridGmsh.msh'
    
    print("Parsing mesh...")
    nodes, hex_elements = parse_msh(msh_file)
    print(f"Loaded {len(nodes)} nodes, {len(hex_elements)} hex elements")
    
    # Analyze regions
    region_stats = analyze_hex_regions(nodes, hex_elements)
    
    # Analyze block structure
    analyze_block_structure(nodes, hex_elements)
    
    # Write full VTK
    write_vtk_unstructured(nodes, hex_elements, 'T1_9_hex_mesh_all.vtk')
    
    # Write individual region VTKs
    for geom_tag in sorted(region_stats.keys()):
        extract_region_vtk(nodes, hex_elements, geom_tag, f'T1_9_hex_region_{geom_tag}.vtk')
    
    # Save stats
    with open('hex_region_stats.json', 'w') as f:
        json.dump(region_stats, f, indent=2)
    print("Stats saved to hex_region_stats.json")

if __name__ == '__main__':
    main()
