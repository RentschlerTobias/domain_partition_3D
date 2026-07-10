#!/usr/bin/env python3
"""
Validate dtOO block corners against hex elements in T1_9 .msh file.

This script:
1. Loads the computed block corners from dtOO_blocks_mapped.json
2. Parses hex elements from T1_9/T1_9_ru_gridGmsh.msh
3. Checks which hex elements fall within each block
4. Reports coverage statistics
"""

import json
import numpy as np
from pathlib import Path

def parse_msh_file(filepath):
    """Parse hex elements from Gmsh 2.2 format."""
    nodes = {}
    hex_elements = []
    
    with open(filepath, 'r') as f:
        line = f.readline()
        while line:
            line = line.strip()
            
            if line == '$Nodes':
                n_nodes = int(f.readline().strip())
                for _ in range(n_nodes):
                    parts = f.readline().strip().split()
                    node_id = int(parts[0])
                    coords = [float(parts[1]), float(parts[2]), float(parts[3])]
                    nodes[node_id] = coords
            
            elif line == '$Elements':
                n_elements = int(f.readline().strip())
                for _ in range(n_elements):
                    parts = f.readline().strip().split()
                    elem_type = int(parts[1])
                    if elem_type == 5:  # Hex element
                        n_tags = int(parts[2])
                        node_ids = [int(p) for p in parts[3+n_tags:]]
                        hex_elements.append(node_ids)
            
            line = f.readline()
    
    return nodes, hex_elements

def point_in_hex(point, hex_corners):
    """
    Check if a point is inside a tri-linear hexahedron.
    Uses a simplified bounding box check first, then more precise test.
    """
    hex_corners = np.array(hex_corners)
    point = np.array(point)
    
    # Bounding box check
    min_coords = np.min(hex_corners, axis=0)
    max_coords = np.max(hex_corners, axis=0)
    
    if np.any(point < min_coords - 1e-6) or np.any(point > max_coords + 1e-6):
        return False
    
    return True

def hex_center(node_ids, nodes):
    """Compute center of a hex element."""
    coords = [nodes[nid] for nid in node_ids]
    return np.mean(coords, axis=0)

def main():
    print("Loading dtOO block definitions...")
    with open('dtOO_blocks_mapped.json', 'r') as f:
        block_data = json.load(f)
    
    blocks = block_data['blocks']
    
    print("Parsing T1_9 mesh file...")
    msh_path = Path('T1_9/T1_9_ru_gridGmsh.msh')
    nodes, hex_elements = parse_msh_file(msh_path)
    
    print(f"Found {len(hex_elements)} hex elements")
    print(f"Found {len(nodes)} nodes")
    
    # Compute hex centers
    print("Computing hex element centers...")
    hex_centers = []
    for elem_nodes in hex_elements:
        center = hex_center(elem_nodes, nodes)
        hex_centers.append(center)
    hex_centers = np.array(hex_centers)
    
    # Check which hex elements fall within each block
    print("\nValidating blocks against hex elements...")
    print("=" * 80)
    
    block_stats = {}
    assigned_elements = set()
    
    for block_name, block_info in blocks.items():
        corners = np.array(block_info['physical'])
        
        # Count hex elements whose centers are inside this block
        inside_count = 0
        inside_indices = []
        
        for i, center in enumerate(hex_centers):
            if point_in_hex(center, corners):
                inside_count += 1
                inside_indices.append(i)
                assigned_elements.add(i)
        
        block_stats[block_name] = {
            'hex_count': inside_count,
            'hex_indices': inside_indices,
            'corner_bounds': {
                'x': [float(np.min(corners[:,0])), float(np.max(corners[:,0]))],
                'y': [float(np.min(corners[:,1])), float(np.max(corners[:,1]))],
                'z': [float(np.min(corners[:,2])), float(np.max(corners[:,2]))]
            }
        }
        
        print(f"\n{block_name}:")
        print(f"  Location: {block_info['location']}, Side: {block_info['side']}")
        print(f"  Bounds: x=[{np.min(corners[:,0]):.4f}, {np.max(corners[:,0]):.4f}], "
              f"y=[{np.min(corners[:,1]):.4f}, {np.max(corners[:,1]):.4f}], "
              f"z=[{np.min(corners[:,2]):.4f}, {np.max(corners[:,2]):.4f}]")
        print(f"  Hex elements inside: {inside_count} ({100*inside_count/len(hex_elements):.1f}%)")
    
    # Report unassigned elements
    unassigned = len(hex_elements) - len(assigned_elements)
    print(f"\n\nSummary:")
    print(f"  Total hex elements: {len(hex_elements)}")
    print(f"  Assigned to blocks: {len(assigned_elements)} ({100*len(assigned_elements)/len(hex_elements):.1f}%)")
    print(f"  Unassigned: {unassigned} ({100*unassigned/len(hex_elements):.1f}%)")
    
    if unassigned > 0:
        print(f"\n  WARNING: {unassigned} hex elements not in any block!")
        print(f"  This may indicate:")
        print(f"    - Blocks are defined for a single blade passage, but mesh has multiple passages")
        print(f"    - The mapping from parametric to physical space needs adjustment")
        print(f"    - Rounding/blade thickness not accounted for")
    
    # Save validation results
    output = {
        'metadata': {
            'total_hex_elements': len(hex_elements),
            'assigned_elements': len(assigned_elements),
            'unassigned_elements': unassigned
        },
        'block_stats': block_stats
    }
    
    with open('block_validation.json', 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"\nSaved validation results to block_validation.json")

if __name__ == '__main__':
    main()
