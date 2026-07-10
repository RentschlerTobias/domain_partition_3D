#!/usr/bin/env python3
"""
Export actual mesh geometry (hex centers + blade surface) to JSON for Three.js visualization.
"""

import json
import numpy as np
from pathlib import Path

def parse_msh_file(filepath):
    """Parse hex elements and blade surface nodes from Gmsh 2.2 format."""
    nodes = {}
    hex_elements = []
    blade_tris = []
    blade_quads = []
    
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
                    n_tags = int(parts[2])
                    tags = [int(p) for p in parts[3:3+n_tags]]
                    node_ids = [int(p) for p in parts[3+n_tags:]]
                    
                    geom_tag = tags[1] if len(tags) > 1 else 0
                    
                    if elem_type == 5:  # Hex element
                        hex_elements.append(node_ids)
                    elif elem_type == 2 and geom_tag in [5, 6]:  # Triangles on blade
                        blade_tris.append(node_ids)
                    elif elem_type == 3 and geom_tag in [5, 6]:  # Quads on blade
                        blade_quads.append(node_ids)
            
            line = f.readline()
    
    return nodes, hex_elements, blade_tris, blade_quads

def main():
    print("Parsing T1_9 mesh file...")
    msh_path = Path('T1_9/T1_9_ru_gridGmsh.msh')
    nodes, hex_elements, blade_tris, blade_quads = parse_msh_file(msh_path)
    
    print(f"Found {len(hex_elements)} hex elements")
    print(f"Found {len(blade_tris)} blade triangles")
    print(f"Found {len(blade_quads)} blade quads")
    
    # Compute hex centers (subsample for performance)
    print("Computing hex centers...")
    hex_centers = []
    for elem_nodes in hex_elements[::20]:  # Every 20th for performance
        coords = [nodes[nid] for nid in elem_nodes]
        center = np.mean(coords, axis=0).tolist()
        hex_centers.append(center)
    
    # Extract blade surface points (unique nodes)
    print("Extracting blade surface...")
    blade_nodes = set()
    for tri in blade_tris:
        blade_nodes.update(tri)
    for quad in blade_quads:
        blade_nodes.update(quad)
    
    blade_coords = [nodes[nid] for nid in blade_nodes]
    
    # Save as JSON for Three.js
    output = {
        'metadata': {
            'n_hex_elements': len(hex_elements),
            'n_hex_centers_sampled': len(hex_centers),
            'n_blade_surface_nodes': len(blade_coords),
            'n_blade_tris': len(blade_tris),
            'n_blade_quads': len(blade_quads)
        },
        'hex_centers': hex_centers,
        'blade_surface': blade_coords
    }
    
    output_path = Path('mesh_geometry.json')
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"\nSaved to {output_path}")
    print(f"File size: {output_path.stat().st_size / 1024:.1f} KB")
    print("\nThis JSON can be loaded by the Three.js visualization.")

if __name__ == '__main__':
    main()
