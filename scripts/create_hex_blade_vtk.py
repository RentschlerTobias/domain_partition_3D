#!/usr/bin/env python3
"""
Create detailed VTK visualization of hex mesh blocks around blade.
Show each region with different colors and add blade surface.
"""

import numpy as np
from collections import defaultdict

def parse_msh(filename):
    nodes = {}
    hex_elements = []
    blade_faces = []  # 2D faces on blade surface
    
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
            nodes[tag] = np.array([x, y, z])
    
    # Parse elements
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
            tag = int(parts[0])
            num_tags = int(parts[2])
            phys_tag = int(parts[3]) if num_tags >= 1 else 0
            geom_tag = int(parts[4]) if num_tags >= 2 else 0
            
            if elem_type == 5:  # Hex
                nodes_idx = [int(p) for p in parts[3 + num_tags:]]
                hex_elements.append({
                    'tag': tag,
                    'phys_tag': phys_tag,
                    'geom_tag': geom_tag,
                    'nodes': nodes_idx
                })
            elif elem_type == 3:  # Quad - potential blade face
                nodes_idx = [int(p) for p in parts[3 + num_tags:]]
                blade_faces.append({
                    'tag': tag,
                    'phys_tag': phys_tag,
                    'geom_tag': geom_tag,
                    'nodes': nodes_idx
                })
    
    return nodes, hex_elements, blade_faces

def write_combined_vtk(nodes, hex_elements, blade_faces, filename):
    """Write VTK with hex blocks and blade surface."""
    # Collect all nodes
    used_node_tags = set()
    for elem in hex_elements:
        used_node_tags.update(elem['nodes'])
    for face in blade_faces:
        used_node_tags.update(face['nodes'])
    
    used_node_tags = sorted(used_node_tags)
    node_index_map = {tag: idx for idx, tag in enumerate(used_node_tags)}
    
    with open(filename, 'w') as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("T1_9 Hex Mesh with Blade\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        
        # Points
        f.write(f"POINTS {len(used_node_tags)} float\n")
        for tag in used_node_tags:
            x, y, z = nodes[tag]
            f.write(f"{x} {y} {z}\n")
        
        # Cells: hex elements + blade quads
        total_cells = len(hex_elements) + len(blade_faces)
        total_ints = len(hex_elements) * 9 + len(blade_faces) * 5
        
        f.write(f"\nCELLS {total_cells} {total_ints}\n")
        
        # Hex elements
        for elem in hex_elements:
            f.write("8 ")
            for nid in elem['nodes']:
                f.write(f"{node_index_map[nid]} ")
            f.write("\n")
        
        # Blade faces (quads)
        for face in blade_faces:
            f.write("4 ")
            for nid in face['nodes']:
                f.write(f"{node_index_map[nid]} ")
            f.write("\n")
        
        # Cell types
        f.write(f"\nCELL_TYPES {total_cells}\n")
        for _ in hex_elements:
            f.write("12\n")  # VTK_HEXAHEDRON
        for _ in blade_faces:
            f.write("9\n")  # VTK_QUAD
        
        # Cell data
        f.write(f"\nCELL_DATA {total_cells}\n")
        
        # Region ID for hex, -1 for blade
        f.write("SCALARS region_id int 1\n")
        f.write("LOOKUP_TABLE default\n")
        for elem in hex_elements:
            f.write(f"{elem['geom_tag']}\n")
        for _ in blade_faces:
            f.write("-1\n")
        
        # Cell type: 0=hex, 1=blade
        f.write("SCALARS cell_type int 1\n")
        f.write("LOOKUP_TABLE default\n")
        for _ in hex_elements:
            f.write("0\n")
        for _ in blade_faces:
            f.write("1\n")
    
    print(f"Combined VTK written: {filename}")
    print(f"  Hex elements: {len(hex_elements)}")
    print(f"  Blade faces: {len(blade_faces)}")
    print(f"  Total nodes: {len(used_node_tags)}")

def write_block_boundary_vtk(nodes, hex_elements, filename):
    """Write only boundary faces of hex blocks as VTK."""
    # Find all boundary faces (faces not shared between two hex elements)
    face_to_elems = defaultdict(list)
    
    for elem in hex_elements:
        n = elem['nodes']
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
    
    # Boundary faces are those shared by only 1 element
    boundary_faces = []
    for face, elems in face_to_elems.items():
        if len(elems) == 1:
            boundary_faces.append({
                'nodes': list(face),
                'elem_tag': elems[0]
            })
    
    # Group by region
    elem_map = {e['tag']: e for e in hex_elements}
    
    used_node_tags = set()
    for face in boundary_faces:
        used_node_tags.update(face['nodes'])
    
    used_node_tags = sorted(used_node_tags)
    node_index_map = {tag: idx for idx, tag in enumerate(used_node_tags)}
    
    with open(filename, 'w') as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("T1_9 Hex Block Boundaries\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        
        f.write(f"POINTS {len(used_node_tags)} float\n")
        for tag in used_node_tags:
            x, y, z = nodes[tag]
            f.write(f"{x} {y} {z}\n")
        
        f.write(f"\nCELLS {len(boundary_faces)} {len(boundary_faces) * 5}\n")
        for face in boundary_faces:
            f.write("4 ")
            for nid in face['nodes']:
                f.write(f"{node_index_map[nid]} ")
            f.write("\n")
        
        f.write(f"\nCELL_TYPES {len(boundary_faces)}\n")
        for _ in boundary_faces:
            f.write("9\n")
        
        # Region ID
        f.write(f"\nCELL_DATA {len(boundary_faces)}\n")
        f.write("SCALARS region_id int 1\n")
        f.write("LOOKUP_TABLE default\n")
        for face in boundary_faces:
            elem = elem_map[face['elem_tag']]
            f.write(f"{elem['geom_tag']}\n")
    
    print(f"Boundary VTK written: {filename}")
    print(f"  Boundary faces: {len(boundary_faces)}")

def main():
    msh_file = 'T1_9/T1_9_ru_gridGmsh.msh'
    
    print("Parsing mesh...")
    nodes, hex_elements, blade_faces = parse_msh(msh_file)
    print(f"Loaded {len(nodes)} nodes, {len(hex_elements)} hex elements, {len(blade_faces)} blade faces")
    
    # Filter blade faces (geom tags 5 and 6 from previous analysis)
    blade_faces_filtered = [f for f in blade_faces if f['geom_tag'] in [5, 6]]
    print(f"Blade faces (geom 5,6): {len(blade_faces_filtered)}")
    
    # Write combined visualization
    write_combined_vtk(nodes, hex_elements, blade_faces_filtered, 'T1_9_hex_and_blade.vtk')
    
    # Write block boundaries only
    write_block_boundary_vtk(nodes, hex_elements, 'T1_9_hex_boundaries.vtk')
    
    # Write individual region VTKs with blade
    regions = defaultdict(list)
    for elem in hex_elements:
        regions[elem['geom_tag']].append(elem)
    
    for geom_tag in sorted(regions.keys()):
        region_elems = regions[geom_tag]
        filename = f'T1_9_region_{geom_tag}_with_blade.vtk'
        write_combined_vtk(nodes, region_elems, blade_faces_filtered, filename)

if __name__ == '__main__':
    main()
