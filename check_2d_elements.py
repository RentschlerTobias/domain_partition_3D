#!/usr/bin/env python3
"""Check all 2D elements in the mesh regardless of physical tag."""
import numpy as np

MESH_PATH = '/root/repos/block_structured_meshing/T1_9/T1_9_ru_gridGmsh.msh'

# Parse nodes
nodes = {}
with open(MESH_PATH, 'r') as f:
    lines = f.readlines()
    in_nodes = False
    for line in lines:
        if '$Nodes' in line:
            in_nodes = True
            continue
        if '$EndNodes' in line:
            break
        if in_nodes:
            parts = line.strip().split()
            if len(parts) == 1:
                continue
            if len(parts) == 4:
                tag = int(parts[0])
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                nodes[tag] = (x, y, z)

print(f"Total nodes: {len(nodes)}")

# Parse all 2D elements (triangles and quads)
triangles = []
quads = []

with open(MESH_PATH, 'r') as f:
    lines = f.readlines()
    in_elements = False
    for line in lines:
        if '$Elements' in line:
            in_elements = True
            continue
        if '$EndElements' in line:
            break
        if in_elements:
            parts = line.strip().split()
            if len(parts) == 1:
                continue
            if len(parts) < 5:
                continue
            elem_type = int(parts[1])
            reg_phys = int(parts[2])
            reg_geom = int(parts[3])
            node_list = [int(p) for p in parts[4:]]
            
            if elem_type == 2:  # Triangle
                triangles.append((reg_phys, reg_geom, node_list))
            elif elem_type == 3:  # Quad
                quads.append((reg_phys, reg_geom, node_list))

print(f"Total triangles: {len(triangles)}")
print(f"Total quads: {len(quads)}")

# Group by physical tag
phys_tri_counts = {}
phys_quad_counts = {}
for reg_phys, reg_geom, nodes in triangles:
    phys_tri_counts[reg_phys] = phys_tri_counts.get(reg_phys, 0) + 1
for reg_phys, reg_geom, nodes in quads:
    phys_quad_counts[reg_phys] = phys_quad_counts.get(reg_phys, 0) + 1

print("\n=== Triangles per physical tag ===")
for phys, count in sorted(phys_tri_counts.items()):
    print(f"  {phys}: {count}")

print("\n=== Quads per physical tag ===")
for phys, count in sorted(phys_quad_counts.items()):
    print(f"  {phys}: {count}")

# Group by geometric tag
geom_tri_counts = {}
geom_quad_counts = {}
for reg_phys, reg_geom, nodes in triangles:
    geom_tri_counts[reg_geom] = geom_tri_counts.get(reg_geom, 0) + 1
for reg_phys, reg_geom, nodes in quads:
    geom_quad_counts[reg_geom] = geom_quad_counts.get(reg_geom, 0) + 1

print("\n=== Triangles per geometric tag ===")
for geom, count in sorted(geom_tri_counts.items()):
    print(f"  {geom}: {count}")

print("\n=== Quads per geometric tag ===")
for geom, count in sorted(geom_quad_counts.items()):
    print(f"  {geom}: {count}")

# Check unique geometric tags
print(f"\nUnique geometric tags with 2D elements: {sorted(set(geom_tri_counts.keys()) | set(geom_quad_counts.keys()))}")

# Check node ranges for each geometric tag
geom_nodes = {}
for reg_phys, reg_geom, node_list in triangles + quads:
    if reg_geom not in geom_nodes:
        geom_nodes[reg_geom] = set()
    geom_nodes[reg_geom].update(node_list)

print("\n=== Node count per geometric tag ===")
for geom, node_set in sorted(geom_nodes.items()):
    print(f"  {geom}: {len(node_set)} nodes")
    
    # Show bounding box
    coords = [nodes[n] for n in node_set if n in nodes]
    if coords:
        xs = [c[0] for c in coords]
        ys = [c[1] for c in coords]
        zs = [c[2] for c in coords]
        print(f"    x=[{min(xs):.4f}, {max(xs):.4f}], y=[{min(ys):.4f}, {max(ys):.4f}], z=[{min(zs):.4f}, {max(zs):.4f}]")
