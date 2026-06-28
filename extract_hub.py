#!/usr/bin/env python3
"""Extract hub surface (z=0) from volume mesh and check for quads."""
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
            if len(parts) == 4:
                tag = int(parts[0])
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                nodes[tag] = (x, y, z)

# Parse all 2D elements
faces = []
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
            node_list = [int(p) for p in parts[4:]]
            
            if elem_type == 2:  # Triangle
                faces.append(('tri', node_list))
            elif elem_type == 3:  # Quad
                faces.append(('quad', node_list))

# Extract hub faces (all nodes at z=0)
hub_faces = []
for ftype, node_list in faces:
    all_z0 = all(abs(nodes[n][2]) < 1e-6 for n in node_list if n in nodes)
    if all_z0:
        hub_faces.append((ftype, node_list))

print(f"Hub faces: {len(hub_faces)}")
tri_count = sum(1 for ftype, _ in hub_faces if ftype == 'tri')
quad_count = sum(1 for ftype, _ in hub_faces if ftype == 'quad')
print(f"  Triangles: {tri_count}")
print(f"  Quads: {quad_count}")

# Extract shroud faces (all nodes at z=2.5)
shroud_faces = []
for ftype, node_list in faces:
    all_z25 = all(abs(nodes[n][2] - 2.5) < 1e-6 for n in node_list if n in nodes)
    if all_z25:
        shroud_faces.append((ftype, node_list))

print(f"\nShroud faces: {len(shroud_faces)}")
tri_count = sum(1 for ftype, _ in shroud_faces if ftype == 'tri')
quad_count = sum(1 for ftype, _ in shroud_faces if ftype == 'quad')
print(f"  Triangles: {tri_count}")
print(f"  Quads: {quad_count}")

# Extract blade faces (all nodes NOT at z=0 or z=2.5)
blade_faces = []
for ftype, node_list in faces:
    zs = [nodes[n][2] for n in node_list if n in nodes]
    if all(z > 1e-6 and z < 2.5 - 1e-6 for z in zs):
        blade_faces.append((ftype, node_list))

print(f"\nBlade faces: {len(blade_faces)}")
tri_count = sum(1 for ftype, _ in blade_faces if ftype == 'tri')
quad_count = sum(1 for ftype, _ in blade_faces if ftype == 'quad')
print(f"  Triangles: {tri_count}")
print(f"  Quads: {quad_count}")

# Check if hub is a single surface
hub_nodes = set()
for ftype, node_list in hub_faces:
    hub_nodes.update(node_list)
print(f"\nHub unique nodes: {len(hub_nodes)}")

# Check hub connectivity
hub_node_neighbors = {}
for ftype, node_list in hub_faces:
    for i in range(len(node_list)):
        n1 = node_list[i]
        n2 = node_list[(i+1) % len(node_list)]
        hub_node_neighbors.setdefault(n1, set()).add(n2)
        hub_node_neighbors.setdefault(n2, set()).add(n1)

# Check valences
valences = {n: len(neighbors) for n, neighbors in hub_node_neighbors.items()}
val4 = sum(1 for v in valences.values() if v == 4)
val_not4 = sum(1 for v in valences.values() if v != 4)
print(f"Hub nodes with valence 4: {val4}")
print(f"Hub nodes with valence != 4: {val_not4}")

# Check if hub is a quad mesh
if quad_count > 0 and tri_count == 0:
    print("\nHub is a pure quad mesh!")
else:
    print(f"\nHub is NOT a pure quad mesh. Needs remeshing.")

# Check if shroud is a quad mesh
if quad_count > 0 and tri_count == 0:
    print("Shroud is a pure quad mesh!")
else:
    print(f"Shroud is NOT a pure quad mesh. Needs remeshing.")
