#!/usr/bin/env python3
"""Parse T1_9 mesh manually to group elements by physical tag."""
import numpy as np

MESH_PATH = '/root/repos/block_structured_meshing/T1_9/T1_9_ru_gridGmsh.msh'

# Parse physical names
phys_names = {}
with open(MESH_PATH, 'r') as f:
    lines = f.readlines()
    in_phys = False
    for line in lines:
        if '$PhysicalNames' in line:
            in_phys = True
            continue
        if '$EndPhysicalNames' in line:
            break
        if in_phys:
            parts = line.strip().split()
            if len(parts) >= 3:
                try:
                    dim = int(parts[0])
                    tag = int(parts[1])
                    name = parts[2].strip('"')
                    if dim == 2:
                        phys_names[tag] = name
                except:
                    pass

# Parse elements
# Group by physical tag
phys_elements = {}  # phys_tag -> list of (elem_type, nodes)

with open(MESH_PATH, 'r') as f:
    lines = f.readlines()
    in_elements = False
    num_elements = 0
    count = 0
    for line in lines:
        if '$Elements' in line:
            in_elements = True
            continue
        if '$EndElements' in line:
            break
        if in_elements:
            parts = line.strip().split()
            if len(parts) == 1 and num_elements == 0:
                num_elements = int(parts[0])
                continue
            if len(parts) < 3:
                continue
            # Format: elm-number elm-type reg-phys reg-geom node1 node2 ...
            elem_num = int(parts[0])
            elem_type = int(parts[1])
            reg_phys = int(parts[2])
            reg_geom = int(parts[3])
            nodes = [int(p) for p in parts[4:]]
            
            if elem_type == 3:  # Quadrilateral
                if reg_phys not in phys_elements:
                    phys_elements[reg_phys] = []
                phys_elements[reg_phys].append(nodes)

print("=== Quadrilateral elements per physical group ===")
for phys_tag, elements in sorted(phys_elements.items()):
    name = phys_names.get(phys_tag, f"unknown_{phys_tag}")
    print(f"  {phys_tag} ({name}): {len(elements)} quads")

# Group by type
hub_phys = []
shroud_phys = []
blade_phys = []
inlet_phys = []
outlet_phys = []
for phys_tag, name in phys_names.items():
    if 'hub' in name.lower():
        hub_phys.append(phys_tag)
    elif 'shroud' in name.lower():
        shroud_phys.append(phys_tag)
    elif 'blade' in name.lower():
        blade_phys.append(phys_tag)
    elif 'inlet' in name.lower():
        inlet_phys.append(phys_tag)
    elif 'outlet' in name.lower():
        outlet_phys.append(phys_tag)

print(f"\nHub phys tags: {hub_phys}")
print(f"Shroud phys tags: {shroud_phys}")
print(f"Blade phys tags: {blade_phys}")

# Collect all hub quads
hub_faces = []
for phys_tag in hub_phys:
    if phys_tag in phys_elements:
        hub_faces.extend(phys_elements[phys_tag])

# Collect all shroud quads
shroud_faces = []
for phys_tag in shroud_phys:
    if phys_tag in phys_elements:
        shroud_faces.extend(phys_elements[phys_tag])

print(f"\nTotal hub quads: {len(hub_faces)}")
print(f"Total shroud quads: {len(shroud_faces)}")

# Check node sharing between hub patches
print("\n=== Shared nodes between hub patches ===")
for i, phys1 in enumerate(hub_phys):
    for j, phys2 in enumerate(hub_phys):
        if i >= j:
            continue
        if phys1 in phys_elements and phys2 in phys_elements:
            nodes1 = set()
            for face in phys_elements[phys1]:
                nodes1.update(face)
            nodes2 = set()
            for face in phys_elements[phys2]:
                nodes2.update(face)
            shared = nodes1 & nodes2
            if shared:
                print(f"  {phys1} ({phys_names.get(phys1)}) <-> {phys2} ({phys_names.get(phys2)}): {len(shared)} shared nodes")

# Check singularities on merged hub
print("\n=== Merged hub singularity analysis ===")
all_nodes = set()
for face in hub_faces:
    all_nodes.update(face)

node_valence = {}
for face in hub_faces:
    for i in range(4):
        n1 = face[i]
        n2 = face[(i+1)%4]
        node_valence.setdefault(n1, set()).add(n2)
        node_valence.setdefault(n2, set()).add(n1)

# Find boundary edges
edge_count = {}
for face in hub_faces:
    for i in range(4):
        e = tuple(sorted([face[i], face[(i+1)%4]]))
        edge_count[e] = edge_count.get(e, 0) + 1

boundary_nodes = set()
for e, c in edge_count.items():
    if c == 1:
        boundary_nodes.update(e)

singularities = []
for node, neighbors in node_valence.items():
    val = len(neighbors)
    if node not in boundary_nodes and val != 4:
        singularities.append((node, val))

print(f"Internal singularities: {len(singularities)}")
for node, val in singularities[:20]:
    print(f"  Node {node}: valence {val}")

# Check distribution by original patch
print("\n=== Singularity distribution by patch ===")
for phys_tag in hub_phys:
    if phys_tag in phys_elements:
        nodes_in_patch = set()
        for face in phys_elements[phys_tag]:
            nodes_in_patch.update(face)
        count = sum(1 for node, val in singularities if node in nodes_in_patch)
        print(f"  {phys_tag} ({phys_names.get(phys_tag)}): {count} singularities")

print("\nAnalysis complete.")
