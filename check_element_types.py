#!/usr/bin/env python3
"""Check all element types for hub and shroud surfaces."""
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

# Parse elements by physical tag
phys_elements = {}  # phys_tag -> list of (elem_type, nodes)

with open(MESH_PATH, 'r') as f:
    lines = f.readlines()
    in_elements = False
    num_elements = 0
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
            elem_type = int(parts[1])
            reg_phys = int(parts[2])
            nodes = [int(p) for p in parts[4:]]
            
            if reg_phys not in phys_elements:
                phys_elements[reg_phys] = []
            phys_elements[reg_phys].append((elem_type, nodes))

# Element type names
etype_names = {
    1: '2-node line',
    2: '3-node triangle',
    3: '4-node quad',
    4: '4-node tetrahedron',
    5: '8-node hexahedron',
    6: '6-node prism',
    7: '5-node pyramid',
    8: '3-node line2',
    9: '6-node triangle2',
    10: '9-node quad2',
    11: '10-node tetrahedron2',
    12: '27-node hexahedron2',
    13: '18-node prism2',
    14: '14-node pyramid2',
    15: '1-node point',
    16: '8-node quad2',
    17: '20-node hexahedron2',
    18: '15-node prism2',
    19: '13-node pyramid2',
}

# Count element types per physical group
print("=== Element types per physical group ===")
for phys_tag in sorted(phys_elements.keys()):
    name = phys_names.get(phys_tag, f"unknown_{phys_tag}")
    type_counts = {}
    for elem_type, nodes in phys_elements[phys_tag]:
        type_counts[elem_type] = type_counts.get(elem_type, 0) + 1
    
    types_str = ", ".join([f"{etype_names.get(t, f'type_{t}')}: {c}" for t, c in sorted(type_counts.items())])
    print(f"  {phys_tag} ({name}): {types_str}")

# Hub and shroud
hub_phys = [6, 7, 8, 9, 10, 11]
shroud_phys = [12, 13, 14, 15, 16, 17]

print("\n=== Hub surface element types ===")
for phys_tag in hub_phys:
    if phys_tag in phys_elements:
        name = phys_names.get(phys_tag, "unknown")
        type_counts = {}
        for elem_type, nodes in phys_elements[phys_tag]:
            type_counts[elem_type] = type_counts.get(elem_type, 0) + 1
        types_str = ", ".join([f"{etype_names.get(t, f'type_{t}')}: {c}" for t, c in sorted(type_counts.items())])
        print(f"  {phys_tag} ({name}): {types_str}")

print("\n=== Shroud surface element types ===")
for phys_tag in shroud_phys:
    if phys_tag in phys_elements:
        name = phys_names.get(phys_tag, "unknown")
        type_counts = {}
        for elem_type, nodes in phys_elements[phys_tag]:
            type_counts[elem_type] = type_counts.get(elem_type, 0) + 1
        types_str = ", ".join([f"{etype_names.get(t, f'type_{t}')}: {c}" for t, c in sorted(type_counts.items())])
        print(f"  {phys_tag} ({name}): {types_str}")

# Total counts
hub_total = 0
for phys_tag in hub_phys:
    if phys_tag in phys_elements:
        hub_total += len(phys_elements[phys_tag])

shroud_total = 0
for phys_tag in shroud_phys:
    if phys_tag in phys_elements:
        shroud_total += len(phys_elements[phys_tag])

print(f"\nTotal hub elements: {hub_total}")
print(f"Total shroud elements: {shroud_total}")

# Check if hub has quads
hub_quads = 0
for phys_tag in hub_phys:
    if phys_tag in phys_elements:
        for elem_type, nodes in phys_elements[phys_tag]:
            if elem_type == 3:
                hub_quads += 1

print(f"Total hub quads: {hub_quads}")

# Check if hub has triangles
hub_triangles = 0
for phys_tag in hub_phys:
    if phys_tag in phys_elements:
        for elem_type, nodes in phys_elements[phys_tag]:
            if elem_type == 2:
                hub_triangles += 1

print(f"Total hub triangles: {hub_triangles}")
