#!/usr/bin/env python3
"""Analyze T1_9 hub surfaces using direct mesh reading."""
import sys
sys.path.insert(0, '/root/environments/blocking/lib/python3.12/site-packages')
import gmsh
import numpy as np

MESH_PATH = '/root/repos/block_structured_meshing/T1_9/T1_9_ru_gridGmsh.msh'

gmsh.initialize()
gmsh.open(MESH_PATH)

# Read physical names from $PhysicalNames section
print("=== All Physical Groups ===")
for dim in [0, 1, 2, 3]:
    groups = gmsh.model.getPhysicalGroups(dim)
    for g in groups:
        tag = g[1]
        name = gmsh.model.getPhysicalName(dim, tag)
        print(f"  dim={dim} tag={tag}: {name}")

# Get all 2D entities
entities = gmsh.model.getEntities(2)
print(f"\n=== Total 2D entities: {len(entities)} ===")

# List all entities with their names
for dim, tag in entities:
    phys_tags = gmsh.model.getPhysicalGroupsForEntity(2, tag)
    names = [gmsh.model.getPhysicalName(2, pt) for pt in phys_tags]
    print(f"  Surface {tag}: physical={phys_tags}, names={names}")

# Get all 3D entities
vol_entities = gmsh.model.getEntities(3)
print(f"\n=== Total 3D entities: {len(vol_entities)} ===")
for dim, tag in vol_entities:
    phys_tags = gmsh.model.getPhysicalGroupsForEntity(3, tag)
    names = [gmsh.model.getPhysicalName(3, pt) for pt in phys_tags]
    print(f"  Volume {tag}: physical={phys_tags}, names={names}")

# Element counts per entity
print("\n=== Element counts per 2D entity ===")
for dim, tag in entities:
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(2, tag)
    total = 0
    for et, etags in zip(elem_types, elem_tags):
        total += len(etags)
    if total > 0:
        print(f"  Surface {tag}: {total} elements")

# Read physical names from file directly
print("\n=== Physical names from file ===")
with open(MESH_PATH, 'r') as f:
    lines = f.readlines()
    in_phys = False
    count = 0
    for line in lines:
        if '$PhysicalNames' in line:
            in_phys = True
            continue
        if '$EndPhysicalNames' in line:
            break
        if in_phys:
            parts = line.strip().split()
            if len(parts) == 1:
                count = int(parts[0])
                continue
            if len(parts) >= 3:
                dim = int(parts[0])
                tag = int(parts[1])
                name = parts[2].strip('"')
                if dim == 2:
                    print(f"  2D {tag}: {name}")

gmsh.finalize()
