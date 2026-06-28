#!/usr/bin/env python3
"""Analyze T1_9 hub surfaces - manual physical mapping."""
import sys
sys.path.insert(0, '/root/environments/blocking/lib/python3.12/site-packages')
import gmsh
import numpy as np

MESH_PATH = '/root/repos/block_structured_meshing/T1_9/T1_9_ru_gridGmsh.msh'

gmsh.initialize()
gmsh.open(MESH_PATH)

# Parse physical names from file
phys_map = {}
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
                        phys_map[tag] = name
                except:
                    pass

print("=== Physical name mapping ===")
for tag, name in phys_map.items():
    print(f"  {tag}: {name}")

# Try to get entities for each physical group
print("\n=== Entities per physical group ===")
for tag, name in phys_map.items():
    try:
        entities = gmsh.model.getEntitiesForPhysicalGroup(2, tag)
        print(f"  {tag} ({name}): entities={entities}")
    except Exception as e:
        print(f"  {tag} ({name}): ERROR {e}")

# Read all nodes
node_tags, coords, _ = gmsh.model.mesh.getNodes()
vertices = coords.reshape(-1, 3)
tag_to_idx = np.full(int(node_tags.max()) + 1, -1)
tag_to_idx[node_tags.astype(int)] = np.arange(len(node_tags))

# Get all 2D entities
entities = gmsh.model.getEntities(2)
print(f"\n=== Total 2D entities: {len(entities)} ===")

# Try to map surfaces by element type
# The mesh file has physical tags stored in the element tags
# Check element tags
for dim, tag in entities[:5]:
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(2, tag)
    for et, etags in zip(elem_types, elem_tags):
        print(f"  Surface {tag}: elem_type={et}, num_elements={len(etags)}, first_tag={etags[0] if len(etags)>0 else None}")

# Check if we can get physical tags from elements
print("\n=== Element physical tags ===")
for dim, tag in entities[:3]:
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(2, tag)
    for et, etags, entags in zip(elem_types, elem_tags, elem_node_tags):
        # Get physical tags for elements
        if len(etags) > 0:
            phys_tags = gmsh.model.getElementPhysicalTags(2, etags[0])
            print(f"  Surface {tag}: elem_type={et}, first_elem_tag={etags[0]}, phys_tags={phys_tags}")

gmsh.finalize()
