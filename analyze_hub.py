#!/usr/bin/env python3
"""Analyze T1_9 hub surfaces to decide merge vs split."""
import sys
sys.path.insert(0, '/root/environments/blocking/lib/python3.12/site-packages')
import gmsh
import numpy as np

MESH_PATH = '/root/repos/block_structured_meshing/T1_9/T1_9_ru_gridGmsh.msh'

gmsh.initialize()
gmsh.open(MESH_PATH)

# Get physical names
phys_groups = gmsh.model.getPhysicalGroups(2)
print("=== 2D Physical Surfaces ===")
for dim, tag in phys_groups:
    name = gmsh.model.getPhysicalName(2, tag)
    print(f"  {tag}: {name}")

# Get all 2D entities
entities = gmsh.model.getEntities(2)
print(f"\n=== Total 2D entities: {len(entities)} ===")

# Collect hub surfaces
hub_tags = []
shroud_tags = []
blade_tags = []
inlet_tags = []
outlet_tags = []
for dim, tag in entities:
    phys_tags = gmsh.model.getPhysicalGroupsForEntity(2, tag)
    for pt in phys_tags:
        name = gmsh.model.getPhysicalName(2, pt)
        if name and 'hub' in name.lower():
            hub_tags.append(tag)
        elif name and 'shroud' in name.lower():
            shroud_tags.append(tag)
        elif name and 'blade' in name.lower():
            blade_tags.append(tag)
        elif name and 'inlet' in name.lower():
            inlet_tags.append(tag)
        elif name and 'outlet' in name.lower():
            outlet_tags.append(tag)

print(f"\nHub surfaces: {hub_tags}")
print(f"Shroud surfaces: {shroud_tags}")
print(f"Blade surfaces: {blade_tags}")

# Nodes per hub surface
node_tags_all, coords_all, _ = gmsh.model.mesh.getNodes()
vertices = coords_all.reshape(-1, 3)
tag_to_idx = np.full(int(node_tags_all.max()) + 1, -1)
tag_to_idx[node_tags_all.astype(int)] = np.arange(len(node_tags_all))

print("\n=== Hub surface details ===")
for tag in hub_tags:
    _, q_tags = gmsh.model.mesh.getElementsByType(3, tag)
    n_quads = len(q_tags) // 4 if q_tags.size else 0
    nodes_in_surface = set()
    if q_tags.size:
        for qt in q_tags:
            nodes_in_surface.add(int(qt))
    print(f"  Surface {tag}: {n_quads} quads, {len(nodes_in_surface)} unique nodes")

# Check for duplicate nodes between hub surfaces
print("\n=== Shared nodes between hub surfaces ===")
for i, tag1 in enumerate(hub_tags):
    for j, tag2 in enumerate(hub_tags):
        if i >= j:
            continue
        _, q1 = gmsh.model.mesh.getElementsByType(3, tag1)
        _, q2 = gmsh.model.mesh.getElementsByType(3, tag2)
        nodes1 = set(int(q) for q in q1) if q1.size else set()
        nodes2 = set(int(q) for q in q2) if q2.size else set()
        shared = nodes1 & nodes2
        if shared:
            print(f"  {tag1} <-> {tag2}: {len(shared)} shared nodes")

# Check edge sharpness between patches
print("\n=== Edge classification between hub patches ===")
for i, tag1 in enumerate(hub_tags):
    for j, tag2 in enumerate(hub_tags):
        if i >= j:
            continue
        # Get boundary edges of each surface
        _, e1 = gmsh.model.mesh.getElementsByType(1, tag1)
        _, e2 = gmsh.model.mesh.getElementsByType(1, tag2)
        if e1.size and e2.size:
            edges1 = set((int(e1[k]), int(e1[k+1])) if int(e1[k]) < int(e1[k+1]) else (int(e1[k+1]), int(e1[k])) for k in range(0, len(e1), 2))
            edges2 = set((int(e2[k]), int(e2[k+1])) if int(e2[k]) < int(e2[k+1]) else (int(e2[k+1]), int(e2[k])) for k in range(0, len(e2), 2))
            shared_edges = edges1 & edges2
            if shared_edges:
                print(f"  {tag1} <-> {tag2}: {len(shared_edges)} shared edges (seam)")

# Check singularities on merged hub
print("\n=== Merged hub analysis ===")
all_hub_nodes = set()
all_hub_faces = []
for tag in hub_tags:
    _, q_tags = gmsh.model.mesh.getElementsByType(3, tag)
    if q_tags.size:
        for k in range(0, len(q_tags), 4):
            face = [int(q_tags[k]), int(q_tags[k+1]), int(q_tags[k+2]), int(q_tags[k+3])]
            all_hub_faces.append(face)
            all_hub_nodes.update(face)

print(f"Total hub faces: {len(all_hub_faces)}")
print(f"Total hub nodes: {len(all_hub_nodes)}")

# Compute valences
node_valence = {}
for face in all_hub_faces:
    for i in range(4):
        n1 = face[i]
        n2 = face[(i+1)%4]
        node_valence.setdefault(n1, set()).add(n2)
        node_valence.setdefault(n2, set()).add(n1)

singularities = []
for node, neighbors in node_valence.items():
    val = len(neighbors)
    if val != 4:
        # Check if on boundary
        # Count edges that appear only once
        edge_count = {}
        for f in all_hub_faces:
            for i in range(4):
                e = tuple(sorted([f[i], f[(i+1)%4]]))
                edge_count[e] = edge_count.get(e, 0) + 1
        boundary_nodes = set()
        for e, c in edge_count.items():
            if c == 1:
                boundary_nodes.update(e)
        if node not in boundary_nodes and val != 4:
            singularities.append((node, val))

print(f"Internal singularities: {len(singularities)}")
for node, val in singularities[:20]:
    print(f"  Node {node}: valence {val}")

# Check if singularities cluster at patch boundaries
if singularities and len(hub_tags) > 1:
    print("\n=== Singularity distribution by surface ===")
    for tag in hub_tags:
        _, q_tags = gmsh.model.mesh.getElementsByType(3, tag)
        nodes_in_surf = set(int(q) for q in q_tags) if q_tags.size else set()
        count = sum(1 for node, val in singularities if node in nodes_in_surf)
        print(f"  Surface {tag}: {count} singularities")

gmsh.finalize()
print("\nAnalysis complete.")
