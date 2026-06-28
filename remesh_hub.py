#!/usr/bin/env python3
"""Remesh hub surface with quasi-structured quad algorithm."""
import sys
sys.path.insert(0, '/root/environments/blocking/lib/python3.12/site-packages')
import gmsh
import numpy as np

HUB_STL = '/root/repos/block_structured_meshing/T1_9_hub.stl'

gmsh.initialize()
gmsh.merge(HUB_STL)

# The hub is a flat surface. Classify and create geometry.
print("Classifying surfaces...")
gmsh.model.mesh.classifySurfaces(120 * np.pi / 180, True, False, 180 * np.pi / 180)
print("Creating geometry...")
gmsh.model.mesh.createGeometry()

# Get surfaces
surfaces = gmsh.model.getEntities(2)
print(f"Surfaces: {len(surfaces)}")
for dim, tag in surfaces:
    print(f"  Surface {tag}")

# Element size
field_id = gmsh.model.mesh.field.add("MathEval")
gmsh.model.mesh.field.setString(field_id, "F", "8")
gmsh.model.mesh.field.setAsBackgroundMesh(field_id)

# Set quasi-structured quad algorithm
gmsh.option.setNumber("Mesh.Algorithm", 11)

# Set mesh size
gmsh.option.setNumber("Mesh.MeshSizeMax", 1.0)
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.5)

# Generate mesh
print("Generating mesh...")
gmsh.model.mesh.generate(2)

# Count elements
node_tags, coords, _ = gmsh.model.mesh.getNodes()
print(f"Nodes: {len(node_tags)}")

for dim, tag in gmsh.model.getEntities(2):
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(2, tag)
    total = 0
    quads = 0
    for et, etags in zip(elem_types, elem_tags):
        total += len(etags)
        if et == 3:
            quads += len(etags)
    if total > 0:
        print(f"  Surface {tag}: {total} elements, {quads} quads")

# Write output
output_path = '/root/repos/block_structured_meshing/T1_9_hub_remeshed.msh'
gmsh.write(output_path)
print(f"\nWrote mesh to {output_path}")

gmsh.finalize()
print("Done.")
