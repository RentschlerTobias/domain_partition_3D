#!/usr/bin/env python3
"""Test STL remeshing with Gmsh quasi-structured quad algorithm."""
import sys
sys.path.insert(0, '/root/environments/blocking/lib/python3.12/site-packages')
import gmsh
import numpy as np

STL_PATH = '/root/repos/block_structured_meshing/T1_9/T1_9_ru_gridGmsh.stl'

print("Initializing Gmsh...")
gmsh.initialize()
print(f"Merging STL: {STL_PATH}")
gmsh.merge(STL_PATH)

# Check initial entities
print("\n=== Initial entities ===")
for dim in [0, 1, 2, 3]:
    entities = gmsh.model.getEntities(dim)
    print(f"  {dim}D: {len(entities)} entities")

# Classify surfaces
print("\nClassifying surfaces...")
angle = 60
include_boundary = 1
force_param = 0
curve_angle = 180
gmsh.model.mesh.classifySurfaces(angle * np.pi / 180, include_boundary, force_param, curve_angle * np.pi / 180)

# Create geometry
print("Creating geometry...")
gmsh.model.mesh.createGeometry()

# Check classified surfaces
print("\n=== Classified 2D surfaces ===")
surfaces = gmsh.model.getEntities(2)
print(f"Total 2D surfaces: {len(surfaces)}")
for dim, tag in surfaces:
    print(f"  Surface {tag}")

# Get nodes
node_tags, coords, _ = gmsh.model.mesh.getNodes()
print(f"\nTotal nodes: {len(node_tags)}")

# Element size
field_id = gmsh.model.mesh.field.add("MathEval")
gmsh.model.mesh.field.setString(field_id, "F", "4")
gmsh.model.mesh.field.setAsBackgroundMesh(field_id)

# Set quasi-structured quad algorithm
print("\nSetting mesh algorithm to 11 (quasi-structured quad)...")
gmsh.option.setNumber("Mesh.Algorithm", 11)

# Generate mesh
print("Generating mesh...")
gmsh.model.mesh.generate(2)

# Check mesh
print("\n=== Mesh statistics ===")
node_tags, coords, _ = gmsh.model.mesh.getNodes()
print(f"Nodes after meshing: {len(node_tags)}")

# Count elements
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
output_path = '/root/repos/block_structured_meshing/T1_9_remeshed.msh'
gmsh.write(output_path)
print(f"\nWrote mesh to {output_path}")

gmsh.finalize()
print("Done.")
