#!/usr/bin/env python3
"""Export T1_9 geometry from dtOO docker container."""
import sys
sys.path.insert(0, '/dtOO-install/tools')
import dtOOPythonSWIG as dtOO
import gmsh

# Initialize
dtOO.dtDefaults.init()

# Build the T1_9 geometry
bC, cV, aF, aG, bV, dC, dP = dtOO.tistos.build()

# Get the Gmsh model
gmsh.model.set_current('ru_gridGmsh')

# Save the geometry
output_path = '/dtOO/test/simpleAxialRunner/T1_9_geometry.msh'
gmsh.write(output_path)
print(f"Geometry saved to {output_path}")

# Print surface info
surfaces = gmsh.model.getEntities(2)
print(f"Surfaces: {len(surfaces)}")
for dim, tag in surfaces:
    name = gmsh.model.getEntityName(dim, tag)
    print(f"  Surface {tag}: {name}")
