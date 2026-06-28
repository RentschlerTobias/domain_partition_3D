import gmsh
import numpy as np

MSH_FILE = "/root/repos/block_structured_meshing/T1_9/T1_9_ru_gridGmsh.msh"

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)

gmsh.open(MSH_FILE)

print("before classify")
print(f"entities: {len(gmsh.model.getEntities())}")
print(f"dim1: {len(gmsh.model.getEntities(1))}")

# classify
# gmsh.model.mesh.classifySurfaces(np.pi/2, True, True, 1e-4)

print("after classify")
print(f"entities: {len(gmsh.model.getEntities())}")
print(f"dim1: {len(gmsh.model.getEntities(1))}")

# create geometry from mesh
# gmsh.model.mesh.createGeometry()

# check if entities persist after clear
print("clearing mesh")
gmsh.model.mesh.clear()

print("after clear")
print(f"entities: {len(gmsh.model.getEntities())}")
print(f"dim1: {len(gmsh.model.getEntities(1))}")

for dim, tag in gmsh.model.getEntities(2)[:3]:
    bnd = gmsh.model.getBoundary([(dim, tag)])
    print(f"surface {tag} boundary: {len(bnd)} curves")

# try to mesh
print("trying to mesh surface 1")
gmsh.option.setNumber("Mesh.Algorithm", 11)
gmsh.model.mesh.generate(2)

print("mesh done")

# info
for dim, tag in gmsh.model.getEntities(2)[:3]:
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(2, tag)
    total = sum(len(t) for t in elem_tags)
    print(f"surface {tag}: {total} elements")

gmsh.finalize()
