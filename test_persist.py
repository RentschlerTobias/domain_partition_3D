import gmsh
import numpy as np

OUT_DIR = "/root/repos/block_structured_meshing"

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)

gmsh.merge(f"{OUT_DIR}/T1_9_hub_raw.stl")

print("before classify")
print(f"entities: {gmsh.model.getEntities()}")

gmsh.model.mesh.classifySurfaces(np.pi/2, True, True, 1e-4)

print("after classify")
print(f"entities: {gmsh.model.getEntities()}")

gmsh.model.mesh.createGeometry()

print("after createGeometry")
print(f"entities: {gmsh.model.getEntities()}")

print("clear mesh")
gmsh.model.mesh.clear()

print("after clear")
print(f"entities: {gmsh.model.getEntities()}")

# try to get boundary of surface
for dim, tag in gmsh.model.getEntities(2):
    bnd = gmsh.model.getBoundary([(dim, tag)])
    print(f"surface {tag} boundary: {bnd}")

gmsh.finalize()
