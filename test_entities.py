import gmsh

OUT_DIR = "/root/repos/block_structured_meshing"

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)

gmsh.merge(f"{OUT_DIR}/T1_9_hub_raw.stl")

print("before classify")
print(f"entities: {gmsh.model.getEntities()}")

gmsh.model.mesh.classifySurfaces(3.14159/2, True, True, 1e-4)

print("after classify")
print(f"dim0: {len(gmsh.model.getEntities(0))}")
print(f"dim1: {len(gmsh.model.getEntities(1))}")
print(f"dim2: {len(gmsh.model.getEntities(2))}")

for d, t in gmsh.model.getEntities(1)[:10]:
    bnd = gmsh.model.getBoundary([(1, t)])
    print(f"curve {t}: boundary={bnd}")

gmsh.model.mesh.createGeometry()

print("after createGeometry")
print(f"dim0: {len(gmsh.model.getEntities(0))}")
print(f"dim1: {len(gmsh.model.getEntities(1))}")
print(f"dim2: {len(gmsh.model.getEntities(2))}")

for d, t in gmsh.model.getEntities(1)[:10]:
    bnd = gmsh.model.getBoundary([(1, t)])
    print(f"curve {t}: boundary={bnd}")

# get boundary of surface
for d, t in gmsh.model.getEntities(2):
    bnd = gmsh.model.getBoundary([(2, t)])
    print(f"surface {t} boundary: {len(bnd)} curves")
    for bd, bt in bnd[:10]:
        print(f"  {bd}/{bt}")

gmsh.finalize()
