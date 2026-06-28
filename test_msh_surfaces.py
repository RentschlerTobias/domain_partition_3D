import gmsh

MSH_FILE = "/root/repos/block_structured_meshing/T1_9/T1_9_ru_gridGmsh.msh"

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)

gmsh.open(MSH_FILE)

# identify surfaces by element count
for dim, tag in gmsh.model.getEntities(2):
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(2, tag)
    total = sum(len(t) for t in elem_tags)
    print(f"surface {tag}: {total} elements")

# identify hub and shroud by physical names (if any)
# check if physical names exist
try:
    for dim, tag in gmsh.model.getPhysicalGroups():
        name = gmsh.model.getPhysicalName(dim, tag)
        ents = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
        print(f"physical {dim} {tag} '{name}': {ents}")
except:
    print("no physical groups")

gmsh.finalize()
