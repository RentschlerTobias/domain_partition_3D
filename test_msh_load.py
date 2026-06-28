import gmsh

MSH_FILE = "/root/repos/block_structured_meshing/T1_9/T1_9_ru_gridGmsh.msh"

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)

print("opening MSH file")
gmsh.open(MSH_FILE)

# check entities
entities = gmsh.model.getEntities()
print(f"entities: {len(entities)}")
for dim in [0, 1, 2, 3]:
    ents = gmsh.model.getEntities(dim)
    print(f"dim {dim}: {len(ents)} entities")
    if ents:
        print(f"  tags: {[t for d,t in ents[:10]]}")

# check elements
for dim in [0, 1, 2, 3]:
    ents = gmsh.model.getEntities(dim)
    total_elems = 0
    for d, t in ents:
        elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(dim, t)
        for et, etags, en in zip(elem_types, elem_tags, elem_node_tags):
            total_elems += len(etags)
    print(f"dim {dim} total elements: {total_elems}")

# check physical groups
physical_groups = gmsh.model.getPhysicalGroups()
print(f"physical groups: {len(physical_groups)}")
for dim, tag in physical_groups[:10]:
    name = gmsh.model.getPhysicalName(dim, tag)
    ents = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
    print(f"  {dim} {tag} '{name}': {len(ents)} entities")

gmsh.finalize()
