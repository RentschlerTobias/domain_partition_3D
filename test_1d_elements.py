import gmsh

MSH_FILE = "/root/repos/block_structured_meshing/T1_9/T1_9_ru_gridGmsh.msh"

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)

gmsh.open(MSH_FILE)

# get 1D elements
elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(1)
print("1D elements:")
for etype, etags, enodes in zip(elem_types, elem_tags, elem_node_tags):
    if etype == 1:  # 2-node line
        print(f"  {len(etags)} lines")
        # get parent entities
        for i in range(min(10, len(etags))):
            tag = etags[i]
            nodes = enodes[i*2:(i+1)*2]
            # get parent entity
            elem_type, elem_node_tags, elem_dim, elem_tag = gmsh.model.mesh.getElement(tag)
            parent_tag = elem_tag
            print(f"    tag={tag}, nodes={list(nodes)}, parent={parent_tag}")

# get entities
print("\nEntities:")
for dim in [0, 1, 2, 3]:
    ents = gmsh.model.getEntities(dim)
    print(f"  dim {dim}: {len(ents)} entities")
    if dim == 1 and ents:
        for d, t in ents[:10]:
            print(f"    tag={t}")

# get boundary of surface 1
print("\nBoundary of surface 1:")
bnd = gmsh.model.getBoundary([(2, 1)])
print(f"  {bnd}")

# get elements for each boundary curve
for dim, tag in bnd:
    if dim == 1:
        elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(1, tag)
        total = sum(len(t) for t in elem_tags)
        print(f"  curve {tag}: {total} elements")

gmsh.finalize()
