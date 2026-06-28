import gmsh

OUT_DIR = "/root/repos/block_structured_meshing"

for name in ["hub", "shroud"]:
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    
    msh_file = f"{OUT_DIR}/T1_9_{name}_blocks.msh"
    gmsh.open(msh_file)
    
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(2)
    n_quads = 0
    for etype, etags, enodes in zip(elem_types, elem_tags, elem_node_tags):
        if etype == 3:
            n_quads = len(etags)
    
    print(f"{name}_blocks: {n_quads} quads")
    
    gmsh.finalize()
