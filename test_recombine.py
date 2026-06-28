import gmsh
import numpy as np

OUT_DIR = "/root/repos/block_structured_meshing"

def test_recombine(msh_file, name):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    
    gmsh.open(msh_file)
    
    # classify
    gmsh.model.mesh.classifySurfaces(np.pi/2, True, True, 1e-4)
    
    # check entities
    print(f"[{name}] entities: {len(gmsh.model.getEntities())}")
    print(f"[{name}] dim1: {len(gmsh.model.getEntities(1))}")
    
    # try to recombine
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    
    # try to mesh
    gmsh.model.mesh.generate(2)
    
    # info
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(2)
    n_quads = 0
    n_tris = 0
    for etype, etags, enodes in zip(elem_types, elem_tags, elem_node_tags):
        if etype == 3:
            n_quads = len(etags)
        elif etype == 2:
            n_tris = len(etags)
    
    print(f"[{name}] quads: {n_quads}, tris: {n_tris}")
    
    gmsh.finalize()

if __name__ == "__main__":
    test_recombine(f"{OUT_DIR}/hub_surface.msh", "hub")
