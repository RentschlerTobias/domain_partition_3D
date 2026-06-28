import gmsh
import numpy as np

OUT_DIR = "/root/repos/block_structured_meshing"

def test_classify(stl_file, name):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 11)
    
    print(f"[{name}] loading {stl_file}")
    gmsh.merge(stl_file)
    
    print(f"[{name}] classifySurfaces")
    gmsh.model.mesh.classifySurfaces(np.pi/2, True, True, 1e-4)
    
    # check what entities exist
    entities = gmsh.model.getEntities()
    print(f"[{name}] entities: {entities}")
    
    # get boundary
    for dim, tag in gmsh.model.getEntities(2):
        bnd = gmsh.model.getBoundary([(dim, tag)])
        print(f"[{name}] surface {tag} boundary: {bnd}")
    
    # try meshing without createGeometry
    print(f"[{name}] generate 2D mesh")
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
    
    print(f"[{name}] remeshed: {n_quads} quads, {n_tris} tris")
    
    gmsh.finalize()

if __name__ == "__main__":
    test_classify(f"{OUT_DIR}/T1_9_hub_raw.stl", "hub")
