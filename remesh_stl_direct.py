import gmsh
import numpy as np

OUT_DIR = "/root/repos/block_structured_meshing"

def remesh_stl_direct(stl_file, name, out_msh, elem_size=0.5):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 11)
    
    gmsh.merge(stl_file)
    
    print(f"[{name}] loaded STL")
    
    # classify surfaces (creates discrete parametrized surface)
    gmsh.model.mesh.classifySurfaces(np.pi/2, True, True, 1e-4)
    
    print(f"[{name}] classified surfaces")
    
    # create geometry from mesh
    gmsh.model.mesh.createGeometry()
    
    print(f"[{name}] created geometry")
    
    # set background field
    field = gmsh.model.mesh.field.add("MathEval")
    gmsh.model.mesh.field.setString(field, "F", str(elem_size))
    gmsh.model.mesh.field.setAsBackgroundMesh(field)
    
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", elem_size)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", elem_size)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromPoints", 0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromCurvature", 0)
    gmsh.option.setNumber("Mesh.QuadqsSizemapMethod", 0)
    
    # generate mesh
    print(f"[{name}] generating mesh...")
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
    
    gmsh.write(out_msh)
    print(f"[{name}] saved to {out_msh}")
    
    gmsh.finalize()

if __name__ == "__main__":
    remesh_stl_direct(f"{OUT_DIR}/T1_9_hub_raw.stl", "hub", f"{OUT_DIR}/T1_9_hub_quad.msh", elem_size=1.0)
    remesh_stl_direct(f"{OUT_DIR}/T1_9_shroud_raw.stl", "shroud", f"{OUT_DIR}/T1_9_shroud_quad.msh", elem_size=1.0)
