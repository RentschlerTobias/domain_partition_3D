import gmsh
import sys

MSH_FILE = "/root/repos/block_structured_meshing/T1_9/T1_9_ru_gridGmsh.msh"
OUT_DIR = "/root/repos/block_structured_meshing"

def remesh_surface(tag, name, out_msh, elem_size=0.5):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 11)
    
    print(f"[{name}] loading MSH")
    gmsh.open(MSH_FILE)
    
    # clear mesh
    print(f"[{name}] clearing mesh")
    gmsh.model.mesh.clear()
    
    # set background field
    print(f"[{name}] setting background field")
    field = gmsh.model.mesh.field.add("MathEval")
    gmsh.model.mesh.field.setString(field, "F", str(elem_size))
    gmsh.model.mesh.field.setAsBackgroundMesh(field)
    
    # generate mesh for specific surface
    print(f"[{name}] generating 2D mesh for surface {tag}")
    gmsh.model.mesh.generate(2)
    
    # info
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(2, tag)
    n_quads = 0
    n_tris = 0
    for etype, etags, enodes in zip(elem_types, elem_tags, elem_node_tags):
        if etype == 3:
            n_quads = len(etags)
        elif etype == 2:
            n_tris = len(etags)
    
    print(f"[{name}] remeshed: {n_quads} quads, {n_tris} tris")
    
    # save
    gmsh.write(out_msh)
    print(f"[{name}] saved to {out_msh}")
    
    gmsh.finalize()

if __name__ == "__main__":
    remesh_surface(1, "hub", f"{OUT_DIR}/T1_9_hub_quad.msh", elem_size=0.5)
    remesh_surface(2, "shroud", f"{OUT_DIR}/T1_9_shroud_quad.msh", elem_size=0.5)
