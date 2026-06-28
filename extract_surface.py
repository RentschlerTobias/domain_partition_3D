import gmsh
import numpy as np

MSH_FILE = "/root/repos/block_structured_meshing/T1_9/T1_9_ru_gridGmsh.msh"
OUT_DIR = "/root/repos/block_structured_meshing"

def extract_surface_with_boundary(tag, name, out_msh):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    
    gmsh.open(MSH_FILE)
    
    # get boundary curves
    bnd = gmsh.model.getBoundary([(2, tag)])
    print(f"[{name}] boundary curves: {bnd}")
    
    # get all nodes and elements for surface and its boundary
    all_nodes = set()
    all_elem_types = []
    all_elem_tags = []
    all_elem_node_tags = []
    
    # get 2D elements
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(2, tag)
    all_elem_types.extend(elem_types)
    all_elem_tags.extend(elem_tags)
    all_elem_node_tags.extend(elem_node_tags)
    
    for enodes in elem_node_tags:
        for n in enodes:
            all_nodes.add(int(n))
    
    # get 1D elements on boundary curves
    for dim, bnd_tag in bnd:
        if dim == 1:
            elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(1, bnd_tag)
            all_elem_types.extend(elem_types)
            all_elem_tags.extend(elem_tags)
            all_elem_node_tags.extend(elem_node_tags)
            for enodes in elem_node_tags:
                for n in enodes:
                    all_nodes.add(int(n))
    
    # get 0D elements on boundary points
    for dim, bnd_tag in bnd:
        if dim == 1:
            curve_bnd = gmsh.model.getBoundary([(1, bnd_tag)])
            for p_dim, p_tag in curve_bnd:
                if p_dim == 0:
                    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(0, p_tag)
                    all_elem_types.extend(elem_types)
                    all_elem_tags.extend(elem_tags)
                    all_elem_node_tags.extend(elem_node_tags)
                    for enodes in elem_node_tags:
                        for n in enodes:
                            all_nodes.add(int(n))
    
    node_tags = sorted(list(all_nodes))
    
    # get coordinates
    coords = []
    for n in node_tags:
        result = gmsh.model.mesh.getNode(n)
        coord = result[0]
        coords.extend(coord)
    
    # create new model
    gmsh.clear()
    
    # add discrete surface
    surf = gmsh.model.addDiscreteEntity(2, 1)
    
    # add nodes
    gmsh.model.mesh.addNodes(2, 1, node_tags, coords)
    
    # add elements
    for etype, etags, enodes in zip(all_elem_types, all_elem_tags, all_elem_node_tags):
        gmsh.model.mesh.addElementsByType(1, etype, etags, enodes)
    
    gmsh.write(out_msh)
    print(f"[{name}] extracted {len(node_tags)} nodes, {sum(len(t) for t in all_elem_tags)} elements")
    
    gmsh.finalize()

if __name__ == "__main__":
    extract_surface_with_boundary(1, "hub", f"{OUT_DIR}/hub_surface.msh")
    extract_surface_with_boundary(2, "shroud", f"{OUT_DIR}/shroud_surface.msh")
