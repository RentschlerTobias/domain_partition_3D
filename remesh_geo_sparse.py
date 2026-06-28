import gmsh
import numpy as np

OUT_DIR = "/root/repos/block_structured_meshing"

def remesh_with_geo_sparse(stl_file, name, out_msh, elem_size=0.5, skip=5):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 11)
    
    gmsh.merge(stl_file)
    
    # get nodes
    tags, coords, _ = gmsh.model.mesh.getNodes()
    coords = np.array(coords).reshape(-1, 3)
    
    # find boundary edges
    tri_nodes = []
    elem_types = gmsh.model.mesh.getElementTypes()
    for et in elem_types:
        if et == 2:  # 3-node triangle
            etags, enodes = gmsh.model.mesh.getElementsByType(et)
            tri_nodes = np.array(enodes).reshape(-1, 3)
            break
    
    # compute boundary edges
    edge_map = {}
    for tri in tri_nodes:
        for i in range(3):
            n1, n2 = sorted([tri[i], tri[(i+1)%3]])
            edge = (n1, n2)
            if edge in edge_map:
                edge_map[edge] += 1
            else:
                edge_map[edge] = 1
    
    boundary_edges = [e for e, count in edge_map.items() if count == 1]
    
    # order boundary edges into a loop
    loop = [boundary_edges[0]]
    used = {0}
    
    while len(loop) < len(boundary_edges):
        last = loop[-1]
        found = False
        for i, edge in enumerate(boundary_edges):
            if i in used:
                continue
            if edge[0] == last[1] or edge[1] == last[1]:
                if edge[0] == last[1]:
                    loop.append(edge)
                else:
                    loop.append((edge[1], edge[0]))
                used.add(i)
                found = True
                break
        if not found:
            break
    
    # extract boundary node tags
    boundary_nodes = [loop[0][0]]
    for edge in loop:
        boundary_nodes.append(edge[1])
    
    if boundary_nodes[-1] == boundary_nodes[0]:
        boundary_nodes = boundary_nodes[:-1]
    
    print(f"[{name}] boundary nodes: {len(boundary_nodes)}")
    
    # sparse subset
    sparse_nodes = [boundary_nodes[i] for i in range(0, len(boundary_nodes), skip)]
    if boundary_nodes[-1] not in sparse_nodes:
        sparse_nodes.append(boundary_nodes[-1])
    
    print(f"[{name}] sparse nodes: {len(sparse_nodes)}")
    
    # create new model with parametric geometry
    gmsh.clear()
    
    # add points
    point_tags = []
    for i, n in enumerate(sparse_nodes):
        idx = np.where(tags == n)[0][0]
        coord = coords[idx]
        tag = gmsh.model.geo.addPoint(coord[0], coord[1], coord[2], elem_size)
        point_tags.append(tag)
    
    # add lines
    line_tags = []
    for i in range(len(point_tags)):
        p1 = point_tags[i]
        p2 = point_tags[(i+1) % len(point_tags)]
        tag = gmsh.model.geo.addLine(p1, p2)
        line_tags.append(tag)
    
    # add curve loop
    loop_tag = gmsh.model.geo.addCurveLoop(line_tags)
    
    # add plane surface
    surf_tag = gmsh.model.geo.addPlaneSurface([loop_tag])
    
    gmsh.model.geo.synchronize()
    
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
    remesh_with_geo_sparse(f"{OUT_DIR}/T1_9_hub_raw.stl", "hub", f"{OUT_DIR}/T1_9_hub_quad.msh", elem_size=2.0, skip=5)
    remesh_with_geo_sparse(f"{OUT_DIR}/T1_9_shroud_raw.stl", "shroud", f"{OUT_DIR}/T1_9_shroud_quad.msh", elem_size=2.0, skip=5)
