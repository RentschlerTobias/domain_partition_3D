import gmsh
import math


def mesh_stl(path_to_stl, output_path="output.msh", element_size=4,
             angle=40, force_parametrizable=False, include_boundary=True,
             curve_angle=180):

    gmsh.initialize()

    # Merge STL
    gmsh.merge(path_to_stl)

    # ClassifySurfaces
    gmsh.model.mesh.classifySurfaces(
        angle * math.pi / 180,
        include_boundary,
        force_parametrizable,
        curve_angle * math.pi / 180
    )

    # CreateGeometry
    gmsh.model.mesh.createGeometry()

    # Surface Loop und Volume
    surfaces = gmsh.model.getEntities(2)
    surface_tags = [s[1] for s in surfaces]
    surface_loop = gmsh.model.geo.addSurfaceLoop(surface_tags)
    gmsh.model.geo.addVolume([surface_loop])
    gmsh.model.geo.synchronize()

    # Size field
    f = gmsh.model.mesh.field.add("MathEval")
    gmsh.model.mesh.field.setString(f, "F", str(element_size))
    gmsh.model.mesh.field.setAsBackgroundMesh(f)

    # Quasi-structured Quad Algorithmus
    gmsh.option.setNumber("Mesh.Algorithm", 11)

    gmsh.model.mesh.generate(2)
    gmsh.write(output_path)
    gmsh.finalize()
