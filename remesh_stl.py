import gmsh
import gmshpy
import sys
// -----------------------------------------------------------------------------
//
// Gmsh GEO tutorial 13
//
// Remeshing an STL file without an underlying CAD model
//
// -----------------------------------------------------------------------------

sys.path.insert(0, "/mnt/opt.net/src/dtOO/ThirdParty~opensuse-leap~15.6/tools")


// Let's merge an STL mesh that we would like to remesh.
Merge "~/strent/Case199_pre_all.stl"


def remesh_stl(
    input_stl: str,
    output_mesh: str = "output.msh",
    angle: float = 40.0,
    mesh_size: float = 4.0,
    force_parametrizable_patches: bool = False,
    include_boundary: bool = True,
    curve_angle: float = 180.0,
    create_volume: bool = True,
):
    """
    Remesht eine STL-Datei mit gmsh.

    Args:
        input_stl: Pfad zur STL-Eingabedatei
        output_mesh: Pfad zur Ausgabedatei (.msh, .vtk, etc.)
        angle: Winkel (Grad) ab dem eine Kante als scharf gilt (20-120)
        mesh_size: Zielgröße der Mesh-Elemente
        force_parametrizable_patches: Erzwingt parametrisierbare Patches (für komplexe Geometrien)
        include_boundary: Randkanten in die Klassifikation einbeziehen
        curve_angle: Winkel zum Aufteilen von Kurven
        create_volume: Volumen aus den Oberflächen erstellen
    """
    gmsh.initialize()
    gmsh.model.add("remesh")

    # STL einlesen
    gmsh.merge(input_stl)

    # Oberflächen klassifizieren (entlang scharfer geometrischer Merkmale aufteilen)
    import math
    gmsh.model.mesh.classifySurfaces(
        angle * math.pi / 180.0,
        include_boundary,
        force_parametrizable_patches,
        curve_angle * math.pi / 180.0,
    )

    # Geometrie für alle diskreten Kurven und Oberflächen erstellen
    gmsh.model.mesh.createGeometry()

    # Volumen erstellen
    if create_volume:
        surfaces = gmsh.model.getEntities(dim=2)
        surface_tags = [s[1] for s in surfaces]

        surface_loop = gmsh.model.geo.addSurfaceLoop(surface_tags)
        gmsh.model.geo.addVolume([surface_loop])
        gmsh.model.geo.synchronize()

    # Mesh-Größenfeld definieren
    field_tag = gmsh.model.mesh.field.add("MathEval")
    gmsh.model.mesh.field.setString(field_tag, "F", str(mesh_size))
    gmsh.model.mesh.field.setAsBackgroundMesh(field_tag)

    # Meshing durchführen
    gmsh.model.mesh.generate(3 if create_volume else 2)

    # Speichern
    gmsh.write(output_mesh)
    print(f"Mesh gespeichert: {output_mesh}")

    gmsh.finalize()


if __name__ == "__main__":
    # Beispielaufruf
    input_file = sys.argv[1] if len(sys.argv) > 1 else "input.stl"
    output_file = sys.argv[2] if len(sys.argv) > 2 else "output.msh"

    remesh_stl(
        input_stl=input_file,
        output_mesh=output_file,
        angle=40.0,          # Schärfe-Winkel für Kantenerkennung
        mesh_size=4.0,        # Elementgröße
        create_volume=True,   # 3D-Volumen erzeugen (False = nur Oberfläche)
    )
