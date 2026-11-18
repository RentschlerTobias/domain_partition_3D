"""tetvtk: a body embedded in a meshed box -> the domain around it and the body alone."""
import numpy as np
import pytest

from hexdp.io import tetvtk


@pytest.fixture(scope="module")
def boxed(tmp_path_factory):
    gmsh = pytest.importorskip("gmsh")
    import meshio
    d = tmp_path_factory.mktemp("tv")
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    box = gmsh.model.occ.addBox(0, 0, 0, 4, 4, 4)
    body = gmsh.model.occ.addBox(1.5, 1.5, 1.5, 1, 1, 1)
    gmsh.model.occ.fragment([(3, box)], [(3, body)])
    gmsh.model.occ.synchronize()
    for dim, tag in gmsh.model.getEntities(3):
        gmsh.model.addPhysicalGroup(3, [tag], tag)
    for dim, tag in gmsh.model.getEntities(2):
        gmsh.model.addPhysicalGroup(2, [tag], tag)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.6)
    gmsh.model.mesh.generate(3)
    gmsh.write(str(d / "boxed.msh"))
    gmsh.finalize()
    m = meshio.read(d / "boxed.msh")
    keep = [i for i, c in enumerate(m.cells) if c.type in ("triangle", "tetra")]
    cells = [m.cells[i] for i in keep]
    ids = [m.cell_data["gmsh:physical"][i] for i in keep]
    mesh = meshio.Mesh(m.points, [(c.type, c.data) for c in cells], cell_data={"CellEntityIds": ids})
    mesh.write(d / "boxed.vtk", file_format="vtk42", binary=False)
    return d / "boxed.vtk"


def test_roles_and_extract(boxed, tmp_path):
    m = tetvtk.read(boxed)
    assert len(m.volumes()) == 2
    roles = tetvtk.volume_roles(m)
    assert len(roles["outer"]) == 1 and len(roles["embedded"]) == 1
    outer = tetvtk.extract_volume(m, roles["outer"])
    body = tetvtk.extract_volume(m, roles["embedded"])
    assert len(outer.tets) + len(body.tets) == len(m.tets)
    assert np.ptp(body.points, 0) == pytest.approx([1, 1, 1])
    assert set(body.tri_ids) <= set(outer.tri_ids)          # the body's faces bound both, with the same ids
    from hexdp.io.gmsh_model import load
    p = tetvtk.write(outer, tmp_path / "outer.vtk")
    q = tetvtk.write(body, tmp_path / "body.vtk")
    assert load(p).boundary_topology() == {"chi": 4, "components": 2, "closed": True}
    assert load(q).boundary_topology() == {"chi": 2, "components": 1, "closed": True}
    assert tetvtk.read(p).volumes() == outer.volumes()


def test_cli(boxed, tmp_path):
    from hexdp.cli import main
    assert main(["extract-volume", str(boxed), str(tmp_path / "f.vtk"), "--role", "outer"]) == 0
    assert tetvtk.read(tmp_path / "f.vtk").volumes()
