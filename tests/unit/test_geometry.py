import numpy as np
import pytest

from hexdp.core.geometry import CurveKind, PatchRole, RoleSource, role_from_name
from hexdp.io.gmsh_model import load


@pytest.mark.parametrize("name, role", [
    ("RU_INLET", PatchRole.INLET), ("ru_outlet", PatchRole.OUTLET), ("RU_PERIOPS", PatchRole.PERIODIC),
    ("hub wall", PatchRole.WALL), ("RU_BLADE", PatchRole.WALL), ("gv_ru_interface", PatchRole.INTERFACE),
    ("Einlass", PatchRole.INLET), ("symmetryPlane", PatchRole.SYMMETRY), ("surface_7", PatchRole.UNKNOWN),
    (None, PatchRole.UNKNOWN),
])
def test_role_from_name(name, role):
    assert role_from_name(name) is role


@pytest.fixture(scope="module")
def cylinder(tmp_path_factory):
    """r = 0.5, h = 1 cylinder as STEP, as .msh with named physical groups and as STL."""
    gmsh = pytest.importorskip("gmsh")
    d = tmp_path_factory.mktemp("cyl")
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, 1, 0.5)
    gmsh.model.occ.synchronize()
    gmsh.write(str(d / "cyl.step"))
    gmsh.model.addPhysicalGroup(2, [1], name="hub wall")
    gmsh.model.addPhysicalGroup(2, [2], name="RU_INLET")
    gmsh.model.addPhysicalGroup(2, [3], name="outlet")
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.08)
    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.SaveAll", 1)
    gmsh.write(str(d / "cyl.msh"))
    gmsh.write(str(d / "cyl.stl"))
    gmsh.finalize()
    return d


def _kinds(geo):
    return sorted(c.kind.value for c in geo.curves.values())


def test_cad_has_exact_entities_and_a_seam(cylinder):
    g = load(cylinder / "cyl.step")
    assert g.exact and len(g.surfaces) == 3
    assert _kinds(g) == ["seam", "sharp", "sharp"]
    assert all(s.role is PatchRole.UNKNOWN for s in g.surfaces.values())   # STEP carries no roles here


def test_msh_roles_come_from_physical_names(cylinder):
    g = load(cylinder / "cyl.msh")
    assert not g.exact
    assert _kinds(g) == ["seam", "sharp", "sharp"]
    roles = {s.role for s in g.surfaces.values()}
    assert roles == {PatchRole.WALL, PatchRole.INLET, PatchRole.OUTLET}
    assert {s.role_source for s in g.surfaces.values()} == {RoleSource.NAME}


def test_stl_is_segmented_and_says_so(cylinder):
    g = load(cylinder / "cyl.stl")
    assert len(g.surfaces) == 3 and _kinds(g) == ["sharp", "sharp"]
    assert any("STL has no entities" in n for n in g.notes)


def test_manual_role(cylinder):
    g = load(cylinder / "cyl.step")
    g.set_role([2, 3], PatchRole.PERIODIC)
    assert g.surfaces[2].role is PatchRole.PERIODIC and g.surfaces[2].role_source is RoleSource.MANUAL


def test_projection_triangulated_and_exact(cylinder):
    g = load(cylinder / "cyl.step")
    Q = np.array([[0.9, 0.0, 0.5], [0.2, 0.2, 0.5], [0.0, 0.1, 1.4]])
    R, d, owner = g.project(Q)
    Re, de, exact = g.project_exact(Q)
    assert exact.all()
    assert np.allclose(np.linalg.norm(Re[:2, :2], axis=1), 0.5, atol=1e-9)   # on the exact mantle
    assert np.allclose(Re[2], [0.0, 0.1, 1.0], atol=1e-9)                     # on the top cap
    assert np.all(np.abs(d - de) < 0.01)   # triangulation close (inside points sit nearer the chords)


def test_exact_projection_respects_trimming(cylinder):
    """getClosestPoint goes to the untrimmed plane; a point beside the cap must not land off the face."""
    g = load(cylinder / "cyl.step")
    cap = next(t for t, s in g.surfaces.items() if np.ptp(s.nodes[:, 2]) < 1e-9 and s.nodes[0, 2] > 0.5)
    R, _, exact = g.project_exact(np.array([[2.0, 0.0, 1.2]]), surfaces=[cap])
    assert np.hypot(R[0, 0], R[0, 1]) <= 0.5 + 1e-6


def test_vtk_tet_mesh_with_entity_ids(tmp_path):
    """An AlgoHex-style input (gmsh VTK with an integer entity field) loads with its entities."""
    gmsh = pytest.importorskip("gmsh")
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.3)
    gmsh.model.mesh.generate(3)
    gmsh.option.setNumber("Mesh.SaveAll", 1)
    gmsh.write(str(tmp_path / "box.msh"))
    gmsh.finalize()
    import meshio
    m = meshio.read(tmp_path / "box.msh")                  # HexMe layout: one integer field per entity
    cells = [(c.type, c.data) for c in m.cells if c.type in ("vertex", "line", "triangle", "tetra")]
    ids = [d for c, d in zip(m.cells, m.cell_data["gmsh:geometrical"]) if c.type in ("vertex", "line", "triangle", "tetra")]
    meshio.write(tmp_path / "box.vtk", meshio.Mesh(m.points, cells, cell_data={"CellEntityIds": ids}), binary=False)
    g = load(tmp_path / "box.vtk")
    assert len(g.surfaces) == 6 and sorted(c.kind.value for c in g.curves.values()) == ["sharp"] * 12
    assert any("VTK tet mesh" in n for n in g.notes)
    from hexdp.strategies.algohex.plugin import AlgoHexStrategy
    assert "no tet_prep" in AlgoHexStrategy().applicable(g, []).reasons[1]
