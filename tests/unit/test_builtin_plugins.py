import pytest

from hexdp.contracts import Context, Feature, describe, resolve_params, validate
from hexdp.export.plugins import FoamExporter, VtkExporter
from hexdp.fill.plugins import B5Filler, KeepDiscretiser, TfiFiller
from hexdp.ops.plugins import BlockCollapseOp
from hexdp.strategies.algohex.plugin import AlgoHexStrategy
from hexdp.strategies.sweep.plugin import CGridStrategy, SweepStrategy

BUILTIN = [("strategies", AlgoHexStrategy()), ("strategies", SweepStrategy()), ("strategies", CGridStrategy()),
           ("ops", BlockCollapseOp()), ("discretisers", KeepDiscretiser()), ("fillers", TfiFiller()),
           ("fillers", B5Filler()), ("exporters", VtkExporter()), ("exporters", FoamExporter())]


@pytest.mark.parametrize("kind, plugin", BUILTIN, ids=[p.name for _, p in BUILTIN])
def test_builtin_plugins_meet_their_contract(kind, plugin):
    assert validate(kind, plugin) == []
    d = describe(plugin)
    assert d["doc"] and all(v["description"] for v in d["params"].values())
    resolve_params(plugin.params, {})          # defaults are consistent with their own bounds


@pytest.fixture(scope="module")
def box(tmp_path_factory):
    gmsh = pytest.importorskip("gmsh")
    d = tmp_path_factory.mktemp("box")
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
    gmsh.model.occ.synchronize()
    gmsh.write(str(d / "box.step"))
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.3)
    gmsh.model.mesh.generate(2)
    gmsh.write(str(d / "box.msh"))
    gmsh.finalize()
    from hexdp.io.gmsh_model import load
    return load(d / "box.step"), load(d / "box.msh")


def test_applicable_says_what_is_missing(box):
    cad, msh = box
    ah = AlgoHexStrategy().applicable(cad, [])
    assert not ah.ok and "dtOO .msh input" in ah.missing
    assert AlgoHexStrategy().applicable(msh, []).ok
    sw = SweepStrategy().applicable(cad, [])
    assert sw.ok and sw.value < 0.2 and any("sweep feature" in m for m in sw.missing)
    assert SweepStrategy().applicable(cad, [Feature("sweep", confidence=1.0, reasons=["similar sections"])]).value == 0.8
    assert not SweepStrategy().applicable(msh, []).ok
    cg = CGridStrategy().applicable(cad, [])
    assert cg.ok and any("thin_wall" in m for m in cg.missing)


def test_sector_angle_is_restricted():
    with pytest.raises(ValueError, match="not in"):
        resolve_params(SweepStrategy.params, {"phi_max": 120.0})


def test_b5_names_its_missing_inputs(tmp_path):
    with pytest.raises(ValueError, match="blocks_b3.vtk"):
        B5Filler().fill(None, None, resolve_params(B5Filler.params, {"run_dir": str(tmp_path)}), Context(tmp_path / "ws"))


def test_cgrid_needs_the_wall_cad(box, tmp_path):
    cad, _ = box
    with pytest.raises(ValueError, match="wall_cad"):
        CGridStrategy().run(cad, [], resolve_params(CGridStrategy.params, {}), Context(tmp_path))
