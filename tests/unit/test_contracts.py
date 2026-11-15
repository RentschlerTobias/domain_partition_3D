import json

import numpy as np
import pytest

from hexdp.contracts import Context, ParamSpec, Score, describe, resolve_params, validate
from hexdp.pipeline.chain import Stage, run_chain
from tests.unit import chain_plugins as cp

ALL = [("detectors", cp.BBoxDetector()), ("strategies", cp.BoxStrategy()), ("strategies", cp.NeverStrategy()),
       ("ops", cp.SplitOp()), ("discretisers", cp.UniformDiscretiser()), ("fillers", cp.LinearFiller()),
       ("exporters", cp.VtkExporter())]


@pytest.mark.parametrize("kind, plugin", ALL, ids=[p.name for _, p in ALL])
def test_example_plugins_meet_their_contract(kind, plugin):
    assert validate(kind, plugin) == []


def test_validate_names_what_is_missing():
    class Half:
        name = "half"
        params = {"x": 1.0}

    problems = validate("ops", Half())
    assert "params not given as ParamSpec: ['x']" in problems
    assert "missing method 'applicable'" in problems and "missing method 'apply'" in problems


def test_resolve_params():
    spec = {"h": ParamSpec(0.1, "size", unit="m", lower=0), "mode": ParamSpec("a", "mode", choices=("a", "b"))}
    assert resolve_params(spec, {"h": 0.2}) == {"h": 0.2, "mode": "a"}
    with pytest.raises(ValueError, match="unknown parameter"):
        resolve_params(spec, {"size": 1})
    with pytest.raises(ValueError, match="outside"):
        resolve_params(spec, {"h": -1})
    with pytest.raises(ValueError, match="not in"):
        resolve_params(spec, {"mode": "c"})


def test_describe_is_json_ready():
    d = describe(cp.UniformDiscretiser())
    assert d["name"] == "uniform" and d["params"]["target_h"] == {"default": 0.25, "description": "target cell size",
                                                                     "unit": "m", "lower": 1e-9}
    json.dumps(d)


def test_score_no():
    s = Score.no("no blade", missing=("blade",))
    assert not s.ok and s.reasons == ["no blade"] and s.missing == ["blade"]


@pytest.fixture(scope="module")
def box_step(tmp_path_factory):
    gmsh = pytest.importorskip("gmsh")
    p = tmp_path_factory.mktemp("box") / "box.step"
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.occ.addBox(0, 0, 0, 2, 1, 1)
    gmsh.model.occ.synchronize()
    gmsh.write(str(p))
    gmsh.finalize()
    return p


def test_whole_chain_on_a_box(box_step, tmp_path):
    from hexdp.io.gmsh_model import load

    geo = load(box_step)
    ctx = Context(tmp_path / "ws")
    res = run_chain(geo, ctx,
                    detectors=[Stage(cp.BBoxDetector())],
                    strategies=[Stage(cp.NeverStrategy()), Stage(cp.BoxStrategy())],
                    ops=[Stage(cp.SplitOp())],
                    discretiser=Stage(cp.UniformDiscretiser(), {"target_h": 0.25}),
                    filler=Stage(cp.LinearFiller()),
                    exporters=[Stage(cp.VtkExporter())])
    assert res.strategy == "box" and not res.scores["never"].ok
    s = res.blocks.summary()
    assert (s["blocks"], s["fine_cells"], s["issues"]) == (2, 2 * (4 * 4 * 4), {})
    assert s["min_sj_fine"] == pytest.approx(1.0)
    log = json.loads((tmp_path / "ws" / "decisions.json").read_text())
    stages = [d["stage"] for d in log]
    assert stages == ["detect", "strategy", "op", "export"]
    strat = log[1]
    assert strat["choice"] == "box" and strat["alternatives"] == {"never": 0.0, "box": 0.5}
    assert (tmp_path / "ws" / "vtk.vtk").exists()


def test_chain_explains_when_nothing_applies(box_step, tmp_path):
    from hexdp.io.gmsh_model import load

    with pytest.raises(RuntimeError, match="no applicable strategy.*no blade found.*needs blade feature"):
        run_chain(load(box_step), Context(tmp_path), strategies=[Stage(cp.NeverStrategy())],
                  discretiser=Stage(cp.UniformDiscretiser()), filler=Stage(cp.LinearFiller()))
