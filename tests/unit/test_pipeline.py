import json
import time

import pytest
import yaml

from hexdp.pipeline.case import CaseError, load_case
from hexdp.pipeline.evaluate import evaluate, markdown
from hexdp.pipeline.runner import run_case, status

P = "tests.unit.chain_plugins:"


@pytest.fixture(scope="module")
def box_step(tmp_path_factory):
    gmsh = pytest.importorskip("gmsh")
    p = tmp_path_factory.mktemp("geo") / "box.step"
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.occ.addBox(0, 0, 0, 2, 1, 1)
    gmsh.model.occ.synchronize()
    gmsh.write(str(p))
    gmsh.finalize()
    return p


def write_case(d, **kw):
    data = {"name": "box", "geometry": "box.step", "strategy": [P + "BoxStrategy"], "detectors": [P + "BBoxDetector"],
            "ops": [P + "SplitOp"], "discretiser": {"name": P + "UniformDiscretiser", "params": {"target_h": 0.25}},
            "filler": P + "LinearFiller", "exporters": [P + "VtkExporter"], "reference": {"blocks": 2, "min_sj": 1.0}}
    data.update(kw)
    f = d / f"{data['name']}.yaml"
    f.write_text(yaml.safe_dump(data))
    return f


@pytest.fixture
def case_dir(tmp_path, box_step):
    (tmp_path / "box.step").write_bytes(box_step.read_bytes())
    return tmp_path


def test_case_paths_and_env(case_dir, monkeypatch):
    monkeypatch.setenv("HEXDP_TEST_DIR", str(case_dir))
    f = write_case(case_dir, geometry="${HEXDP_TEST_DIR}/box.step", params={P + "SplitOp": {}},
                   resources={"cores": 4, "cpu_list": "0-3"})
    c = load_case(f)
    assert c.geometry == case_dir / "box.step" and c.workspace == case_dir / "runs" / "box"
    assert (c.cores, c.cpu_list) == (4, "0-3") and c.discretiser.params == {"target_h": 0.25}
    c2 = load_case(write_case(case_dir, name="p", ops=[{"name": P + "SplitOp", "params": {"d": "${HEXDP_TEST_DIR}/x"}}]))
    assert c2.ops[0].params == {"d": f"{case_dir}/x"}


def test_invalid_yaml_is_a_case_error(case_dir):
    f = case_dir / "bad.yaml"
    f.write_text("name: x\nreuse: [${A}/b]\n")
    with pytest.raises(CaseError, match="not valid YAML"):
        load_case(f)


@pytest.mark.parametrize("change, msg", [
    ({"strategie": "x"}, "unknown key"),
    ({"geometry": "${HEXDP_NOT_SET_XYZ}/box.step"}, "not set"),
    ({"geometry": "missing.step"}, "geometry not found"),
    ({"ops": [{"params": {}}]}, "expected a plugin name"),
    ({"params": {P + "SplitOp": {"x": "${HEXDP_NOT_SET_XYZ}"}}, "ops": [P + "SplitOp"]}, "not set"),
])
def test_case_errors(case_dir, change, msg):
    with pytest.raises(CaseError, match=msg):
        load_case(write_case(case_dir, **change))


def test_run_case_writes_report_status_and_decisions(case_dir):
    rep = run_case(load_case(write_case(case_dir)))
    ws = case_dir / "runs" / "box"
    assert rep["state"] == "done" and rep["strategy"] == "box"
    assert rep["blocks"]["blocks"] == 2 and rep["reference"]["blocks"]["ok"] and rep["reference"]["min_sj"]["ok"]
    assert status(ws)["state"] == "done"
    assert json.loads((ws / "report.json").read_text())["timings"]["strategy:box"] >= 0
    assert (ws / "decisions.json").exists() and (ws / "blocks.json").exists()
    assert rep["outputs"]["vtk"][0].endswith("out/vtk.vtk")


def test_failed_run_is_reported_not_raised(case_dir):
    rep = run_case(load_case(write_case(case_dir, name="never", strategy=[P + "NeverStrategy"])))
    assert rep["state"] == "failed" and "no applicable strategy" in rep["error"]
    assert status(case_dir / "runs" / "never")["state"] == "failed"


def test_roles_from_the_case(case_dir):
    f = write_case(case_dir, name="roles", roles={1: "inlet"})
    rep = run_case(load_case(f))
    assert rep["geometry_summary"]["roles"]["inlet"] == [1]
    with pytest.raises(CaseError, match="unknown key"):
        load_case(write_case(case_dir, name="bad", role={1: "inlet"}))


def test_eval_reuses_finished_reports(case_dir):
    f1 = write_case(case_dir)
    f2 = write_case(case_dir, name="wrongref", reference={"blocks": 3})
    rows = evaluate([f1, f2])
    assert [r["reference"] for r in rows] == ["ok", "FAIL blocks"]
    t = (case_dir / "runs" / "box" / "report.json").stat().st_mtime
    evaluate([f1])                                      # finished: not run again
    assert (case_dir / "runs" / "box" / "report.json").stat().st_mtime == t
    assert "| box | done | box | 2 |" in markdown(rows)


def test_background_job(case_dir):
    from hexdp.pipeline.jobs import start

    ws = start(write_case(case_dir, name="bg"))
    for _ in range(120):
        if status(ws)["state"] in ("done", "failed"):
            break
        time.sleep(0.5)
    assert status(ws)["state"] == "done", (ws / "run.out").read_text()
