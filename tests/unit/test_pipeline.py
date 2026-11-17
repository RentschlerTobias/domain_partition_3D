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


def test_workflow_parameters_and_accept(case_dir):
    (case_dir / "s1").mkdir()
    (case_dir / "s1" / "box.step").write_bytes((case_dir / "box.step").read_bytes())
    f = write_case(case_dir, name="wf_${sample}", geometry="${sample}/box.step", parameters={"sample": "s1"},
                   accept={"valid": True, "blocks": [1, 2], "min_sj": 0.5},
                   reuse=[{"path": "${sample}/hex_hex_x.ovm", "optional": True}], notes=["why"], origin={"job": "j"})
    c = load_case(f)
    assert c.name == "wf_s1" and c.geometry == case_dir / "s1" / "box.step" and c.reuse == []
    with pytest.raises(CaseError, match="unknown parameter"):
        load_case(f, {"nope": "x"})
    with pytest.raises(CaseError, match="geometry not found"):
        load_case(f, {"sample": "s2"})
    rep = run_case(c)
    assert rep["state"] == "done" and rep["accepted"] and rep["parameters"] == {"sample": "s1"}
    assert rep["blocks"]["valid"] is True
    f2 = write_case(case_dir, name="strict", accept={"blocks": 1})
    rep2 = run_case(load_case(f2))
    assert rep2["accepted"] is False and not rep2["accept"]["blocks"]["ok"]
    req = write_case(case_dir, name="req_${sample}", parameters={"sample": None})
    with pytest.raises(CaseError, match="no default"):
        load_case(req)


def test_build_workflow_reproduces_the_case_and_batch(case_dir, tmp_path):
    from hexdp.pipeline.evaluate import markdown
    from hexdp.pipeline.workflow import build_workflow, dump_workflow, parse_each, run_batch

    for s in ("sample_1", "sample_2"):
        (case_dir / "fam" / s).mkdir(parents=True)
        (case_dir / "fam" / s / "box.step").write_bytes((case_dir / "box.step").read_bytes())
    job_case = yaml.safe_load(write_case(case_dir, name="job_1", geometry=str(case_dir / "fam/sample_1/box.step"),
                                         ops=[], workspace=str(case_dir / "jobs/job_1")).read_text())
    wf = build_workflow(job_case, name="boxes", ops=[{"name": P + "SplitOp", "params": {}}],
                        placeholders={"sample": "sample_1"}, data_root=case_dir, notes=["split: one block"])
    assert wf["name"] == "boxes_${sample}" and wf["geometry"] == "${HEXDP_DATA}/fam/${sample}/box.step"
    assert "workspace" not in wf and wf["accept"] == {"valid": True} and wf["notes"] == ["split: one block"]
    f = dump_workflow(wf, tmp_path / "boxes.yaml")
    assert f.read_text().startswith("# hexdp workflow 'boxes'")
    import os
    os.environ["HEXDP_DATA"] = str(case_dir)
    try:
        c = load_case(f)
        assert c.name == "boxes_sample_1" and [o.name for o in c.ops] == [P + "SplitOp"]
        rows = run_batch(f, parse_each("sample=sample_1,sample_2,sample_3"), workspace_root=tmp_path / "v",
                         status_dir=tmp_path / "v")
    finally:
        del os.environ["HEXDP_DATA"]
    assert [r["state"] for r in rows] == ["done", "done", "invalid"]
    assert [r["accept"] for r in rows[:2]] == ["ok", "ok"] and rows[0]["blocks"] == 2
    rep = json.loads((tmp_path / "v" / "report.json").read_text())
    assert rep["kind"] == "batch" and rep["accepted"] == 2 and status(tmp_path / "v")["state"] == "done"
    assert "| boxes_sample_2 | done |" in markdown(rows)
