"""The MCP server through a real stdio client: tool listing, a full job on a box with the example
plugins, images, and the path guards."""
import asyncio
import json
import os
import sys
from pathlib import Path

import pytest

mcp = pytest.importorskip("mcp")
from mcp.client.session import ClientSession  # noqa: E402
from mcp.client.stdio import StdioServerParameters, stdio_client  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
P = "tests.unit.chain_plugins:"


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    gmsh = pytest.importorskip("gmsh")
    d = tmp_path_factory.mktemp("data")
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.occ.addBox(0, 0, 0, 2, 1, 1)
    gmsh.model.occ.synchronize()
    gmsh.write(str(d / "box.step"))
    gmsh.finalize()
    for name in ("sample_1", "sample_2"):                       # a geometry family for the workflow test
        (d / "fam" / name).mkdir(parents=True)
        (d / "fam" / name / "box.step").write_bytes((d / "box.step").read_bytes())
    (d.parent / "elsewhere").mkdir(exist_ok=True)
    (d.parent / "elsewhere" / "ring.stl").write_text("solid x\nendsolid x\n")
    (d / "linked").symlink_to(d.parent / "elsewhere")          # data roots are often symlinked directories
    (d / "charite_secret").mkdir()
    (d / "charite_secret" / "lv.stl").write_text("solid x\nendsolid x\n")
    return d


def call(session, name, **args):
    async def go():
        return await session.call_tool(name, args)
    return go()


async def _session_run(data, fn, preset_cpus="1"):
    env = dict(os.environ, HEXDP_DATA=str(data), HEXDP_MCP_ROOT=str(data / "runs"),
               PYTHONPATH=os.pathsep.join([str(REPO / "src"), str(REPO)]))
    env.pop("HEXDP_MCP_CPUS", None)
    if preset_cpus:
        env["HEXDP_MCP_CPUS"] = preset_cpus
    params = StdioServerParameters(command=sys.executable, args=["-m", "hexdp.mcp.server"], env=env, cwd=str(REPO))
    async with stdio_client(params) as (r, w):
        async with ClientSession(r, w) as s:
            await s.initialize()
            return await fn(s)


def payload(res):
    assert not res.is_error, res.content
    if res.structured_content is not None:
        return res.structured_content.get("result", res.structured_content)
    return json.loads(res.content[0].text)


def test_full_session_on_a_box(data):
    async def fn(s):
        tools = {t.name for t in (await s.list_tools()).tools}
        assert {"open_geometry", "detect_features", "list_strategies", "run_strategy", "job_status", "job_result",
                "evaluate", "render_blocks", "export", "record_decision"} <= tools
        files = payload(await s.call_tool("list_inputs", {}))["files"]
        assert [f["path"] for f in files] == ["box.step", "fam/sample_1/box.step", "fam/sample_2/box.step",
                                              "linked/ring.stl"]   # charite is not listed
        g = payload(await s.call_tool("open_geometry", {"path": "box.step"}))
        gid = g["geo_id"]
        assert g["summary"]["surfaces"] == 6
        img = await s.call_tool("render_geometry", {"geo_id": gid})
        assert img.content[0].type == "image" and len(img.content[0].data) > 1000
        feats = payload(await s.call_tool("detect_features", {"geo_id": gid}))["features"]
        assert feats == []                                              # a box: no rotation axis, no sweep
        strat = payload(await s.call_tool("list_strategies", {"geo_id": gid}))["strategies"]
        assert {r["name"] for r in strat} >= {"algohex", "sweep", "cgrid"}
        job = payload(await s.call_tool("run_strategy", {
            "geo_id": gid, "strategy": P + "BoxStrategy", "params": {"ref_file": "box.step"},   # relative to the data root
            "ops": [P + "SplitOp"], "detectors": [P + "BBoxDetector"],
            "discretiser": P + "UniformDiscretiser", "discretiser_params": {"target_h": 0.25},
            "filler": P + "LinearFiller"}))
        jid = job["job_id"]
        assert jid.startswith("boxstrategy_")
        assert job["case"]["strategy"][0]["params"]["ref_file"] == str(data / "box.step")
        second = await s.call_tool("run_strategy", {"geo_id": gid, "strategy": P + "BoxStrategy"})
        assert second.is_error and "still running" in second.content[0].text   # one computation at a time
        st = payload(await s.call_tool("wait_job", {"job_id": jid, "max_s": 60}))
        assert st["state"] == "done", payload(await s.call_tool("job_log", {"job_id": jid}))
        res = payload(await s.call_tool("job_result", {"job_id": jid}))
        assert res["blocks"]["blocks"] == 2 and res["bs_id"]
        ev = payload(await s.call_tool("evaluate", {"bs_id": res["bs_id"]}))
        assert ev["issues_total"] == 0 and ev["summary"]["min_sj_fine"] == pytest.approx(1.0)
        for tool in ("render_blocks", "render_quality"):
            im = await s.call_tool(tool, {"bs_id": res["bs_id"]})
            assert im.content[0].type == "image"
        ex = payload(await s.call_tool("export", {"bs_id": res["bs_id"], "format": "vtk"}))
        assert Path(ex["written"][0]).exists()
        assert res["valid"] is True
        dup = payload(await s.call_tool("block_op", {"bs_id": res["bs_id"], "op": P + "DuplicateOp"}))
        ev2 = payload(await s.call_tool("evaluate", {"bs_id": dup["bs_id"]}))
        assert ev2["valid"] is False and any(i["kind"] == "nonmanifold_face" for i in ev2["issues"])
        refused = await s.call_tool("export", {"bs_id": dup["bs_id"]})
        assert refused.is_error and "invalid" in refused.content[0].text      # no export of an invalid structure
        forced = payload(await s.call_tool("export", {"bs_id": dup["bs_id"], "allow_invalid": True}))
        assert Path(forced["written"][0]).exists()
        note = payload(await s.call_tool("record_decision", {"job_or_geo_id": jid, "text": "box: one block, split"}))
        assert "split" in Path(note["file"]).read_text()
        return True

    assert asyncio.run(_session_run(data, fn))


def test_guards(data):
    async def fn(s):
        r1 = await s.call_tool("open_geometry", {"path": "charite_secret/lv.stl"})
        r2 = await s.call_tool("open_geometry", {"path": "/etc/hostname"})
        r3 = await s.call_tool("evaluate", {"bs_id": "bs_nope"})
        r4 = await s.call_tool("open_geometry", {"path": "../../../../etc/hostname"})
        r5 = await s.call_tool("list_inputs", {"subdir": "linked"})
        gid = payload(await s.call_tool("open_geometry", {"path": "box.step"}))["geo_id"]
        r6 = await s.call_tool("run_strategy", {"geo_id": gid, "strategy": P + "BoxStrategy",
                                                "params": {"ref_file": "charite_secret/lv.stl"}})
        r7 = await s.call_tool("run_strategy", {"geo_id": gid, "strategy": P + "BoxStrategy",
                                                "params": {"ref_file": "nope.stp"}})
        return r1, r2, r3, r4, r5, r6, r7

    r1, r2, r3, r4, r5, r6, r7 = asyncio.run(_session_run(data, fn))
    assert r6.is_error and "confidential" in r6.content[0].text      # file parameters go through the guard
    assert r7.is_error and "FileNotFoundError" in r7.content[0].text
    assert r4.is_error and "outside the allowed roots" in r4.content[0].text
    assert not r5.is_error                                   # a symlinked data directory is readable
    assert r1.is_error and "confidential" in r1.content[0].text
    assert r2.is_error and "outside the allowed roots" in r2.content[0].text
    assert r3.is_error and "unknown bs_id" in r3.content[0].text


def test_workflow_from_a_session(data):
    async def fn(s):
        gid = payload(await s.call_tool("open_geometry", {"path": "fam/sample_1/box.step"}))["geo_id"]
        job = payload(await s.call_tool("run_strategy", {
            "geo_id": gid, "strategy": P + "BoxStrategy", "detectors": [P + "BBoxDetector"],
            "discretiser": P + "UniformDiscretiser", "discretiser_params": {"target_h": 0.5},
            "filler": P + "LinearFiller"}))["job_id"]
        assert payload(await s.call_tool("wait_job", {"job_id": job, "max_s": 60}))["state"] == "done"
        bs = payload(await s.call_tool("job_result", {"job_id": job}))["bs_id"]
        split = payload(await s.call_tool("block_op", {"bs_id": bs, "op": P + "SplitOp"}))["bs_id"]
        await s.call_tool("record_decision", {"job_or_geo_id": job, "text": "one block, then split_i"})
        dup = payload(await s.call_tool("block_op", {"bs_id": split, "op": P + "DuplicateOp"}))["bs_id"]
        bad = await s.call_tool("save_workflow", {"bs_id": dup, "name": "boxes"})
        assert bad.is_error and "invalid" in bad.content[0].text            # only from a valid structure
        nope = await s.call_tool("save_workflow", {"bs_id": split, "name": "boxes", "placeholders": {"sample": "zz"}})
        assert nope.is_error and "does not occur" in nope.content[0].text
        wf = payload(await s.call_tool("save_workflow", {"bs_id": split, "name": "boxes",
                                                         "placeholders": {"sample": "sample_1"},
                                                         "accept": {"valid": True, "blocks": 2}}))
        assert "${HEXDP_DATA}/fam/${sample}/box.step" in wf["text"] and P + "SplitOp" in wf["text"]
        assert "one block, then split_i" in wf["text"]
        again = await s.call_tool("save_workflow", {"bs_id": split, "name": "boxes"})
        assert again.is_error and "overwrite" in again.content[0].text
        guard = await s.call_tool("validate_workflow", {"workflow": "boxes", "runs": [{"sample": "../charite_secret"}]})
        assert guard.is_error
        v = payload(await s.call_tool("validate_workflow", {"workflow": "boxes",
                                                            "runs": [{"sample": "sample_1"}, {"sample": "sample_2"}]}))
        st = payload(await s.call_tool("wait_job", {"job_id": v["job_id"], "max_s": 120}))
        assert st["state"] == "done" and st["runs_done"] == 2, st
        res = payload(await s.call_tool("job_result", {"job_id": v["job_id"]}))
        return wf, res

    wf, res = asyncio.run(_session_run(data, fn))
    assert res["accepted"] == 2 and [r["blocks"] for r in res["rows"]] == [2, 2], res
    assert [r["parameters"]["sample"] for r in res["rows"]] == ["sample_1", "sample_2"]
    assert all(r["accept"] == "ok" and r["valid"] for r in res["rows"])
    # the saved workflow runs without the MCP server
    from hexdp.pipeline.case import load_case
    os.environ["HEXDP_DATA"] = str(data)
    try:
        c = load_case(wf["workflow"], {"sample": "sample_2"})
    finally:
        del os.environ["HEXDP_DATA"]
    assert c.geometry == data / "fam" / "sample_2" / "box.step" and c.accept == {"valid": True, "blocks": 2}


def test_cpus_must_be_given_without_preset(data):
    """No CPU default: a computing job without cpus is refused with the hint to ask the user; check_environment
    reports the hardware and every backend."""
    async def fn(s):
        env = payload(await s.call_tool("check_environment", {}))
        assert env["hardware"]["cpus_available"] >= 1 and "algohex" in env["backends"]
        assert all("ok" in b and "needed_for" in b for b in env["backends"].values())
        gid = payload(await s.call_tool("open_geometry", {"path": "box.step"}))["geo_id"]
        r = await s.call_tool("run_strategy", {"geo_id": gid, "strategy": P + "BoxStrategy"})
        assert r.is_error and "Ask the user" in r.content[0].text
        bad = await s.call_tool("run_strategy", {"geo_id": gid, "strategy": P + "BoxStrategy", "cpus": "100000"})
        assert bad.is_error and "available" in bad.content[0].text
        job = payload(await s.call_tool("run_strategy", {"geo_id": gid, "strategy": P + "BoxStrategy", "cpus": "1", "params": {"ref_file": "box.step"}, "detectors": [P + "BBoxDetector"],
                                                         "discretiser": P + "UniformDiscretiser",
                                                         "discretiser_params": {"target_h": 0.5},
                                                         "filler": P + "LinearFiller"}))
        assert job["case"]["resources"]["cores"] == 1
        st = payload(await s.call_tool("wait_job", {"job_id": job["job_id"], "max_s": 60}))
        assert st["state"] == "done", payload(await s.call_tool("job_log", {"job_id": job["job_id"]}))
    asyncio.run(_session_run(data, fn, preset_cpus=None))


def test_parse_cpus():
    from hexdp.backends import compact, cpus_available, parse_cpus
    avail = cpus_available()
    assert parse_cpus("1") == avail[:1]
    assert parse_cpus(compact(avail[:2])) == avail[:2]
    assert compact([0, 1, 2, 5, 7, 8]) == "0-2,5,7-8"
    with pytest.raises(ValueError):
        parse_cpus(str(len(avail) + 1))
