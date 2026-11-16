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


async def _session_run(data, fn):
    env = dict(os.environ, HEXDP_DATA=str(data), HEXDP_MCP_ROOT=str(data / "runs"), HEXDP_MCP_CPUS="0-3",
               HEXDP_MCP_CORES="4", PYTHONPATH=os.pathsep.join([str(REPO / "src"), str(REPO)]))
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
        assert [f["path"] for f in files] == ["box.step", "linked/ring.stl"]   # charite is not listed
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
