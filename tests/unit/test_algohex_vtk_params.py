from types import SimpleNamespace

import pytest

from hexdp.contracts import Context, resolve_params
from hexdp.strategies.algohex import plugin as ah


@pytest.fixture
def reuse_tet(tmp_path, monkeypatch):
    """A run whose tet.vtk already exists: tet_prep is skipped and the AlgoHex stages are stubbed."""
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "tet.vtk").write_text("")
    (ws / "blocks.vtk").write_text("")
    calls = []
    monkeypatch.setattr(ah, "run_module", lambda ctx, mod, args, log, **kw: calls.append(mod))
    monkeypatch.setattr(ah, "read_blocks_vtk", lambda p: SimpleNamespace(meta={}))
    return SimpleNamespace(ws=ws, geo=SimpleNamespace(path="in.vtk"), ctx=Context(ws), calls=calls)


def _run(reuse_tet, given):
    s = ah.AlgoHexStrategy()
    return s.run(reuse_tet.geo, [], resolve_params(s.params, given), reuse_tet.ctx)


def test_size_max_and_remesh_all_warn_when_tet_prep_is_skipped(reuse_tet, capsys):
    _run(reuse_tet, {"size_max": 0.3, "remesh_all": 0.2})
    warned = [d for d in reuse_tet.ctx.decisions if "ignored" in d.choice]
    assert len(warned) == 1
    assert "size_max=0.3" in warned[0].reasons[0] and "remesh_all=0.2" in warned[0].reasons[0]
    assert "tet_prep" in warned[0].reasons[0]
    assert "[algohex] WARNING:" in capsys.readouterr().out


def test_defaults_do_not_warn_when_tet_prep_is_skipped(reuse_tet):
    _run(reuse_tet, {})
    assert not [d for d in reuse_tet.ctx.decisions if "ignored" in d.choice]
