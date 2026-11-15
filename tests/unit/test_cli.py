import pytest

from hexdp import __version__, cli


def test_version(capsys):
    with pytest.raises(SystemExit) as stop:
        cli.main(["--version"])
    assert stop.value.code == 0
    assert capsys.readouterr().out.strip() == f"hexdp {__version__}"


def test_plugins_lists_every_kind(capsys):
    assert cli.main(["plugins"]) == 0
    out = capsys.readouterr().out
    for kind in ("detectors", "strategies", "ops", "discretisers", "fillers", "exporters", "interfaces"):
        assert f"{kind}:" in out


def test_no_command_prints_help(capsys):
    assert cli.main([]) == 0
    assert "usage: hexdp" in capsys.readouterr().out


def test_describe_prints_json(capsys):
    import json

    assert cli.main(["plugins", "ops", "--describe"]) == 0
    ops = json.loads(capsys.readouterr().out)["ops"]
    assert all({"name", "doc", "params"} <= set(o) for o in ops)
