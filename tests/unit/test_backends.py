import os
from pathlib import Path

import pytest

from hexdp import backends


@pytest.fixture(autouse=True)
def _fresh_probe():
    backends._probe.cache_clear()
    yield
    backends._probe.cache_clear()


def _binary(directory, name="HexMeshing"):
    directory.mkdir(parents=True, exist_ok=True)
    p = directory / name
    p.write_text("#!/bin/sh\n")
    p.chmod(0o755)
    return p


def test_hexmeshing_bin_wins_and_its_path_is_reported(tmp_path, monkeypatch):
    exe = _binary(tmp_path / "patched")
    monkeypatch.setenv("HEXMESHING_BIN", str(exe))
    r = backends._probe("algohex")
    assert r["ok"] and Path(r["path"]) == exe and r["how"] == "$HEXMESHING_BIN"
    assert "notes" not in r


def test_path_find_is_reported_with_the_hexmeshing_bin_hint(tmp_path, monkeypatch):
    exe = _binary(tmp_path / "bin")
    monkeypatch.delenv("HEXMESHING_BIN", raising=False)
    monkeypatch.setenv("PATH", os.pathsep.join([str(exe.parent), "/usr/bin", "/bin"]))
    r = backends._probe("algohex")
    assert r["ok"] and Path(r["path"]) == exe and r["how"] == "PATH"
    assert "PATH" in r["notes"] and "$HEXMESHING_BIN" in r["notes"]


def test_a_set_but_missing_hexmeshing_bin_is_named(monkeypatch):
    monkeypatch.setenv("HEXMESHING_BIN", "/nonexistent/HexMeshing")
    r = backends._probe("algohex")
    assert not r["ok"] and "$HEXMESHING_BIN=/nonexistent/HexMeshing" in r["why"]
