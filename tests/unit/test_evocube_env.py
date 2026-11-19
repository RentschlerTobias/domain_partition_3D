import os
import subprocess

from hexdp.contracts import Context
from hexdp.strategies.polycube import plugin as pc


def _build(tmp_path):
    build = tmp_path / "evocube"
    (build / "lib" / "gco").mkdir(parents=True)
    (build / "lib" / "OpenVolumeMesh" / "src").mkdir(parents=True)
    (build / "lib" / "libevocube.so").write_text("")
    (build / "lib" / "gco" / "libgco.so").write_text("")
    (build / "lib" / "OpenVolumeMesh" / "src" / "libOpenVolumeMesh.so.2.2").write_text("")
    for name in ("prep_cli", "evolabel_cli", "polycube_withHexEx"):
        (build / name).write_text("")
    return build


def test_lib_dirs_are_found_next_to_the_build(tmp_path):
    build = _build(tmp_path)
    libs = pc._lib_dirs(build)
    assert build / "lib" in libs and build / "lib" / "gco" in libs
    assert build / "lib" / "OpenVolumeMesh" / "src" in libs


def test_lib_dirs_accept_a_bin_subdirectory(tmp_path):
    build = _build(tmp_path)
    bin_dir = build / "bin"
    bin_dir.mkdir()
    libs = pc._lib_dirs(bin_dir)
    assert build / "lib" in libs and build / "lib" / "gco" in libs


def test_run_puts_the_build_libraries_on_the_library_path(tmp_path, monkeypatch):
    build = _build(tmp_path)
    seen = {}

    def fake_run(cmd, **kw):
        seen.update(cmd=cmd, **kw)
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(pc.subprocess, "run", fake_run)
    monkeypatch.setenv("LD_LIBRARY_PATH", "/somewhere/else")
    ctx = Context(tmp_path / "ws")
    pc._run(ctx, [build / "prep_cli", "."], "prep.log", 2, pc._lib_dirs(build))
    path = seen["env"]["LD_LIBRARY_PATH"].split(os.pathsep)
    assert str(build / "lib") in path and str(build / "lib" / "gco") in path
    assert path[-1] == "/somewhere/else"
    assert seen["env"]["OMP_NUM_THREADS"] == "2"
