"""The C4Hex stage: it finds the .ovm c4hex_cli writes and reports failures with the log tail."""
import stat

import pytest

from hexdp.strategies.algohex import extract_c4hex as cx

C4HEX = cx.DEFAULT_BIN
SAMPLE = "/opt/stack/ext/C4HexMeshing/extern/QGP3D/extern/MC3D/tests/resources/minimal.hexex"


def _fake(tmp_path, body):
    p = tmp_path / "c4hex_cli"
    p.write_text("#!/bin/sh\n" + body)
    p.chmod(p.stat().st_mode | stat.S_IEXEC)
    return p


def test_the_written_hex_mesh_is_returned(tmp_path):
    # $4 is the "--output-hex-file" value; the real cli appends _<n>_hex.ovm
    b = _fake(tmp_path, 'touch "${4}_100_hex.ovm"\n')
    ovm = cx.run(tmp_path / "map.hexex", tmp_path, "hex", b)
    assert ovm.name == "hex_100_hex.ovm" and ovm.exists()


def test_a_failed_run_names_the_exit_code_and_the_log_tail(tmp_path, capsys):
    b = _fake(tmp_path, 'echo "IGMUntangler: Boundary hf 542174 of block 12 inverted!" >&2\nexit 22\n')
    with pytest.raises(SystemExit) as exc:
        cx.run(tmp_path / "map.hexex", tmp_path, "hex", b)
    assert exc.value.code == 3
    err = capsys.readouterr().err
    assert "exit 22" in err and "Boundary hf 542174" in err


def test_two_hex_meshes_for_one_run_are_an_error(tmp_path):
    b = _fake(tmp_path, 'touch "${4}_100_hex.ovm" "${4}_200_hex.ovm"\n')
    with pytest.raises(SystemExit) as exc:
        cx.run(tmp_path / "map.hexex", tmp_path, "hex", b)
    assert "2 hex meshes" in str(exc.value)


def test_a_missing_binary_says_how_to_set_it(tmp_path):
    with pytest.raises(SystemExit) as exc:
        cx.run(tmp_path / "map.hexex", tmp_path, "hex", tmp_path / "nope")
    assert "HEXDP_C4HEX_BIN" in str(exc.value)


@pytest.mark.skipif(not __import__("os").path.exists(C4HEX) or not __import__("os").path.exists(SAMPLE),
                    reason="C4HexMeshing not built here")
def test_the_real_binary_turns_a_seamless_map_into_readable_hexes(tmp_path):
    from hexdp.io.ovm import ovm_to_cells, read_ovm

    ovm = cx.run(SAMPLE, tmp_path, "minimal")
    P, E, F, polys = read_ovm(ovm)
    H, skipped = ovm_to_cells(P, E, F, polys)
    assert H.shape[1] == 8 and len(H) > 0 and skipped == 0
