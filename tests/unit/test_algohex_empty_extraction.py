"""HexMeshing exits 0 for a folded integer-grid map and leaves an .ovm without cells."""
from hexdp.strategies.algohex.run_algohex import _hex_cell_count


def test_an_ovm_without_polyhedra_counts_zero(tmp_path):
    pts = "".join(f"{x} {y} {z}\n" for x in (0, 1) for y in (0, 1) for z in (0, 1))
    p = tmp_path / "empty.ovm"
    p.write_text("OVM ASCII\nVertices\n8\n" + pts + "Edges\n0\nFaces\n0\nPolyhedra\n0\n")
    assert _hex_cell_count(p) == 0


def test_a_missing_file_counts_zero(tmp_path):
    assert _hex_cell_count(tmp_path / "does_not_exist.ovm") == 0
