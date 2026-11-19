"""The polycube strategy explains an empty HexEx result instead of writing a mesh with zero cells."""
from hexdp.strategies.polycube.plugin import _label_note


def test_the_note_counts_the_used_directions(tmp_path):
    (tmp_path / "labeling.txt").write_text("tets 4 vertices 8\n1\n2\n4\n")
    note = _label_note(tmp_path)
    assert "3 of 6 directions" in note and "[1, 2, 4]" in note


def test_the_note_names_labels_outside_the_six_directions(tmp_path):
    (tmp_path / "labeling.txt").write_text("tets 4 vertices 8\n0\n1\n2\n3\n4\n5\n")
    note = _label_note(tmp_path)
    assert "invalid labels [0]" in note


def test_the_note_is_empty_without_a_labeling(tmp_path):
    assert _label_note(tmp_path) == ""
