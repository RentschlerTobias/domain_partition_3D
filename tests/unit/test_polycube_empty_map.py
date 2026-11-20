"""The polycube strategy explains an empty HexEx result instead of writing a mesh with zero cells."""
from hexdp.strategies.polycube.plugin import _label_note


def test_the_note_counts_the_used_directions(tmp_path):
    (tmp_path / "labeling.txt").write_text("0\n2\n4\n")
    note = _label_note(tmp_path)
    assert "3 of 6 directions" in note and "[0, 2, 4]" in note


def test_label_zero_is_a_direction_not_an_invalid_label(tmp_path):
    # the labels are 2*axis + dir, so 0..5 are the six cube directions (evocube flagging_utils.cpp)
    (tmp_path / "labeling.txt").write_text("0\n1\n2\n3\n4\n5\n")
    note = _label_note(tmp_path)
    assert "6 of 6 directions" in note and "invalid" not in note


def test_the_note_names_labels_outside_the_six_directions(tmp_path):
    (tmp_path / "labeling.txt").write_text("0\n2\n9\n")
    note = _label_note(tmp_path)
    assert "invalid labels [9]" in note


def test_the_note_ignores_the_unlabelled_marker_of_the_tet_file(tmp_path):
    (tmp_path / "labeling.txt").write_text("-1\n0\n3\n")
    note = _label_note(tmp_path)
    assert "2 of 6 directions" in note and "invalid" not in note


def test_the_note_is_empty_without_a_labeling(tmp_path):
    assert _label_note(tmp_path) == ""
