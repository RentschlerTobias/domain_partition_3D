"""A valid block structure can still be unusable: the scaled-Jacobian floor is a second question."""
from types import SimpleNamespace

import numpy as np

from hexdp.ops import beam_collapse as bc


def _fake(monkeypatch, per_block):
    def fake_read(path):
        return SimpleNamespace(fine_min_jacobian=lambda: np.array(per_block, float))
    monkeypatch.setattr("hexdp.io.blocks_vtk.read_blocks_vtk", fake_read)


def test_a_healthy_structure_has_no_quality_issue(monkeypatch):
    _fake(monkeypatch, [0.9, 0.7])
    assert bc.quality_issues("blocks.vtk", 0.05) == []


def test_a_nearly_degenerate_block_is_reported(monkeypatch):
    _fake(monkeypatch, [0.9, 2.7e-06])
    issues = bc.quality_issues("blocks.vtk", 0.05)
    assert len(issues) == 1 and issues[0].kind == "low_jacobian"
    assert "block 1" in issues[0].message and "< 0.05" in issues[0].message


def test_a_floor_of_zero_switches_the_gate_off(monkeypatch):
    _fake(monkeypatch, [2.7e-06])
    assert bc.quality_issues("blocks.vtk", 0.0) == []


def test_a_structure_without_fine_grids_is_not_judged(monkeypatch):
    monkeypatch.setattr("hexdp.io.blocks_vtk.read_blocks_vtk",
                        lambda p: SimpleNamespace(fine_min_jacobian=lambda: None))
    assert bc.quality_issues("blocks.vtk", 0.05) == []


def test_a_structure_that_drifted_to_the_tolerance_limit_is_reported():
    issues = bc.hausdorff_issue({"hausdorff_rel": 0.99}, 0.5)
    assert len(issues) == 1 and issues[0].kind == "off_surface" and "0.99" in issues[0].message


def test_a_structure_near_the_surface_passes_the_hausdorff_gate():
    assert bc.hausdorff_issue({"hausdorff_rel": 0.07}, 0.5) == []


def test_a_limit_of_zero_switches_the_hausdorff_gate_off():
    assert bc.hausdorff_issue({"hausdorff_rel": 0.99}, 0.0) == []
