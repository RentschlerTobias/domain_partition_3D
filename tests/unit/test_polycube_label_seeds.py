"""The labeling search tries fixed seeds and keeps the best one."""
import stat
from types import SimpleNamespace

from hexdp.strategies.polycube.plugin import final_score, label_search


def _ctx(tmp_path):
    return SimpleNamespace(workspace=tmp_path, path=lambda p: tmp_path / p, cpu_list="", timeout_s=60)


def _bin(tmp_path):
    d = tmp_path / "bin"
    d.mkdir()
    p = d / "evolabel_cli"
    p.write_text('#!/bin/sh\necho "labeling $3" > labeling.txt\necho "tets $3" > labeling_on_tets.txt\n'
                 'echo "FINAL_SCORE $((10 - $3))"\n')
    p.chmod(p.stat().st_mode | stat.S_IEXEC)
    return d


def test_the_best_seed_wins_and_is_restored(tmp_path):
    ctx = _ctx(tmp_path)
    best = label_search(ctx, _bin(tmp_path), tmp_path, 3)
    assert best == (3, 7.0)
    assert (tmp_path / "labeling.txt").read_text().strip() == "labeling 3"
    assert (tmp_path / "labeling_on_tets.txt").read_text().strip() == "tets 3"
    assert (tmp_path / "label_seeds" / "s1" / "labeling.txt").read_text().strip() == "labeling 1"
    assert (tmp_path / "evocube_label.log").exists()


def test_one_seed_is_a_single_run(tmp_path):
    ctx = _ctx(tmp_path)
    best = label_search(ctx, _bin(tmp_path), tmp_path, 1)
    assert best == (1, 9.0) and (tmp_path / "labeling.txt").read_text().strip() == "labeling 1"


def test_the_score_line_is_read_from_the_end(tmp_path):
    p = tmp_path / "l.log"
    p.write_text("noise\nFINAL_SCORE 1.5 invalid 0 flip 0\nmore\nFINAL_SCORE 2.5 invalid 3 flip 1\n")
    assert final_score(p) == 2.5


def test_a_log_without_the_line_gives_none(tmp_path):
    p = tmp_path / "l.log"
    p.write_text("Flipped or Degenerate Tet\n")
    assert final_score(p) is None
    assert final_score(tmp_path / "missing.log") is None
