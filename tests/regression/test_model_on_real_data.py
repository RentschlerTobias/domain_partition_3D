"""The data model on real inputs. Inputs are not versioned; tests skip where a file is absent.
Override the locations with HEXDP_TISTOS_RUN and HEXDP_CANADA_DIR."""
import os
from pathlib import Path

import pytest

from hexdp.io.blocks_vtk import read_blocks_vtk
from hexdp.io.gmsh_model import load

TISTOS = Path(os.environ.get("HEXDP_TISTOS_RUN", "/work/datagen_keepbl/runs/tistos_s1_0131"))
CANADA = Path(os.environ.get("HEXDP_CANADA_DIR", "/work/canada"))


def need(p: Path) -> Path:
    if not p.exists():
        pytest.skip(f"input not available: {p}")
    return p


def test_tistos_0131_blocks():
    bs = read_blocks_vtk(need(TISTOS / "blocks.vtk"))
    s = bs.summary()
    assert (s["blocks"], s["fine_cells"], s["sheets"]) == (12, 582, 10)
    assert s["issues"] == {}
    assert s["min_sj_fine"] == pytest.approx(0.3814, abs=1e-3)


def test_tistos_0131_dtoo_mesh_entities():
    g = load(need(TISTOS / "case" / "mesh.msh"))
    s = g.summary()
    assert (s["surfaces"], s["curves"]) == (31, 52)
    assert s["curve_kinds"] == {"sharp": 42, "smooth": 10}
    assert any("derived from shared mesh nodes" in n for n in g.notes)


def test_canada_domain_brep():
    g = load(need(CANADA / "cad" / "domain_blk.brep"))
    s = g.summary()
    assert g.exact and (s["surfaces"], s["curves"]) == (25, 71)


SWEEP = CANADA / "blocking_study" / "sweep_phi"


def test_canada_draft_tube_sweep():
    s = read_blocks_vtk(need(SWEEP / "sweep_tfi_blocks.vtk")).summary()
    assert (s["blocks"], s["fine_cells"], s["boundary_faces"]) == (72, 217908, 152)
    assert s["min_sj_fine"] == pytest.approx(0.0588, abs=1e-4)
    # fine mesh conforming, but 15 block edges are split by a vertex of the neighbouring blocks
    assert s["issues"] == {"t_junction": 15}


def test_canada_draft_tube_cgrid():
    s = read_blocks_vtk(need(SWEEP / "cgrid" / "sweep_tfi_blocks.vtk")).summary()
    assert (s["blocks"], s["fine_cells"], s["boundary_faces"]) == (100, 292212, 184)
    assert s["min_sj_fine"] == pytest.approx(0.0353, abs=1e-4)
    assert s["issues"] == {}
