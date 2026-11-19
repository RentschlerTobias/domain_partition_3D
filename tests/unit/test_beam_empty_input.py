"""An empty hex mesh is an input error with a message, not a ValueError from np.vstack([])."""
import meshio
import numpy as np
import pytest

from hexdp.ops.beam_collapse import hexahedra


def test_an_empty_hex_mesh_is_an_error_not_a_crash(capsys):
    mesh = meshio.Mesh(np.zeros((8, 3)), [])
    with pytest.raises(SystemExit) as exc:
        hexahedra(mesh, "hexes.vtk")
    assert exc.value.code == 2
    err = capsys.readouterr().err
    assert "no hexahedra" in err and "hexes.vtk" in err


def test_the_cells_are_returned_as_int():
    cells = np.arange(16).reshape(2, 8)
    mesh = meshio.Mesh(np.zeros((16, 3)), [("hexahedron", cells)])
    H = hexahedra(mesh)
    assert H.shape == (2, 8) and H.dtype == np.int64 and (H == cells).all()
