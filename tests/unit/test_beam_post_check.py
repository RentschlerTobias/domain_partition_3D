import numpy as np

from hexdp.core.blockstructure import INVALID_KINDS, invalid
from hexdp.ops.beam_collapse import post_check


def _l_shaped_blocks_vtk(path):
    """Three hexes in a 2x2 square with one cell missing: face-connected, but no structured grid."""
    import meshio

    def vid(i, j, k):
        return (i * 3 + j) * 2 + k

    cells = []
    for cx, cy in ((0, 0), (1, 0), (0, 1)):
        cells.append([vid(cx, cy, 0), vid(cx + 1, cy, 0), vid(cx + 1, cy + 1, 0), vid(cx, cy + 1, 0),
                      vid(cx, cy, 1), vid(cx + 1, cy, 1), vid(cx + 1, cy + 1, 1), vid(cx, cy + 1, 1)])
    points = np.array([[x, y, z] for x in range(3) for y in range(3) for z in range(2)], float)
    meshio.write(path, meshio.Mesh(points, [("hexahedron", np.array(cells))],
                                   cell_data={"block": [np.zeros(len(cells), int)]}))


def test_a_block_that_is_no_structured_grid_is_an_issue_not_a_crash(tmp_path):
    p = tmp_path / "blocks.vtk"
    _l_shaped_blocks_vtk(p)
    bad = post_check(p)
    assert [i.kind for i in bad] == ["not_structured"]
    assert "structured grid" in bad[0].message
    assert bad == invalid(bad) and "not_structured" in INVALID_KINDS
