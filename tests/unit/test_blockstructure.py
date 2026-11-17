import numpy as np
import pytest

from hexdp.core.blockstructure import OFF, BlockStructure, invalid
from hexdp.io.blocks_vtk import from_grids, read_blocks_vtk, write_blocks_vtk


def box_blocks(boxes):
    """Blocks from axis-aligned boxes (lo, hi); shared corners get one vertex."""
    pts, idx, blocks = [], {}, []
    for lo, hi in boxes:
        lo, hi = np.asarray(lo, float), np.asarray(hi, float)
        c8 = []
        for o in OFF:
            p = tuple(np.where(o, hi, lo).round(12))
            if p not in idx:
                idx[p] = len(pts)
                pts.append(p)
            c8.append(idx[p])
        blocks.append(c8)
    return BlockStructure(np.array(pts), np.array(blocks))


def kinds(bs):
    return sorted({i.kind for i in bs.check()})


def test_two_cubes_conform():
    bs = box_blocks([((0, 0, 0), (1, 1, 1)), ((1, 0, 0), (2, 1, 1))])
    assert kinds(bs) == []
    s = bs.summary()
    assert (s["blocks"], s["vertices"], s["boundary_faces"]) == (2, 12, 10)
    assert len(bs.faces()) == 11 and len(bs.edges()) == 20
    # i-edges of both blocks are separate sheets? no: block i-edges are not shared, j/k sheets are
    assert sorted(len(sh["members"]) for sh in bs.sheets()) == [1, 1, 2, 2]


def test_flipped_block_is_reported():
    bs = box_blocks([((0, 0, 0), (1, 1, 1))])
    bs.blocks = bs.blocks[:, [4, 5, 6, 7, 0, 1, 2, 3]]
    assert kinds(bs) == ["orientation"]


def test_t_junction_is_reported():
    bs = box_blocks([((0, 0, 0), (1, 1, 1)), ((1, 0, 0), (2, 0.5, 1)), ((1, 0.5, 0), (2, 1, 1))])
    issues = [i for i in bs.check() if i.kind == "hanging_vertex"]
    assert len(issues) == 2       # the two mid vertices of x = 1 lie on block 0's face


def test_sheet_counts_must_agree():
    bs = box_blocks([((0, 0, 0), (1, 1, 1)), ((1, 0, 0), (2, 1, 1))])
    for e, users in bs.edges().items():
        bs.edge_cells[e] = 4
    j_edge = next(e for e, users in bs.edges().items() if users[0][1] == 1 and len(users) == 2)
    bs.edge_cells[j_edge] = 5
    assert "sheet_counts" in kinds(bs)


def test_json_roundtrip(tmp_path):
    bs = box_blocks([((0, 0, 0), (1, 1, 1)), ((1, 0, 0), (2, 1, 1))])
    e = next(iter(bs.edges()))
    bs.edge_cells[e] = 7
    bs.edge_grading[e] = 1.2
    bs.edge_curves[e] = np.array([[0, 0, 0], [0.5, 0.1, 0], [1, 0, 0]], float)
    bs.save(tmp_path / "bs.json")
    back = BlockStructure.load(tmp_path / "bs.json")
    assert np.array_equal(back.blocks, bs.blocks) and back.edge_cells[e] == 7 and back.edge_grading[e] == 1.2
    assert np.allclose(back.edge_curves[e], bs.edge_curves[e])


def _grid(lo, hi, n, P, index):
    """structured node-id grid of a box with n cells per axis, sharing nodes through `index`."""
    axes = [np.linspace(lo[a], hi[a], n[a] + 1) for a in range(3)]
    G = np.empty([k + 1 for k in n], np.int64)
    for i, x in enumerate(axes[0]):
        for j, y in enumerate(axes[1]):
            for k, z in enumerate(axes[2]):
                key = (round(x, 12), round(y, 12), round(z, 12))
                if key not in index:
                    index[key] = len(P)
                    P.append(key)
                G[i, j, k] = index[key]
    return G


def test_vtk_roundtrip_keeps_grids_and_counts(tmp_path):
    P, index = [], {}
    grids = [_grid((0, 0, 0), (1, 1, 1), (2, 3, 4), P, index), _grid((1, 0, 0), (2, 1, 1), (5, 3, 4), P, index)]
    bs = from_grids(np.array(P, float), grids)
    assert kinds(bs) == [] and bs.summary()["fine_cells"] == 2 * 3 * 4 + 5 * 3 * 4
    write_blocks_vtk(bs, tmp_path / "b.vtk")
    back = read_blocks_vtk(tmp_path / "b.vtk")
    assert back.n_blocks == 2 and sorted(back.edge_cells.values()) == sorted(bs.edge_cells.values())
    assert back.summary()["min_sj_fine"] == pytest.approx(1.0)


def test_block_level_t_junction_with_conforming_fine_mesh():
    """Two half-blocks next to one block whose fine grid has a node at the split: a 1:2 block layout
    that is still conforming cell by cell."""
    P, index = [], {}
    grids = [_grid((0, 0, 0), (1, 1, 1), (2, 2, 2), P, index),
             _grid((1, 0, 0), (2, 0.5, 1), (2, 1, 2), P, index),
             _grid((1, 0.5, 0), (2, 1, 1), (2, 1, 2), P, index)]
    bs = from_grids(np.array(P, float), grids)
    found = [i.kind for i in bs.check()]
    assert found.count("t_junction") == 2 and "hanging_vertex" not in found


def test_validity():
    ok = box_blocks([((0, 0, 0), (1, 1, 1)), ((1, 0, 0), (2, 1, 1))])
    assert ok.summary()["valid"] and invalid(ok.check()) == []
    dup = BlockStructure(ok.points, np.vstack([ok.blocks, ok.blocks[:1]]))   # block 0 twice
    assert "nonmanifold_face" in kinds(dup) and not dup.summary()["valid"]


def test_t_junction_alone_is_valid():
    P, index = [], {}
    grids = [_grid((0, 0, 0), (1, 1, 1), (2, 2, 2), P, index),
             _grid((1, 0, 0), (2, 0.5, 1), (2, 1, 2), P, index),
             _grid((1, 0.5, 0), (2, 1, 1), (2, 1, 2), P, index)]
    assert from_grids(np.array(P, float), grids).summary()["valid"]


def test_beam_topology_issues_of_a_clean_hex_mesh():
    from hexdp.ops.beam_collapse import topology_issues
    n = 3                                  # 3x3x3 unit hexes of one cube: one block, no defect
    ids = np.arange((n + 1) ** 3).reshape(n + 1, n + 1, n + 1)
    P = np.array([(i, j, k) for i in range(n + 1) for j in range(n + 1) for k in range(n + 1)], float)
    H = np.array([[ids[i + o[0], j + o[1], k + o[2]] for o in OFF]
                  for i in range(n) for j in range(n) for k in range(n)])
    assert topology_issues(P, H) == {}


def _grid_blocks(nx, ny, skip=()):
    """nx x ny x 1 unit blocks on a lattice, without the cells in skip."""
    from hexdp.core.blockstructure import OFF
    idx = lambda i, j, k: (i * (ny + 1) + j) * 2 + k  # noqa: E731
    P = np.array([(i, j, k) for i in range(nx + 1) for j in range(ny + 1) for k in range(2)], float)
    B = [[idx(i + o[0], j + o[1], o[2]) for o in OFF] for i in range(nx) for j in range(ny) if (i, j) not in skip]
    return BlockStructure(P, np.array(B))


class _Geo:
    def __init__(self, chi, comps=1, closed=True):
        self.t = {"chi": chi, "components": comps, "closed": closed}

    def boundary_topology(self):
        return self.t


def test_boundary_topology_finds_a_hole():
    from hexdp.core.blockstructure import invalid
    full, ring = _grid_blocks(3, 3), _grid_blocks(3, 3, skip={(1, 1)})
    assert full.boundary_topology() == {"chi": 2, "components": 1}
    assert ring.boundary_topology() == {"chi": 0, "components": 1}          # a torus: the middle block is missing
    assert not invalid(full.check(geo=_Geo(2)))
    bad = invalid(ring.check(geo=_Geo(2)))
    assert [i.kind for i in bad] == ["boundary_topology"]
    assert not invalid(ring.check(geo=_Geo(2, closed=False)))                # open geometry boundary: not tested
    assert not invalid(ring.check())                                         # without geometry: not tested


def test_geometry_boundary_topology_of_a_box(tmp_path):
    gmsh = pytest.importorskip("gmsh")
    from hexdp.io.gmsh_model import load
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.occ.addBox(0, 0, 0, 2, 1, 1)
    gmsh.model.occ.synchronize()
    gmsh.write(str(tmp_path / "box.step"))
    gmsh.finalize()
    assert load(tmp_path / "box.step").boundary_topology() == {"chi": 2, "components": 1, "closed": True}


def test_msh_export_groups(tmp_path):
    gmsh = pytest.importorskip("gmsh")
    from hexdp.core.blockstructure import FaceInfo
    from hexdp.core.geometry import PatchRole
    from hexdp.export.msh import write_msh
    bs = _grid_blocks(2, 1)
    for key, users in bs.faces().items():
        if len(users) == 1 and users[0][1] == 0 and users[0][2] == 0:          # the x = 0 face: a wall
            bs.face_info[key] = FaceInfo(PatchRole.WALL, 7)
    for block_groups, nvol in ((False, 1), (True, 2)):
        f = write_msh(bs, tmp_path / f"m{nvol}", block_groups=block_groups)
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(str(f))
        groups = {gmsh.model.getPhysicalName(d, t): (d, t) for d, t in gmsh.model.getPhysicalGroups()}
        n_hex = sum(len(e) for t, e in zip(*gmsh.model.mesh.getElements(3)[:2]) if t == 5)
        quads = {}
        for name, (d, t) in groups.items():
            if d == 2:
                ent = int(gmsh.model.getEntitiesForPhysicalGroup(d, t)[0])
                types, tags, _ = gmsh.model.mesh.getElements(2, ent)
                quads[name] = sum(len(e) for ty, e in zip(types, tags) if ty == 3)
        gmsh.finalize()
        assert n_hex == 2 and sum(1 for d, _ in groups.values() if d == 3) == nvol
        assert quads == {"wall_7": 1, "boundary": 9}


def test_canonical_hex_mesh_ignores_entity_order():
    from hexdp.io.ovm import _hex_rotations, canonical_hex_mesh
    bs = _grid_blocks(3, 2)
    P, H = bs.points, bs.blocks
    rng = np.random.default_rng(0)
    perm = rng.permutation(len(P))                      # renumber the vertices
    inv = np.empty_like(perm)
    inv[perm] = np.arange(len(P))
    rot = _hex_rotations()
    H2 = np.array([inv[h][rot[rng.integers(0, 24)]] for h in H])[rng.permutation(len(H))]
    a, b = canonical_hex_mesh(P, H), canonical_hex_mesh(P[perm], H2)
    assert np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1])
    from hexdp.core.blockstructure import BlockStructure
    assert BlockStructure(*b).corner_jacobians().min() > 0          # rotations keep the orientation
