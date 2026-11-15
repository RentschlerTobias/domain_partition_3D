"""canada machine as a ParaView multiblock (.vtm): every component its own block.

out: /work/canada/vtk/canada/canada.vtm (+ parts/)
"""
import collections, sys
from pathlib import Path

import numpy as np
import meshio

sys.path.insert(0, "/opt/stack/domain_partition_3D/experimentell/hex3d_algohex")
sys.path.insert(0, "/opt/stack/domain_partition_3D")
import envelope_cut as ec                     # noqa: E402
from dp3d.extraction import parse_msh         # noqa: E402

C = Path("/work/canada"); E = C / "envelope"
out = C / "vtk" / "canada"
(out / "parts").mkdir(parents=True, exist_ok=True)
tree = collections.OrderedDict()
VOLT = ("tetra", "pyramid", "wedge", "hexahedron")


def add(group, name, mesh):
    f = f"parts/{group}__{name}.vtu".replace(" ", "_")
    mesh.write(out / f)
    tree.setdefault(group, []).append((name, f))


def sub(points, cells, R=None):
    cells = [(t, c) for t, c in cells if len(c)]
    used = np.unique(np.concatenate([c.ravel() for _, c in cells]))
    mp = -np.ones(len(points), int); mp[used] = np.arange(len(used))
    Q = points[used] if R is None else points[used] @ R.T
    return meshio.Mesh(Q, [(t, mp[c]) for t, c in cells])


def rot(a):
    return np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])


# ------------------------------------------- blk: inlet, duct outside, diffuser, outlet
b = meshio.read(C / "blk_phys.msh")
names = {v[0]: k for k, v in b.field_data.items()}
vc = [(c.type, c.data) for c in b.cells if c.type in VOLT]
add("1_blk_inlet_mantle_drafttube_outlet", "volume_quarter", sub(b.points, vc))
mir = [np.diag(d) for d in ((1, 1, 1), (-1, 1, 1), (1, -1, 1), (-1, -1, 1))]
patches = collections.defaultdict(list)
for c, ph in zip(b.cells, b.cell_data["gmsh:physical"]):
    if c.type in ("triangle", "quad"):
        ph = np.asarray(ph)
        for p in np.unique(ph):
            patches[names.get(int(p), str(p))].append((c.type, c.data[ph == p]))
for n, cells in sorted(patches.items()):
    if "SYMN" in n:
        add("1_blk_patches_quarter", n, sub(b.points, cells))
        continue
    for i, M in enumerate(mir):                       # full 360 by mirroring the quarter
        add(f"1_blk_patch_{n}", f"quarter{i}", sub(b.points, cells, M))

# ------------------------------------------------ blade passages: guide vane and runner
for comp, msh, nb in (("2_guide_vane", "gv_gridGmsh.msh", 5), ("3_runner", "ru_gridGmsh.msh", 4)):
    m = meshio.read(C / msh)
    P = m.points
    vol, surf = collections.defaultdict(list), collections.defaultdict(list)
    for c, g in zip(m.cells, m.cell_data["gmsh:geometrical"]):
        g = np.asarray(g)
        for e in np.unique(g):
            if c.type in VOLT:
                vol[c.type].append((int(e), c.data[g == e]))
            elif c.type in ("triangle", "quad"):
                surf[int(e)].append((c.type, c.data[g == e]))
    add(f"{comp}_passage_volume", "core_tets", sub(P, [(t, d) for t in ("tetra", "pyramid") for _, d in vol[t]]))
    pr = np.vstack([d for _, d in vol["wedge"]])
    r = np.hypot(*P[pr].mean(1)[:, :2].T); mid = 0.5 * (r.min() + r.max())
    add(f"{comp}_passage_volume", "BL_prisms_hub", sub(P, [("wedge", pr[r < mid])]))
    add(f"{comp}_passage_volume", "BL_prisms_shroud", sub(P, [("wedge", pr[r >= mid])]))
    add(f"{comp}_passage_volume", "dtoo_ogrid_hexes", sub(P, [("hexahedron", np.vstack([d for _, d in vol["hexahedron"]]))]))
    nodes, els = parse_msh(C / msh)
    _, _, body, _ = ec.body_wall(nodes, els)
    sn = {1: "hub_wall", 2: "shroud_wall", 3: "inlet", 4: "outlet", 5: "periodic_B", 6: "periodic_A"}
    for e, n in sn.items():
        add(f"{comp}_passage_surface", n, sub(P, surf[e]))
    add(f"{comp}_passage_surface", "blade_wall", sub(P, [x for e in body for x in surf[e]]))
    walls = [x for e in (1, 2, *body) for x in surf[e]]
    for i in range(nb):
        add(f"{comp}_360_walls", f"passage{i}", sub(P, walls, rot(2 * np.pi * i / nb)))

# ------------------------------------------------ AlgoHex inputs and block structures
LAB = {1: "inlet", 2: "outlet", 3: "periodic_A", 4: "periodic_B", 5: "hub_BL_cut", 6: "shroud_BL_cut",
       7: "blade_cut_envelope"}


def algohex_input(group, src, labnames):
    t = meshio.read(src / "tet.vtk")
    key = next(iter(t.cell_data))
    T = np.vstack([c.data for c in t.cells if c.type == "triangle"])
    L = np.concatenate([np.asarray(d).ravel() for c, d in zip(t.cells, t.cell_data[key]) if c.type == "triangle"])
    for l in np.unique(L):
        add(group, f"surface_{int(l)}_{labnames.get(int(l), 'label')}", sub(t.points, [("triangle", T[L == l])]))
    add(group, "feature_curves", sub(t.points, [("line", np.vstack([c.data for c in t.cells if c.type == "line"]))]))


def blocks(group, src):
    b = meshio.read(src / "blocks.vtk")
    H = np.vstack([c.data for c in b.cells if c.type == "hexahedron"])
    bid = np.concatenate([np.asarray(v).ravel() for k, vs in b.cell_data.items() if k == "block_id" for v in vs])
    for i in np.unique(bid):
        add(group, f"block_{int(i):02d}", sub(b.points, [("hexahedron", H[bid == i])]))
    e = meshio.read(src / "blocks_edges.vtk")
    add(group, "block_edges", sub(e.points, [("line", np.vstack([c.data for c in e.cells if c.type == "line"]))]))


algohex_input("4_algohex_input_runner_reduced", E / "e2", LAB)
algohex_input("4_algohex_input_guide_vane_fluid", E / "gv10_d050", LAB)
dtl = {int(l.split()[0]): l.split()[1] for l in open(C / "dt/b2/labels.txt")}
algohex_input("4_algohex_input_blk_drafttube", C / "dt/b2", dtl)
blocks("5_blocks_runner_reduced_12", E / "e2")
blocks("5_blocks_guide_vane_fluid_12", E / "gv10_d050")

with open(out / "canada.vtm", "w") as fh:
    fh.write('<?xml version="1.0"?>\n<VTKFile type="vtkMultiBlockDataSet" version="1.0">\n<vtkMultiBlockDataSet>\n')
    for gi, (g, items) in enumerate(tree.items()):
        fh.write(f'  <Block index="{gi}" name="{g}">\n')
        for ii, (n, f) in enumerate(items):
            fh.write(f'    <DataSet index="{ii}" name="{n}" file="{f}"/>\n')
        fh.write('  </Block>\n')
    fh.write('</vtkMultiBlockDataSet>\n</VTKFile>\n')
print("wrote", out / "canada.vtm", {g: len(v) for g, v in tree.items()})
