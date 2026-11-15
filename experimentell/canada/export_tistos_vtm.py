"""tistos sample as a ParaView multiblock (.vtm): every component its own block,
so it can be switched on/off in the Pipeline Browser / Multiblock Inspector.

usage: export_tistos_vtm.py <sample_dir> <out_dir>
"""
import collections, sys
from pathlib import Path

import numpy as np
import meshio

src, out = Path(sys.argv[1]), Path(sys.argv[2])
(out / "parts").mkdir(parents=True, exist_ok=True)
tree = collections.OrderedDict()          # group -> [(name, file)]


def add(group, name, mesh):
    f = f"parts/{group}__{name}.vtu".replace(" ", "_")
    mesh.write(out / f)
    tree.setdefault(group, []).append((name, f))


def sub(points, cells):
    """mesh restricted to the used points"""
    used = np.unique(np.concatenate([c.ravel() for _, c in cells]))
    mp = -np.ones(len(points), int); mp[used] = np.arange(len(used))
    return meshio.Mesh(points[used], [(t, mp[c]) for t, c in cells])


# ------------------------------------------------- dtOO passage (ru_gridGmsh)
m = meshio.read(src / "case" / "mesh.msh")
P = m.points
geo = m.cell_data["gmsh:geometrical"]
vol = collections.defaultdict(list)
surf = collections.defaultdict(list)
for c, g in zip(m.cells, geo):
    g = np.asarray(g)
    for e in np.unique(g):
        if c.type in ("tetra", "pyramid", "wedge", "hexahedron"):
            vol[(int(e), c.type)].append(c.data[g == e])
        elif c.type in ("triangle", "quad"):
            surf[int(e)].append((c.type, c.data[g == e]))
# volumes: core (tets + pyramids), boundary-layer prisms split hub/shroud, O-grid hexes per entity
core = [(t, np.vstack(v)) for (e, t), v in vol.items() if t in ("tetra", "pyramid")]
add("1_dtoo_volume", "core_tets", sub(P, core))
pr = np.vstack([np.vstack(v) for (e, t), v in vol.items() if t == "wedge"])
r = np.hypot(*P[pr].mean(1)[:, :2].T)
mid = 0.5 * (r.min() + r.max())
add("1_dtoo_volume", "BL_prisms_hub", sub(P, [("wedge", pr[r < mid])]))
add("1_dtoo_volume", "BL_prisms_shroud", sub(P, [("wedge", pr[r >= mid])]))
for (e, t), v in sorted(vol.items()):
    if t == "hexahedron":
        add("1_dtoo_volume", f"ogrid_hex_entity{e}", sub(P, [(t, np.vstack(v))]))
# surfaces: dtOO geometric entities 1 hub, 2 shroud, 3 inlet, 4 outlet, 5/6 periodic, others = O-grid
names = {1: "hub_wall", 2: "shroud_wall", 3: "inlet", 4: "outlet", 5: "periodic_B", 6: "periodic_A"}
other = []
for e, lst in sorted(surf.items()):
    cells = [(t, d) for t, d in lst]
    if e in names:
        add("2_dtoo_surface", names[e], sub(P, cells))
    else:
        other.append((e, cells))
for e, cells in other:
    add("2_dtoo_surface_ogrid_entities", f"entity{e}", sub(P, cells))

# blade wall = O-grid exterior entities that touch neither the core nor hub/shroud
sys.path.insert(0, "/opt/stack/domain_partition_3D/experimentell/hex3d_algohex")
sys.path.insert(0, "/opt/stack/domain_partition_3D")
import envelope_cut as ec
from dp3d.extraction import parse_msh
_n, _e = parse_msh(src / "case" / "mesh.msh")
_, _, body, _ = ec.body_wall(_n, _e)
add("2_dtoo_surface", "blade_wall", sub(P, [(t, d) for e in body for t, d in surf[e]]))
# full runner, 360 deg: hub, shroud and blade walls (4 passages)
nb = 4
walls = [(t, d) for e in (1, 2, *body) for t, d in surf[e]]
W = sub(P, walls)
for i in range(nb):
    a = 2 * np.pi * i / nb
    R = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
    add("3_runner_360_walls", f"passage{i}", meshio.Mesh(W.points @ R.T, W.cells))

# ------------------------------------------------------ AlgoHex input
t = meshio.read(src / "tet.vtk")
key = next(iter(t.cell_data))
lab_names = {1: "inlet", 2: "outlet", 3: "periodic_A", 4: "periodic_B", 5: "hub_BL_cut",
             6: "shroud_BL_cut", 7: "ogrid_cut_blade"}
T = np.vstack([c.data for c in t.cells if c.type == "triangle"])
L = np.concatenate([np.asarray(d).ravel() for c, d in zip(t.cells, t.cell_data[key]) if c.type == "triangle"])
for l in np.unique(L):
    add("4_algohex_input", f"surface_{int(l)}_{lab_names.get(int(l), 'label')}", sub(t.points, [("triangle", T[L == l])]))
add("4_algohex_input", "feature_curves", sub(t.points, [("line", np.vstack([c.data for c in t.cells if c.type == "line"]))]))
add("4_algohex_input", "tet_volume", sub(t.points, [("tetra", np.vstack([c.data for c in t.cells if c.type == "tetra"]))]))

# ----------------------------------------------------- block structure
b = meshio.read(src / "blocks.vtk")
H = np.vstack([c.data for c in b.cells if c.type == "hexahedron"])
bid = np.concatenate([np.asarray(v).ravel() for k, vs in b.cell_data.items() if k == "block_id" for v in vs])
for i in np.unique(bid):
    add("5_blocks", f"block_{int(i):02d}", sub(b.points, [("hexahedron", H[bid == i])]))
e = meshio.read(src / "blocks_edges.vtk")
add("5_blocks", "block_edges", sub(e.points, [("line", np.vstack([c.data for c in e.cells if c.type == "line"]))]))
if (src / "tfi.vtk").exists():
    f = meshio.read(src / "tfi.vtk")
    add("6_tfi_refill", "tfi_mesh", meshio.Mesh(f.points, f.cells, cell_data=f.cell_data))

# ------------------------------------------------------------- .vtm
with open(out / "tistos.vtm", "w") as fh:
    fh.write('<?xml version="1.0"?>\n<VTKFile type="vtkMultiBlockDataSet" version="1.0">\n<vtkMultiBlockDataSet>\n')
    for gi, (g, items) in enumerate(tree.items()):
        fh.write(f'  <Block index="{gi}" name="{g}">\n')
        for ii, (n, f) in enumerate(items):
            fh.write(f'    <DataSet index="{ii}" name="{n}" file="{f}"/>\n')
        fh.write('  </Block>\n')
    fh.write('</vtkMultiBlockDataSet>\n</VTKFile>\n')
print("wrote", out / "tistos.vtm", {g: len(v) for g, v in tree.items()})
