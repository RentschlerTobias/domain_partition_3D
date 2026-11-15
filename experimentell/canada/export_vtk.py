"""Export the canada machine, the blade cut-outs and the block structures as VTU for ParaView.

out: /work/canada/vtk/{01_machine,02_blade_cut,03_blocks}
"""
import shutil, sys, collections
from pathlib import Path
import numpy as np
import meshio

HEX = Path("/opt/stack/domain_partition_3D/experimentell/hex3d_algohex")
sys.path.insert(0, str(HEX.parent.parent)); sys.path.insert(0, str(HEX))
from dp3d.extraction import parse_msh          # noqa: E402
import envelope_cut as ec                      # noqa: E402
import tet_prep as tp                          # noqa: E402

C = Path("/work/canada"); E = C / "envelope"; O = C / "vtk"
for d in ("01_machine", "02_blade_cut", "03_blocks"):
    (O / d).mkdir(parents=True, exist_ok=True)

COMP = {"blk (inlet, duct outside, diffuser, outlet)": 1, "guide vane passage": 2, "runner passage": 3}
VOL = {"tetra", "hexahedron", "wedge", "pyramid"}


def volume_part(m, comp):
    cells, geo = [], []
    for c, g in zip(m.cells, m.cell_data["gmsh:geometrical"]):
        if c.type in VOL:
            cells.append((c.type, c.data)); geo.append(np.asarray(g))
    return m.points, cells, geo, comp


def rot(P, ang):
    c, s = np.cos(ang), np.sin(ang)
    Q = P.copy(); Q[:, 0] = c * P[:, 0] - s * P[:, 1]; Q[:, 1] = s * P[:, 0] + c * P[:, 1]
    return Q


def write(parts, path):
    pts, cells, cd = [], collections.defaultdict(list), collections.defaultdict(lambda: collections.defaultdict(list))
    off = 0
    for P, cl, geo, comp, copy in parts:
        pts.append(P)
        for (t, d), g in zip(cl, geo):
            cells[t].append(d + off)
            cd[t]["component"].append(np.full(len(d), comp))
            cd[t]["entity"].append(g)
            cd[t]["copy"].append(np.full(len(d), copy))
        off += len(P)
    types = list(cells)
    m = meshio.Mesh(np.vstack(pts), [(t, np.vstack(cells[t])) for t in types],
                    cell_data={k: [np.concatenate(cd[t][k]) for t in types] for k in ("component", "entity", "copy")})
    m.write(path)
    print("wrote", path, sum(len(np.vstack(cells[t])) for t in types), "cells")


# ---------------------------------------------------------------- 01 machine
blk = volume_part(meshio.read(C / "full_blk.msh"), 1)
gv = volume_part(meshio.read(C / "gv_gridGmsh.msh"), 2)
ru = volume_part(meshio.read(C / "ru_gridGmsh.msh"), 3)
write([(*blk, 0)], O / "01_machine/blk_quarter.vtu")
write([(*gv, 0)], O / "01_machine/guide_vane_passage.vtu")
write([(*ru, 0)], O / "01_machine/runner_passage.vtu")
# full 360: blk quarter mirrored (symmetry planes x=0, y=0), 5 gv and 4 ru passages rotated
P = blk[0]
mir = [P, P * [-1, 1, 1], P * [1, -1, 1], P * [-1, -1, 1]]
parts = [(Q, blk[1], blk[2], 1, i) for i, Q in enumerate(mir)]
parts += [(rot(gv[0], 2 * np.pi * i / 5), gv[1], gv[2], 2, i) for i in range(5)]
parts += [(rot(ru[0], 2 * np.pi * i / 4), ru[1], ru[2], 3, i) for i in range(4)]
write(parts, O / "01_machine/machine_360.vtu")

# ------------------------------------------------------------ 02 blade cut
NAMES = {1: "inlet", 2: "outlet", 3: "periodic_A", 4: "periodic_B", 5: "hub_BL_cut",
         6: "shroud_BL_cut", 7: "ogrid_or_envelope"}


def tet_vtk_surfaces(src, dst_prefix):
    m = meshio.read(src)
    key = next(iter(m.cell_data))
    tri = [(c.data, np.asarray(d).ravel()) for c, d in zip(m.cells, m.cell_data[key]) if c.type == "triangle"]
    T = np.vstack([t for t, _ in tri]); L = np.concatenate([l for _, l in tri])
    meshio.Mesh(m.points, [("triangle", T)], cell_data={"label": [L]}).write(f"{dst_prefix}_boundary.vtu")
    meshio.Mesh(m.points, [("triangle", T[L == 7])]).write(f"{dst_prefix}_label7.vtu")
    lines = [c.data for c in m.cells if c.type == "line"]
    if lines:
        meshio.Mesh(m.points, [("line", np.vstack(lines))]).write(f"{dst_prefix}_feature_curves.vtu")
    print("wrote", dst_prefix, "boundary/label7/feature curves")


for comp, msh, red, flu in (("runner", "ru_gridGmsh.msh", E / "e2", E / "ru_fluid"),
                            ("guide_vane", "gv_gridGmsh.msh", E / "gv1", E / "gv10_d050")):
    d = O / "02_blade_cut" / comp; d.mkdir(exist_ok=True)
    nodes, els = parse_msh(C / msh)
    BV, BF, body, caps = ec.body_wall(nodes, els)
    meshio.Mesh(BV, [("triangle", BF)]).write(d / "blade_wall.vtu")
    # dtOO O-grid region (hexes of the O-grid volumes) and its interface to the rest
    og = [e for e in els if e[0] == 5 and e[2] in ec.v5.OGRID_GEOM]
    ids = sorted({v for e in og for v in e[3]}); mp = {v: i for i, v in enumerate(ids)}
    meshio.Mesh(np.array([nodes[v] for v in ids]), [("hexahedron", np.array([[mp[v] for v in e[3]] for e in og]))]
                ).write(d / "dtoo_ogrid_volume.vtu")
    tet_vtk_surfaces(red / "tet.vtk", d / "reduced")
    tet_vtk_surfaces(flu / "tet.vtk", d / "fluid")

# ------------------------------------------------------------- 03 blocks
for name, src in (("runner_reduced_e2", E / "e2"), ("guide_vane_reduced_gv1_artefact_5blocks", E / "gv1"),
                  ("guide_vane_reduced_gv1_gapguard_12blocks", E / "gv1/gap"),
                  ("guide_vane_fluid_gv10_delta050", E / "gv10_d050")):
    d = O / "03_blocks" / name; d.mkdir(exist_ok=True)
    for f in ("blocks.vtk", "blocks_edges.vtk", "blocks_quality.vtk", "tfi.vtk"):
        if (src / f).exists():
            shutil.copy(src / f, d / f)
    print("copied", name, sorted(p.name for p in d.iterdir()))
