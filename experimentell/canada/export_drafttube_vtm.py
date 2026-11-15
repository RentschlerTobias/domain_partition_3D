"""AlgoHex input of the canada blade-free domain (draft tube + inlet + duct outside) as a multiblock,
plus the raw AlgoHex hex meshes of the two failed runs (n=2000, n=20000).

out: /work/canada/vtk/drafttube/drafttube.vtm
"""
import collections, sys
from pathlib import Path

import numpy as np
import meshio
from scipy.spatial import cKDTree

HEX = "/opt/stack/domain_partition_3D/experimentell/hex3d_algohex"
sys.path.insert(0, HEX); sys.path.insert(0, "/opt/stack/domain_partition_3D")
import ovm_io                                        # noqa: E402
sys.path.insert(0, "/work/canada")
import tet_prep_named as tn                          # noqa: E402

C = Path("/work/canada"); D = C / "dt"
out = C / "vtk" / "drafttube"
(out / "parts").mkdir(parents=True, exist_ok=True)
tree = collections.OrderedDict()


def add(group, name, mesh):
    f = f"parts/{group}__{name}.vtu"
    mesh.write(out / f)
    tree.setdefault(group, []).append((name, f))


def sub(points, cells, cell_data=None):
    cells = [(t, c) for t, c in cells if len(c)]
    used = np.unique(np.concatenate([c.ravel() for _, c in cells]))
    mp = -np.ones(len(points), int); mp[used] = np.arange(len(used))
    return meshio.Mesh(points[used], [(t, mp[c]) for t, c in cells], cell_data=cell_data)


# ------------------------------------------------------------------ AlgoHex input
src = D / "b2"
t = meshio.read(src / "tet.vtk")
key = next(iter(t.cell_data))
names = {int(l.split()[0]): l.split()[1] for l in open(src / "labels.txt")}
T = np.vstack([c.data for c in t.cells if c.type == "triangle"])
L = np.concatenate([np.asarray(d).ravel() for c, d in zip(t.cells, t.cell_data[key]) if c.type == "triangle"])
for l in np.unique(L):
    add("1_input_surfaces", f"{int(l):02d}_{names[int(l)]}", sub(t.points, [("triangle", T[L == l])]))
# feature curves: label boundaries vs sharp edges inside one label
lines = np.vstack([c.data for c in t.cells if c.type == "line"])
e2l = collections.defaultdict(set)
for tri, l in zip(T, L):
    for i in range(3):
        e2l[tuple(sorted((int(tri[i]), int(tri[(i + 1) % 3]))))].add(int(l))
is_lab = np.array([len(e2l[tuple(sorted(map(int, e)))]) > 1 for e in lines])
add("2_feature_curves", "label_boundaries", sub(t.points, [("line", lines[is_lab])]))
add("2_feature_curves", "sharp_edges_inside_a_label", sub(t.points, [("line", lines[~is_lab])]))
fv = np.concatenate([c.data.ravel() for c in t.cells if c.type == "vertex"])
add("2_feature_curves", "feature_vertices", meshio.Mesh(t.points[fv], [("vertex", np.arange(len(fv))[:, None])]))
# sharp-edge fragments that the anchoring rule dropped (b1 kept them)
t1 = meshio.read(D / "b1" / "tet.vtk")
P1 = t1.points
T1 = np.vstack([c.data for c in t1.cells if c.type == "triangle"])
L1 = np.concatenate([np.asarray(d).ravel() for c, d in zip(t1.cells, t1.cell_data[next(iter(t1.cell_data))])
                     if c.type == "triangle"])
sh_all = tn.sharp_edges(P1, [tuple(x) for x in T1], L1, 40)
e2l1 = collections.defaultdict(set)
for tri, l in zip(T1, L1):
    for i in range(3):
        e2l1[tuple(sorted((int(tri[i]), int(tri[(i + 1) % 3]))))].add(int(l))
lab_feat = [e for e, ls in e2l1.items() if len(ls) > 1]
kept = set(tn.anchored_chains(sh_all, lab_feat))
dropped = np.array([e for e in sh_all if e not in kept])
if len(dropped):
    add("2_feature_curves", "dropped_sharp_fragments_duct_lip", sub(P1, [("line", dropped)]))
add("3_input_volume", "tets", sub(t.points, [("tetra", np.vstack([c.data for c in t.cells if c.type == "tetra"]))]))

# ------------------------------------------------- raw AlgoHex hex meshes (failed runs)
S = t.points[T].mean(1)
tree_s = cKDTree(S)
for run, n in (("b2", 2000), ("b2_n20k", 20000)):
    ovm = D / run / "hex_hex_x.ovm"
    if not ovm.exists():
        continue
    P, e, f, poly = ovm_io.read_ovm(str(ovm))
    H, _ = ovm_io.ovm_to_cells(P, e, f, poly)
    H = np.asarray(H)
    sj = np.array([ovm_io.scaled_jacobian(P, h) for h in H])
    # boundary faces of the hex mesh and their distance to the input surface
    HF = ((0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7))
    cnt = collections.Counter(frozenset(h[list(q)]) for h in H for q in HF)
    bq = np.array([h[list(q)] for h in H for q in HF if cnt[frozenset(h[list(q)])] == 1])
    dist = tree_s.query(P[bq].mean(1))[0]
    add(f"4_algohex_raw_n{n}", "hex_cells", meshio.Mesh(P, [("hexahedron", H)],
                                                         cell_data={"scaled_jacobian": [sj]}))
    add(f"4_algohex_raw_n{n}", "boundary_quads_dist_to_input", meshio.Mesh(
        P, [("quad", bq)], cell_data={"dist_to_input_surface": [dist]}))

with open(out / "drafttube.vtm", "w") as fh:
    fh.write('<?xml version="1.0"?>\n<VTKFile type="vtkMultiBlockDataSet" version="1.0">\n<vtkMultiBlockDataSet>\n')
    for gi, (g, items) in enumerate(tree.items()):
        fh.write(f'  <Block index="{gi}" name="{g}">\n')
        for ii, (n, f) in enumerate(items):
            fh.write(f'    <DataSet index="{ii}" name="{n}" file="{f}"/>\n')
        fh.write('  </Block>\n')
    fh.write('</vtkMultiBlockDataSet>\n</VTKFile>\n')
print("wrote", out / "drafttube.vtm", {g: [n for n, _ in v] for g, v in tree.items()})
