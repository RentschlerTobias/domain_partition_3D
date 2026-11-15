"""Export inspection VTKs of the extracted geometry (ParaView).

Writes output/hex3d_algohex/vtk/:
  01_boundary_from_2D_elements.vtk  the MSH's own tagged 2D elements, as-is,
                                    cell data = geom id  ("just take the 2D
                                    elements" hypothesis -- see the printed
                                    manifold report for why it is not enough)
  02_boundary_topological.vtk       skin derived from the volume cells,
                                    cell data = physical surface id 1..7
  03_feature_curves.vtk             feature edges + feature vertices
  04_tet_volume.vtk                 the pure tet mesh fed to AlgoHex
  05_hexmesh_v1.vtk                 v1 AlgoHex result (incomplete, 31.4%)
"""

import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "experimentell" / "hex3d_algohex"))

from dp3d.extraction import parse_msh                                 # noqa: E402
import tet_prep as tp                                                 # noqa: E402

MSH = REPO / "data" / "T1_9" / "T1_9_ru_gridGmsh.msh"
OUT = REPO / "output" / "hex3d_algohex" / "vtk"


def write_vtk(path, points, cells, cell_types, cell_data, name, title):
    with open(path, "w") as f:
        f.write("# vtk DataFile Version 2.0\n")
        f.write(f"{title}\n")
        f.write("ASCII\nDATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {len(points)} float\n")
        for p in points:
            f.write(f"{p[0]} {p[1]} {p[2]}\n")
        total = sum(1 + len(c) for c in cells)
        f.write(f"CELLS {len(cells)} {total}\n")
        for c in cells:
            f.write(f"{len(c)} " + " ".join(str(int(v)) for v in c) + "\n")
        f.write(f"CELL_TYPES {len(cells)}\n")
        for t in cell_types:
            f.write(f"{t}\n")
        f.write(f"CELL_DATA {len(cells)}\n")
        f.write(f"SCALARS {name} int 1\nLOOKUP_TABLE default\n")
        for d in cell_data:
            f.write(f"{int(d)}\n")
    print(f"wrote {path}  ({len(points)} pts, {len(cells)} cells)")


def manifold_report(tris, label):
    """A boundary usable as a meshing input must be closed and 2-manifold:
    every edge shared by exactly 2 faces."""
    cnt = Counter()
    for t in tris:
        n = len(t)
        for a in range(n):
            p, q = t[a], t[(a + 1) % n]
            cnt[(min(p, q), max(p, q))] += 1

    h = Counter(cnt.values())
    print(f"  [{label}] edge-multiplicity histogram: {dict(sorted(h.items()))}")
    ok = set(h) == {2}
    print(f"  [{label}] closed 2-manifold: {'YES' if ok else 'NO'}"
          + ("" if ok else
             f"  ({h.get(1,0)} open/border edges, "
             f"{sum(v for k,v in h.items() if k>2)} non-manifold edges)"))
    return ok


def export_2d_elements(nodes, elements):
    """The MSH's own tagged 2D elements (etype 2 tri / 3 quad), untouched.

    Also flags, per geom id, whether that patch lies on the true OUTER skin
    or is an INTERIOR face set (the hex O-grid block interfaces are tagged
    as 2D elements too, which is what makes the raw 2D set non-manifold)."""
    # true outer skin, as face-vertex sets, for the outer/interior test
    outer = {frozenset(f) for f, _cell in tp.volume_exterior_faces(elements)}

    used = sorted({t for et, _p, _g, nds in elements if et in (2, 3)
                   for t in nds})
    nmap = {t: i for i, t in enumerate(used)}
    P = np.array([nodes[t] for t in used])
    cells, types, data, raw, isout = [], [], [], [], []
    per_geom = defaultdict(lambda: [0, 0])   # geom -> [n_outer, n_interior]
    for et, _p, geom, nds in elements:
        if et == 2:
            nd = nds[:3]; vt = 5
        elif et == 3:
            nd = nds[:4]; vt = 9
        else:
            continue
        # a tagged quad matches the skin either directly or as its 2 tris
        o = frozenset(nd) in outer
        if not o and len(nd) == 4:
            o = (frozenset((nd[0], nd[1], nd[2])) in outer
                 or frozenset((nd[0], nd[2], nd[3])) in outer)
        cells.append([nmap[t] for t in nd]); types.append(vt)
        data.append(geom); raw.append(tuple(nd)); isout.append(1 if o else 0)
        per_geom[geom][0 if o else 1] += 1

    print("\n[01] MSH tagged 2D elements (the 'just take the 2D elements' route)")
    print(f"  {len(cells)} 2D elements over {len(per_geom)} geom ids")
    manifold_report(raw, "all 2D elements")
    print("  geom id -> (on outer skin / interior):")
    outer_ids, inner_ids = [], []
    for g in sorted(per_geom):
        no, ni = per_geom[g]
        tagline = "OUTER" if no > ni else "interior"
        (outer_ids if no > ni else inner_ids).append(g)
        print(f"     geom {g:3d}: {no:5d} outer / {ni:5d} interior   -> {tagline}")
    print(f"  => outer geom ids   : {outer_ids}")
    print(f"  => interior geom ids: {inner_ids}  (O-grid block interfaces)")

    keep = [i for i, o in enumerate(isout) if o]
    manifold_report([raw[i] for i in keep], "2D elements, outer only")

    write_vtk(OUT / "01_boundary_from_2D_elements.vtk", P, cells, types, data,
              "geom_id", "T1_9 MSH tagged 2D elements (cell data = geom id)")
    write_vtk(OUT / "01b_2D_elements_outer_flag.vtk", P, cells, types, isout,
              "on_outer_skin",
              "T1_9 tagged 2D elements (1 = on true outer skin, 0 = interior)")


def export_topological(nodes, elements):
    ext = tp.volume_exterior_faces(elements)
    tris = tp.orient_and_triangulate(ext, nodes)
    used = sorted({t for tri in tris for t in tri})
    nmap = {t: i for i, t in enumerate(used)}
    P = np.array([nodes[t] for t in used])
    T = [[nmap[t] for t in tri] for tri in tris]
    print("\n[02] boundary derived topologically from ALL volume cells")
    print(f"  {len(tris)} triangles")
    manifold_report(T, "topological skin")
    ids = tp.classify_boundary(P, T)
    write_vtk(OUT / "02_boundary_topological.vtk", P, T, [5] * len(T), ids,
              "surface_id",
              "T1_9 topological skin (cell data = physical surface id 1..7)")
    return P, T, ids


def export_tet_and_features(vtk_in):
    """Re-export the AlgoHex input split into volume + feature curves."""
    import plot_stages as ps
    P, cells, ty, co = ps.read_vtk(vtk_in)
    tets = [(cells[c], co[c]) for c in range(len(cells)) if ty[c] == 10]
    lines = [(cells[c], co[c]) for c in range(len(cells)) if ty[c] == 3]
    verts = [(cells[c], co[c]) for c in range(len(cells)) if ty[c] == 1]

    print(f"\n[03/04] AlgoHex input {Path(vtk_in).name}")
    print(f"  {len(tets)} tets, {len(lines)} feature edges, "
          f"{len(verts)} feature vertices")
    fc = [c for c, _ in lines] + [c for c, _ in verts]
    ft = [3] * len(lines) + [1] * len(verts)
    fd = [1] * len(lines) + [2] * len(verts)
    write_vtk(OUT / "03_feature_curves.vtk", P, fc, ft, fd, "kind",
              "T1_9 feature graph (1 = feature edge, 2 = feature vertex)")
    write_vtk(OUT / "04_tet_volume.vtk", P, [c for c, _ in tets],
              [10] * len(tets), [0] * len(tets), "zero",
              "T1_9 pure tet mesh fed to AlgoHex")


def export_hexmesh(ovm, out_name, title):
    import ovm_io
    if not Path(ovm).exists():
        return
    P, e, f, poly = ovm_io.read_ovm(ovm)
    hexes, _ = ovm_io.ovm_to_cells(P, e, f, poly)
    vols = np.array([ovm_io._hex_volume(P[c]) for c in hexes])
    print(f"\n[05] {Path(ovm).name}: {len(hexes)} hexes, "
          f"volume {vols.sum():.4f} of 6.5253 "
          f"({100*vols.sum()/6.52525:.1f}%), {(vols<=0).sum()} inverted")
    write_vtk(OUT / out_name, P, hexes, [12] * len(hexes),
              (vols <= 0).astype(int), "inverted", title)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    nodes, elements = parse_msh(MSH)
    export_2d_elements(nodes, elements)
    export_topological(nodes, elements)
    tetvtk = REPO / "data" / "T1_9" / "T1_9_tet.vtk"
    if tetvtk.exists():
        export_tet_and_features(tetvtk)
    export_hexmesh(REPO / "output" / "hex3d_algohex"
                   / "T1_9_hex_v1_classifysurfaces.ovm",
                   "05_hexmesh_v1.vtk",
                   "AlgoHex v1 result (INCOMPLETE, 31.4% of domain)")
    print(f"\nAll VTKs in {OUT}")


if __name__ == "__main__":
    main()
