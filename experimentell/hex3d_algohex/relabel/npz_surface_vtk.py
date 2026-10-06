"""Write the labelled boundary surface of a sample.npz as AlgoHex input VTK,
in the format clean_blocks.read_input_surface reads (ASCII, CELL_DATA color).
Replaces tet.vtk for dataset samples that only have sample.npz + blocks.vtk.

    python npz_surface_vtk.py <sample.npz> <surface.vtk>
"""
import sys
import numpy as np

s = np.load(sys.argv[1])
P, T, L = s["surface_points"], s["surface_tris"], s["surface_tri_label"]
with open(sys.argv[2], "w") as f:
    f.write("# vtk DataFile Version 2.0\nsurface from sample.npz\nASCII\nDATASET UNSTRUCTURED_GRID\n")
    f.write(f"POINTS {len(P)} double\n")
    f.writelines(f"{x:.12g} {y:.12g} {z:.12g}\n" for x, y, z in P)
    f.write(f"CELLS {len(T)} {4 * len(T)}\n")
    f.writelines(f"3 {a} {b} {c}\n" for a, b, c in T)
    f.write(f"CELL_TYPES {len(T)}\n" + "5\n" * len(T))
    f.write(f"CELL_DATA {len(T)}\nSCALARS color int 1\nLOOKUP_TABLE default\n")
    f.writelines(f"{int(l)}\n" for l in L)
