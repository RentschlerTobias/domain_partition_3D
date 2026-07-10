#!/usr/bin/env python3
"""
Erstelle VTK-Datei mit den dtOO Hexa-Bloecken.
Exportiert die 28 Block-Instanzen als Hexaeder (Voxels).
"""

import json
import numpy as np
from pathlib import Path

BASE = Path("/home/t1dde/Duty/projects/domain_partition/domain_partition_3D")
blocks_json = BASE / "block_instances_28.json"
out_vtk = BASE / "dtOO_hexa_blocks.vtk"

# Laden
with open(blocks_json, "r") as f:
    blocks_data = json.load(f)

instances = blocks_data["instances"]

# Sammle alle einzigartigen Knoten und Hexaeder
all_nodes = {}  # node_id -> (x, y, z)
hexahedra = []   # Liste von 8 node indices
node_counter = 0

for inst in instances:
    corners = inst["corners"]
    
    # Pruefe ob wir 8 einzigartige Ecken haben
    # (in der Praxis sollten die 8 Ecken eines Hexaeder sein)
    node_indices = []
    for c in corners:
        key = tuple(np.round(c, 6))  # Runden fuer Gleichheit
        if key not in all_nodes:
            all_nodes[key] = node_counter
            node_counter += 1
        node_indices.append(all_nodes[key])
    
    # Die 8 Indizes bilden ein Hexaeder
    # dtOO Reihenfolge: Bottom face (0,1,2,3), Top face (4,5,6,7)
    # VTK Hexaeder Reihenfolge: (0,1,2,3,4,5,6,7)
    if len(node_indices) == 8:
        hexahedra.append(node_indices)

# Konvertiere zu Arrays
node_coords = np.zeros((len(all_nodes), 3))
for key, idx in all_nodes.items():
    node_coords[idx] = key

n_nodes = len(all_nodes)
n_cells = len(hexahedra)

# VTK Legacy ASCII Format schreiben
with open(out_vtk, "w") as f:
    f.write("# vtk DataFile Version 3.0\n")
    f.write("dtOO Hexa Blocks\n")
    f.write("ASCII\n")
    f.write("DATASET UNSTRUCTURED_GRID\n")
    
    # Punkte
    f.write(f"POINTS {n_nodes} double\n")
    for coord in node_coords:
        f.write(f"{coord[0]:.10f} {coord[1]:.10f} {coord[2]:.10f}\n")
    
    # Zellen
    f.write(f"\nCELLS {n_cells} {n_cells * 9}\n")
    for hex_nodes in hexahedra:
        f.write(f"8 {' '.join(map(str, hex_nodes))}\n")
    
    # Zelltypen (12 = Hexahedron)
    f.write(f"\nCELL_TYPES {n_cells}\n")
    for _ in hexahedra:
        f.write("12\n")
    
    # Zell-Daten: Block-Name und Region
    f.write(f"\nCELL_DATA {n_cells}\n")
    f.write("SCALARS block_id int 1\n")
    f.write("LOOKUP_TABLE default\n")
    for i, inst in enumerate(instances):
        f.write(f"{i}\n")
    
    f.write("\nSCALARS region int 1\n")
    f.write("LOOKUP_TABLE default\n")
    for inst in instances:
        # Extrahiere Rotations-Index aus dem Namen
        name = inst["name"]
        rot = inst["rotation"]
        f.write(f"{rot}\n")

print(f"VTK gespeichert: {out_vtk}")
print(f"  Knoten: {n_nodes}")
print(f"  Hexaeder: {n_cells}")
print(f"  Bloecke: {len(instances)}")
