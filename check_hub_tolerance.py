#!/usr/bin/env python3
"""Extract hub surface with relaxed tolerance."""
import numpy as np

MESH_PATH = '/root/repos/block_structured_meshing/T1_9/T1_9_ru_gridGmsh.msh'

# Parse nodes
nodes = {}
with open(MESH_PATH, 'r') as f:
    lines = f.readlines()
    in_nodes = False
    for line in lines:
        if '$Nodes' in line:
            in_nodes = True
            continue
        if '$EndNodes' in line:
            break
        if in_nodes:
            parts = line.strip().split()
            if len(parts) == 4:
                tag = int(parts[0])
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                nodes[tag] = (x, y, z)

# Parse all 2D elements
faces = []
with open(MESH_PATH, 'r') as f:
    lines = f.readlines()
    in_elements = False
    for line in lines:
        if '$Elements' in line:
            in_elements = True
            continue
        if '$EndElements' in line:
            break
        if in_elements:
            parts = line.strip().split()
            if len(parts) == 1:
                continue
            if len(parts) < 5:
                continue
            elem_type = int(parts[1])
            node_list = [int(p) for p in parts[4:]]
            
            if elem_type == 2:  # Triangle
                faces.append(('tri', node_list))
            elif elem_type == 3:  # Quad
                faces.append(('quad', node_list))

# Check z values for all face nodes
face_z_stats = []
for ftype, node_list in faces:
    zs = [nodes[n][2] for n in node_list if n in nodes]
    face_z_stats.append((ftype, min(zs), max(zs), sum(zs)/len(zs)))

# Count faces by z range
print("=== Faces by z range ===")
for tol in [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1]:
    hub_count = sum(1 for ftype, zmin, zmax, zavg in face_z_stats if zmax < tol)
    shroud_count = sum(1 for ftype, zmin, zmax, zavg in face_z_stats if zmin > 2.5 - tol)
    print(f"  tol={tol}: hub={hub_count}, shroud={shroud_count}")

# Check if any faces have zmin close to 0
print("\n=== Faces with zmin < 0.01 ===")
for ftype, zmin, zmax, zavg in face_z_stats:
    if zmin < 0.01:
        print(f"  {ftype}: zmin={zmin:.6f}, zmax={zmax:.6f}, zavg={zavg:.6f}")
        break

# Check all face z ranges
z_ranges = sorted(set((round(zmin, 4), round(zmax, 4)) for ftype, zmin, zmax, zavg in face_z_stats))
print(f"\nUnique face z ranges: {len(z_ranges)}")
for zmin, zmax in z_ranges[:10]:
    count = sum(1 for ft, zm, zM, za in face_z_stats if abs(zm - zmin) < 1e-6 and abs(zM - zmax) < 1e-6)
    print(f"  [{zmin:.4f}, {zmax:.4f}]: {count} faces")
