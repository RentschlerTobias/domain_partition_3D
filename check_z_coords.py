#!/usr/bin/env python3
"""Find hub/shroud nodes by z coordinate."""
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

# Analyze z coordinates
zs = [c[2] for c in nodes.values()]
print(f"Z range: {min(zs):.6f} to {max(zs):.6f}")

# Count nodes at z=0
z0_nodes = [tag for tag, c in nodes.items() if abs(c[2]) < 1e-6]
print(f"Nodes at z=0: {len(z0_nodes)}")

# Count nodes at z=2.5
z25_nodes = [tag for tag, c in nodes.items() if abs(c[2] - 2.5) < 1e-6]
print(f"Nodes at z=2.5: {len(z25_nodes)}")

# Count nodes near z=0 (within 0.1)
near_z0 = [tag for tag, c in nodes.items() if abs(c[2]) < 0.1]
print(f"Nodes near z=0 (within 0.1): {len(near_z0)}")

# Count nodes near z=2.5
near_z25 = [tag for tag, c in nodes.items() if abs(c[2] - 2.5) < 0.1]
print(f"Nodes near z=2.5 (within 0.1): {len(near_z25)}")

# Check z histogram
hist, bins = np.histogram(zs, bins=20)
print("\nZ histogram:")
for i in range(len(hist)):
    print(f"  [{bins[i]:.4f}, {bins[i+1]:.4f}): {hist[i]} nodes")
