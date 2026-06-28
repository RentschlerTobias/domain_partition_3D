#!/usr/bin/env python3
"""Extract hub surface from STL file."""
import numpy as np

STL_PATH = '/root/repos/block_structured_meshing/T1_9/T1_9_ru_gridGmsh.stl'

# Parse STL
facets = []
with open(STL_PATH, 'r') as f:
    lines = f.readlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('facet normal'):
            # Read normal
            parts = line.split()
            normal = [float(parts[2]), float(parts[3]), float(parts[4])]
            
            # Read vertices
            i += 1
            if 'outer loop' not in lines[i].strip():
                i += 1
                continue
            
            vertices = []
            for j in range(3):
                i += 1
                vline = lines[i].strip()
                parts = vline.split()
                if parts[0] == 'vertex':
                    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                    vertices.append((x, y, z))
            
            # Check if all vertices are at z=0
            all_z0 = all(abs(v[2]) < 1e-6 for v in vertices)
            if all_z0:
                facets.append((normal, vertices))
        i += 1

print(f"Total hub facets: {len(facets)}")

# Write hub STL
HUB_STL_PATH = '/root/repos/block_structured_meshing/T1_9_hub.stl'
with open(HUB_STL_PATH, 'w') as f:
    f.write("solid hub\n")
    for normal, vertices in facets:
        f.write(f"  facet normal {normal[0]} {normal[1]} {normal[2]}\n")
        f.write("    outer loop\n")
        for v in vertices:
            f.write(f"      vertex {v[0]} {v[1]} {v[2]}\n")
        f.write("    endloop\n")
        f.write("  endfacet\n")
    f.write("endsolid hub\n")

print(f"Wrote hub STL to {HUB_STL_PATH}")

# Check if hub is flat
if facets:
    zs = [v[2] for _, vertices in facets for v in vertices]
    print(f"Hub z range: {min(zs):.6f} to {max(zs):.6f}")
    
    # Check normals
    normals = [n for n, _ in facets]
    avg_normal = np.mean(normals, axis=0)
    print(f"Average normal: {avg_normal}")

# Also check for shroud facets
shroud_facets = []
with open(STL_PATH, 'r') as f:
    lines = f.readlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('facet normal'):
            parts = line.split()
            normal = [float(parts[2]), float(parts[3]), float(parts[4])]
            
            i += 1
            if 'outer loop' not in lines[i].strip():
                i += 1
                continue
            
            vertices = []
            for j in range(3):
                i += 1
                vline = lines[i].strip()
                parts = vline.split()
                if parts[0] == 'vertex':
                    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                    vertices.append((x, y, z))
            
            all_z25 = all(abs(v[2] - 2.5) < 1e-6 for v in vertices)
            if all_z25:
                shroud_facets.append((normal, vertices))
        i += 1

print(f"\nTotal shroud facets: {len(shroud_facets)}")

# Write shroud STL
SHROUD_STL_PATH = '/root/repos/block_structured_meshing/T1_9_shroud.stl'
with open(SHROUD_STL_PATH, 'w') as f:
    f.write("solid shroud\n")
    for normal, vertices in shroud_facets:
        f.write(f"  facet normal {normal[0]} {normal[1]} {normal[2]}\n")
        f.write("    outer loop\n")
        for v in vertices:
            f.write(f"      vertex {v[0]} {v[1]} {v[2]}\n")
        f.write("    endloop\n")
        f.write("  endfacet\n")
    f.write("endsolid shroud\n")

print(f"Wrote shroud STL to {SHROUD_STL_PATH}")
