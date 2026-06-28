#!/usr/bin/env python3
"""Split T1_9 STL into separate surface files."""
import os

STL_PATH = '/root/repos/block_structured_meshing/T1_9/T1_9_ru_gridGmsh.stl'
OUTPUT_DIR = '/root/repos/block_structured_meshing/T1_9/stl_parts'

# Parse STL
facets = []
with open(STL_PATH, 'r') as f:
    lines = f.readlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('facet normal'):
            normal = [float(x) for x in line.split()[2:5]]
            i += 1
            vertices = []
            for _ in range(3):
                i += 1
                parts = lines[i].strip().split()
                vertices.append([float(parts[1]), float(parts[2]), float(parts[3])])
            avg_z = sum(v[2] for v in vertices) / 3
            facets.append((normal, vertices, avg_z))
        i += 1

print(f"Total facets: {len(facets)}")

# Categorize - order matters to avoid overlap
# 1. Inlet and outlet (most specific)
inlet = [f for f in facets if f[2] < 0.05]
outlet = [f for f in facets if f[2] > 2.45]

# 2. Periodic boundaries (vertical walls)
periodic = [f for f in facets if abs(f[0][2]) < 0.1 and f[2] >= 0.05 and f[2] <= 2.45]

# 3. Hub and shroud (horizontal, but not inlet/outlet)
hub = [f for f in facets if f[0][2] < -0.9 and f[2] > 0.05 and f[2] < 2.45]
shroud = [f for f in facets if f[0][2] > 0.9 and f[2] > 0.05 and f[2] < 2.45]

# 4. Blade is the rest
used = set(id(f) for f in hub + shroud + inlet + outlet + periodic)
blade = [f for f in facets if id(f) not in used]

# Write STL files
def write_stl(facets, filename):
    with open(filename, 'w') as f:
        f.write(f'solid {os.path.basename(filename)}\n')
        for normal, vertices, _ in facets:
            f.write(f'  facet normal {normal[0]} {normal[1]} {normal[2]}\n')
            f.write('    outer loop\n')
            for v in vertices:
                f.write(f'      vertex {v[0]} {v[1]} {v[2]}\n')
            f.write('    endloop\n')
            f.write('  endfacet\n')
        f.write(f'endsolid {os.path.basename(filename)}\n')

os.makedirs(OUTPUT_DIR, exist_ok=True)
write_stl(hub, f'{OUTPUT_DIR}/hub.stl')
write_stl(shroud, f'{OUTPUT_DIR}/shroud.stl')
write_stl(blade, f'{OUTPUT_DIR}/blade.stl')
write_stl(inlet, f'{OUTPUT_DIR}/inlet.stl')
write_stl(outlet, f'{OUTPUT_DIR}/outlet.stl')
write_stl(periodic, f'{OUTPUT_DIR}/periodic.stl')

print(f"Hub: {len(hub)} facets")
print(f"Shroud: {len(shroud)} facets")
print(f"Blade: {len(blade)} facets")
print(f"Inlet: {len(inlet)} facets")
print(f"Outlet: {len(outlet)} facets")
print(f"Periodic: {len(periodic)} facets")
print(f"Total check: {len(hub) + len(shroud) + len(blade) + len(inlet) + len(outlet) + len(periodic)}")
print(f"\nFiles written to {OUTPUT_DIR}")
