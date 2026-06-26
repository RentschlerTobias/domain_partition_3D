#!/usr/bin/env python3
"""
Extract surfaces with reg_geom=1 and reg_geom=2 from T1_9 MSH file
and write them as STL files.

MSH format 2.2:
  element_num element_type num_tags reg_phys reg_geom ...

2D elements:
  type 2 = triangle (3 nodes)
  type 3 = quadrilateral (4 nodes)

reg_geom=1 -> hub
reg_geom=2 -> shroud
"""
import struct
import os

INPUT = "/root/repos/block_structured_meshing/T1_9/T1_9_ru_gridGmsh.msh"
OUT_DIR = "/root/repos/block_structured_meshing"

def read_msh_nodes(lines):
    """Read nodes section and return dict {tag: (x,y,z)}."""
    nodes = {}
    in_nodes = False
    for line in lines:
        line = line.strip()
        if line == "$Nodes":
            in_nodes = True
            continue
        if line == "$EndNodes":
            break
        if not in_nodes:
            continue
        parts = line.split()
        if len(parts) == 4:
            tag = int(parts[0])
            x, y, z = map(float, parts[1:])
            nodes[tag] = (x, y, z)
    return nodes

def read_msh_2d_elements(lines, target_reg_geom):
    """
    Read elements section and extract 2D elements (triangles=2, quads=3)
    matching target_reg_geom.
    Returns list of node index lists.
    """
    elements = []
    in_elements = False
    for line in lines:
        line = line.strip()
        if line == "$Elements":
            in_elements = True
            continue
        if line == "$EndElements":
            break
        if not in_elements:
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        # Format: elem_num elem_type num_tags reg_phys reg_geom [nodes...]
        elem_num = int(parts[0])
        elem_type = int(parts[1])
        num_tags = int(parts[2])
        reg_phys = int(parts[3])
        reg_geom = int(parts[4])
        
        if reg_geom != target_reg_geom:
            continue
        
        if elem_type == 2:  # Triangle
            # 3 nodes after 2 tags (reg_phys, reg_geom)
            nodes = list(map(int, parts[5:8]))
            elements.append(nodes)
        elif elem_type == 3:  # Quad
            # 4 nodes after 2 tags
            nodes = list(map(int, parts[5:9]))
            # Triangulate quad: (0,1,2) and (0,2,3)
            elements.append([nodes[0], nodes[1], nodes[2]])
            elements.append([nodes[0], nodes[2], nodes[3]])
    return elements

def compute_normal(p1, p2, p3):
    """Compute unit normal for triangle."""
    ux = p2[0] - p1[0]
    uy = p2[1] - p1[1]
    uz = p2[2] - p1[2]
    vx = p3[0] - p1[0]
    vy = p3[1] - p1[1]
    vz = p3[2] - p1[2]
    nx = uy * vz - uz * vy
    ny = uz * vx - ux * vz
    nz = ux * vy - uy * vx
    length = (nx**2 + ny**2 + nz**2) ** 0.5
    if length == 0:
        return (0.0, 0.0, 0.0)
    return (nx/length, ny/length, nz/length)

def write_binary_stl(elements, nodes, filename):
    """Write binary STL file."""
    header = b"dtOO surface extraction" + b" " * (80 - 23)
    num_triangles = len(elements)
    with open(filename, "wb") as f:
        f.write(header)
        f.write(struct.pack("<I", num_triangles))
        for tri in elements:
            p1 = nodes[tri[0]]
            p2 = nodes[tri[1]]
            p3 = nodes[tri[2]]
            normal = compute_normal(p1, p2, p3)
            f.write(struct.pack("<fff", *normal))
            f.write(struct.pack("<fff", *p1))
            f.write(struct.pack("<fff", *p2))
            f.write(struct.pack("<fff", *p3))
            f.write(struct.pack("<H", 0))  # attribute byte count
    print(f"Wrote {num_triangles} triangles to {filename}")

def main():
    print("Reading MSH file...")
    with open(INPUT, "r") as f:
        lines = f.readlines()
    
    nodes = read_msh_nodes(lines)
    print(f"Read {len(nodes)} nodes")
    
    hub_elements = read_msh_2d_elements(lines, target_reg_geom=1)
    print(f"Hub (reg_geom=1): {len(hub_elements)} triangles")
    
    shroud_elements = read_msh_2d_elements(lines, target_reg_geom=2)
    print(f"Shroud (reg_geom=2): {len(shroud_elements)} triangles")
    
    write_binary_stl(hub_elements, nodes, os.path.join(OUT_DIR, "T1_9_hub_raw.stl"))
    write_binary_stl(shroud_elements, nodes, os.path.join(OUT_DIR, "T1_9_shroud_raw.stl"))
    
    print("Done.")

if __name__ == "__main__":
    main()
