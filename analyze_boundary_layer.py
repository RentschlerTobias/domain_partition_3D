#!/usr/bin/env python3
"""
Analyze the boundary-layer structure of T1_9_ru_gridGmsh.msh.

The boundary layer is formed by tetrahedra and wedges connecting the hub (z=0)
and shroud (z=2.5) surfaces to the interior hexahedral core.
"""
import numpy as np
from collections import defaultdict, Counter
import sys

MESH_PATH = '/home/t1dde/Duty/projects/domain_partition/domain_partition_3D/T1_9/T1_9_ru_gridGmsh.msh'

# Element type names (Gmsh 2.2 format)
ETYPE_NAMES = {
    1: '2-node line',
    2: '3-node triangle',
    3: '4-node quad',
    4: '4-node tetrahedron',
    5: '8-node hexahedron',
    6: '6-node prism (wedge)',
    7: '5-node pyramid',
    15: '1-node point',
}


def parse_msh(path):
    """Parse Gmsh 2.2 ASCII format file."""
    print(f"Parsing {path}...")
    
    physical_names = {}
    nodes = {}
    elements = []
    
    with open(path, 'r') as f:
        section = None
        for line in f:
            line = line.strip()
            if line.startswith('$'):
                if line == '$PhysicalNames':
                    section = 'physical_names'
                elif line == '$Nodes':
                    section = 'nodes'
                elif line == '$Elements':
                    section = 'elements'
                elif line == '$EndPhysicalNames':
                    section = None
                elif line == '$EndNodes':
                    section = None
                elif line == '$EndElements':
                    section = None
                continue
            
            if section == 'physical_names':
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        dim = int(parts[0])
                        tag = int(parts[1])
                        name = parts[2].strip('"')
                        if dim not in physical_names:
                            physical_names[dim] = {}
                        physical_names[dim][tag] = name
                    except ValueError:
                        pass
            
            elif section == 'nodes':
                parts = line.split()
                if len(parts) == 4:
                    try:
                        tag = int(parts[0])
                        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                        nodes[tag] = (x, y, z)
                    except ValueError:
                        pass
            
            elif section == 'elements':
                parts = line.split()
                if len(parts) >= 5:
                    try:
                        elem_id = int(parts[0])
                        elem_type = int(parts[1])
                        n_tags = int(parts[2])
                        tags = [int(t) for t in parts[3:3+n_tags]]
                        node_list = [int(n) for n in parts[3+n_tags:]]
                        phys_tag = tags[0] if len(tags) > 0 else 0
                        geom_tag = tags[1] if len(tags) > 1 else 0
                        elements.append((elem_id, elem_type, phys_tag, geom_tag, node_list))
                    except (ValueError, IndexError):
                        pass
    
    print(f"  Nodes: {len(nodes)}")
    print(f"  Elements: {len(elements)}")
    
    return physical_names, nodes, elements


def get_boundary_faces_for_element(elem_type, node_list, nodes, z_tol=1e-6):
    """
    Return list of faces that lie at z=0 or z=2.5.
    Each face is a tuple (z_level, node_list).
    """
    faces = []
    
    if elem_type == 4:  # Tetrahedron - 4 triangular faces
        tet_faces = [(0,1,2), (0,1,3), (0,2,3), (1,2,3)]
        for face_idx in tet_faces:
            face_nodes = [node_list[i] for i in face_idx]
            zs = [nodes[n][2] for n in face_nodes if n in nodes]
            if len(zs) == 3:
                if all(abs(z) < z_tol for z in zs):
                    faces.append(('z0', face_nodes))
                elif all(abs(z - 2.5) < z_tol for z in zs):
                    faces.append(('z25', face_nodes))
    
    elif elem_type == 5:  # Hexahedron - 6 quad faces
        hex_faces = [
            (0,1,2,3),  # bottom
            (4,5,6,7),  # top
            (0,1,5,4),  # front
            (1,2,6,5),  # right
            (2,3,7,6),  # back
            (3,0,4,7),  # left
        ]
        for face_idx in hex_faces:
            face_nodes = [node_list[i] for i in face_idx]
            zs = [nodes[n][2] for n in face_nodes if n in nodes]
            if len(zs) == 4:
                if all(abs(z) < z_tol for z in zs):
                    faces.append(('z0', face_nodes))
                elif all(abs(z - 2.5) < z_tol for z in zs):
                    faces.append(('z25', face_nodes))
    
    elif elem_type == 6:  # Wedge - 2 triangular + 3 quad faces
        wedge_faces = [
            (0,1,2),     # bottom triangle
            (3,4,5),     # top triangle
            (0,1,4,3),  # quad 1
            (1,2,5,4),  # quad 2
            (2,0,3,5),  # quad 3
        ]
        for face_idx in wedge_faces:
            face_nodes = [node_list[i] for i in face_idx]
            zs = [nodes[n][2] for n in face_nodes if n in nodes]
            if len(zs) == len(face_idx):
                if all(abs(z) < z_tol for z in zs):
                    faces.append(('z0', face_nodes))
                elif all(abs(z - 2.5) < z_tol for z in zs):
                    faces.append(('z25', face_nodes))
    
    elif elem_type == 7:  # Pyramid - 1 quad base + 4 triangular sides
        pyr_faces = [
            (0,1,2,3),  # base (quad)
            (0,1,4),    # side 1
            (1,2,4),    # side 2
            (2,3,4),    # side 3
            (3,0,4),    # side 4
        ]
        for face_idx in pyr_faces:
            face_nodes = [node_list[i] for i in face_idx]
            zs = [nodes[n][2] for n in face_nodes if n in nodes]
            if len(zs) == len(face_idx):
                if all(abs(z) < z_tol for z in zs):
                    faces.append(('z0', face_nodes))
                elif all(abs(z - 2.5) < z_tol for z in zs):
                    faces.append(('z25', face_nodes))
    
    return faces


def main():
    physical_names, nodes, elements = parse_msh(MESH_PATH)
    
    # Print physical names
    print("\n=== Physical Names ===")
    print("2D surfaces:")
    for tag, name in sorted(physical_names.get(2, {}).items()):
        print(f"  {tag}: {name}")
    print("3D volumes:")
    for tag, name in sorted(physical_names.get(3, {}).items()):
        print(f"  {tag}: {name}")
    
    # Separate elements by type
    elems_by_type = defaultdict(list)
    for elem_id, elem_type, phys_tag, geom_tag, node_list in elements:
        elems_by_type[elem_type].append((elem_id, phys_tag, geom_tag, node_list))
    
    print("\n=== Element Type Summary ===")
    for et in sorted(elems_by_type.keys()):
        print(f"  Type {et} ({ETYPE_NAMES.get(et, 'unknown')}): {len(elems_by_type[et])}")
    
    # Find boundary faces at z=0 and z=2.5
    print("\n=== Boundary Face Collection ===")
    
    z0_faces = []  # (elem_id, elem_type, phys_tag, face_nodes)
    z25_faces = []
    
    for elem_id, elem_type, phys_tag, geom_tag, node_list in elements:
        if elem_type not in (4, 5, 6, 7):
            continue
        faces = get_boundary_faces_for_element(elem_type, node_list, nodes)
        for level, face_nodes in faces:
            if level == 'z0':
                z0_faces.append((elem_id, elem_type, phys_tag, face_nodes))
            elif level == 'z25':
                z25_faces.append((elem_id, elem_type, phys_tag, face_nodes))
    
    print(f"  Total faces at z=0: {len(z0_faces)}")
    print(f"  Total faces at z=2.5: {len(z25_faces)}")
    
    # Breakdown by element type
    print("\n  z=0 faces by element type:")
    z0_by_type = Counter(f[1] for f in z0_faces)
    for et, cnt in sorted(z0_by_type.items()):
        print(f"    Type {et} ({ETYPE_NAMES.get(et, 'unknown')}): {cnt}")
    
    print("\n  z=2.5 faces by element type:")
    z25_by_type = Counter(f[1] for f in z25_faces)
    for et, cnt in sorted(z25_by_type.items()):
        print(f"    Type {et} ({ETYPE_NAMES.get(et, 'unknown')}): {cnt}")
    
    # Group boundary faces into patches using union-find
    print("\n=== Boundary Patch Grouping ===")
    
    def group_faces_into_patches(faces):
        """Group faces that share nodes into patches."""
        parent = {}
        
        def find(x):
            while parent.get(x, x) != x:
                parent[x] = parent.get(parent[x], parent[x])
                x = parent[x]
            return x
        
        def union(x, y):
            rx, ry = find(x), find(y)
            if rx != ry:
                parent[rx] = ry
        
        for face in faces:
            face_nodes = face[3]
            for n in face_nodes:
                if n not in parent:
                    parent[n] = n
            for i in range(1, len(face_nodes)):
                union(face_nodes[0], face_nodes[i])
        
        groups = defaultdict(list)
        for face in faces:
            root = find(face[3][0])
            groups[root].append(face)
        
        return groups
    
    z0_groups = group_faces_into_patches(z0_faces)
    z25_groups = group_faces_into_patches(z25_faces)
    
    print(f"\n  Distinct patches at z=0: {len(z0_groups)}")
    print(f"  Distinct patches at z=2.5: {len(z25_groups)}")
    
    # Analyze each patch
    def analyze_patch(groups, label, nodes):
        print(f"\n  --- {label} Patches (sorted by size) ---")
        patch_list = []
        for i, (root, faces) in enumerate(sorted(groups.items(), key=lambda x: -len(x[1]))):
            all_nodes = set()
            for f in faces:
                all_nodes.update(f[3])
            
            # Count face types
            face_types = Counter(f[1] for f in faces)
            
            # Count face shapes (tri vs quad)
            n_tri_faces = sum(1 for f in faces if len(f[3]) == 3)
            n_quad_faces = sum(1 for f in faces if len(f[3]) == 4)
            
            # Compute centroid
            coords = np.array([nodes[n] for n in all_nodes if n in nodes])
            if len(coords) == 0:
                continue
            centroid = coords.mean(axis=0)
            r = np.sqrt(centroid[0]**2 + centroid[1]**2)
            theta = np.degrees(np.arctan2(centroid[1], centroid[0]))
            
            patch_list.append({
                'id': i,
                'n_nodes': len(all_nodes),
                'n_faces': len(faces),
                'n_tri_faces': n_tri_faces,
                'n_quad_faces': n_quad_faces,
                'face_types': dict(face_types),
                'centroid': centroid,
                'r': r,
                'theta': theta,
            })
            
            print(f"    Patch {i}: {len(all_nodes)} nodes, {len(faces)} faces "
                  f"(tri={n_tri_faces}, quad={n_quad_faces}), "
                  f"types={dict(face_types)}, r={r:.3f}, theta={theta:.1f}°")
        
        return patch_list
    
    patch_info_z0 = analyze_patch(z0_groups, "z=0", nodes)
    patch_info_z25 = analyze_patch(z25_groups, "z=2.5", nodes)
    
    # Find corner nodes: nodes with high valence in boundary faces
    print("\n=== Corner Node Analysis ===")
    
    # Count how many boundary faces each node belongs to
    node_face_count = defaultdict(int)
    for face in z0_faces + z25_faces:
        for n in face[3]:
            node_face_count[n] += 1
    
    # Find nodes at z=0 and z=2.5
    z0_nodes = set()
    z25_nodes = set()
    for tag, (x, y, z) in nodes.items():
        if abs(z) < 1e-6:
            z0_nodes.add(tag)
        elif abs(z - 2.5) < 1e-6:
            z25_nodes.add(tag)
    
    print(f"\n  Nodes at z=0: {len(z0_nodes)}")
    print(f"  Nodes at z=2.5: {len(z25_nodes)}")
    
    # Valence distribution
    print("\n  Valence distribution at z=0:")
    z0_valence = Counter()
    for n in z0_nodes:
        z0_valence[node_face_count[n]] += 1
    for v in sorted(z0_valence.keys()):
        print(f"    Valence {v}: {z0_valence[v]} nodes")
    
    print("\n  Valence distribution at z=2.5:")
    z25_valence = Counter()
    for n in z25_nodes:
        z25_valence[node_face_count[n]] += 1
    for v in sorted(z25_valence.keys()):
        print(f"    Valence {v}: {z25_valence[v]} nodes")
    
    # High-valence nodes (corners)
    print("\n  High-valence nodes at z=0 (valence >= 3):")
    high_val_z0 = [(n, node_face_count[n]) for n in z0_nodes if node_face_count[n] >= 3]
    high_val_z0.sort(key=lambda x: -x[1])
    print(f"    Count: {len(high_val_z0)}")
    for n, v in high_val_z0[:30]:
        if n in nodes:
            x, y, z = nodes[n]
            r = np.sqrt(x**2 + y**2)
            theta = np.degrees(np.arctan2(y, x))
            print(f"      Node {n}: valence={v}, r={r:.3f}, theta={theta:.1f}°, pos=({x:.3f}, {y:.3f})")
    
    print("\n  High-valence nodes at z=2.5 (valence >= 3):")
    high_val_z25 = [(n, node_face_count[n]) for n in z25_nodes if node_face_count[n] >= 3]
    high_val_z25.sort(key=lambda x: -x[1])
    print(f"    Count: {len(high_val_z25)}")
    for n, v in high_val_z25[:30]:
        if n in nodes:
            x, y, z = nodes[n]
            r = np.sqrt(x**2 + y**2)
            theta = np.degrees(np.arctan2(y, x))
            print(f"      Node {n}: valence={v}, r={r:.3f}, theta={theta:.1f}°, pos=({x:.3f}, {y:.3f})")
    
    # Analyze how hex/wedge elements connect boundary to interior
    print("\n=== Boundary-to-Interior Connection Analysis ===")
    
    # For each boundary face, find the 3D element it belongs to
    # and trace how that element connects to the interior
    
    # Count elements that have a boundary face
    elems_with_z0_face = set()
    elems_with_z25_face = set()
    for elem_id, elem_type, phys_tag, face_nodes in z0_faces:
        elems_with_z0_face.add(elem_id)
    for elem_id, elem_type, phys_tag, face_nodes in z25_faces:
        elems_with_z25_face.add(elem_id)
    
    print(f"\n  Elements with z=0 face: {len(elems_with_z0_face)}")
    print(f"  Elements with z=2.5 face: {len(elems_with_z25_face)}")
    
    # For these elements, find their z-range
    elem_z_range = {}
    for elem_id, elem_type, phys_tag, geom_tag, node_list in elements:
        zs = [nodes[n][2] for n in node_list if n in nodes]
        if zs:
            elem_z_range[elem_id] = (min(zs), max(zs), elem_type)
    
    # Distribution of z-thickness for elements touching boundary
    print("\n  Z-thickness of elements touching z=0:")
    thicknesses_z0 = []
    for eid in elems_with_z0_face:
        if eid in elem_z_range:
            zmin, zmax, et = elem_z_range[eid]
            thicknesses_z0.append((zmax - zmin, et))
    
    thick_by_type_z0 = defaultdict(list)
    for t, et in thicknesses_z0:
        thick_by_type_z0[et].append(t)
    for et in sorted(thick_by_type_z0.keys()):
        ts = thick_by_type_z0[et]
        print(f"    Type {et}: {len(ts)} elements, thickness range [{min(ts):.4f}, {max(ts):.4f}], mean={np.mean(ts):.4f}")
    
    print("\n  Z-thickness of elements touching z=2.5:")
    thicknesses_z25 = []
    for eid in elems_with_z25_face:
        if eid in elem_z_range:
            zmin, zmax, et = elem_z_range[eid]
            thicknesses_z25.append((zmax - zmin, et))
    
    thick_by_type_z25 = defaultdict(list)
    for t, et in thicknesses_z25:
        thick_by_type_z25[et].append(t)
    for et in sorted(thick_by_type_z25.keys()):
        ts = thick_by_type_z25[et]
        print(f"    Type {et}: {len(ts)} elements, thickness range [{min(ts):.4f}, {max(ts):.4f}], mean={np.mean(ts):.4f}")
    
    # Analyze the interior hex layer
    print("\n=== Interior Hex Layer Analysis ===")
    
    hex_elems = elems_by_type.get(5, [])
    hex_z_pairs = []
    for elem_id, phys_tag, geom_tag, node_list in hex_elems:
        zs = [nodes[n][2] for n in node_list if n in nodes]
        if len(zs) == 8:
            hex_z_pairs.append((min(zs), max(zs)))
    
    # Histogram of hex z-centers
    hex_z_centers = [(zmin + zmax) / 2 for zmin, zmax in hex_z_pairs]
    print(f"\n  Hex z-center range: [{min(hex_z_centers):.4f}, {max(hex_z_centers):.4f}]")
    print(f"  Hex z-min range: [{min(z[0] for z in hex_z_pairs):.4f}, {max(z[0] for z in hex_z_pairs):.4f}]")
    print(f"  Hex z-max range: [{min(z[1] for z in hex_z_pairs):.4f}, {max(z[1] for z in hex_z_pairs):.4f}]")
    
    # Histogram of hex z-centers
    hist, bin_edges = np.histogram(hex_z_centers, bins=20)
    print("\n  Hex z-center histogram:")
    for i in range(len(hist)):
        if hist[i] > 0:
            print(f"    [{bin_edges[i]:.4f}, {bin_edges[i+1]:.4f}): {hist[i]} hex")
    
    # Pyramid analysis
    print("\n=== Pyramid (Transition) Element Analysis ===")
    pyr_elems = elems_by_type.get(7, [])
    pyr_z_pairs = []
    for elem_id, phys_tag, geom_tag, node_list in pyr_elems:
        zs = [nodes[n][2] for n in node_list if n in nodes]
        if len(zs) == 5:
            pyr_z_pairs.append((min(zs), max(zs)))
    
    if pyr_z_pairs:
        pyr_z_centers = [(zmin + zmax) / 2 for zmin, zmax in pyr_z_pairs]
        print(f"\n  Pyramid count: {len(pyr_elems)}")
        print(f"  Pyramid z-center range: [{min(pyr_z_centers):.4f}, {max(pyr_z_centers):.4f}]")
        print(f"  Pyramid z-min range: [{min(z[0] for z in pyr_z_pairs):.4f}, {max(z[0] for z in pyr_z_pairs):.4f}]")
        print(f"  Pyramid z-max range: [{min(z[1] for z in pyr_z_pairs):.4f}, {max(z[1] for z in pyr_z_pairs):.4f}]")
    
    # Summary
    print("\n" + "="*60)
    print("=== SUMMARY ===")
    print("="*60)
    print(f"Total nodes: {len(nodes)}")
    print(f"Total elements: {len(elements)}")
    print(f"  Tet (type 4): {len(elems_by_type.get(4, []))}")
    print(f"  Hex (type 5): {len(elems_by_type.get(5, []))}")
    print(f"  Wedge (type 6): {len(elems_by_type.get(6, []))}")
    print(f"  Pyramid (type 7): {len(elems_by_type.get(7, []))}")
    print(f"\nBoundary nodes at z=0: {len(z0_nodes)}")
    print(f"Boundary nodes at z=2.5: {len(z25_nodes)}")
    print(f"\nBoundary patches at z=0: {len(patch_info_z0)}")
    print(f"Boundary patches at z=2.5: {len(patch_info_z25)}")
    print(f"\nHigh-valence nodes at z=0: {len(high_val_z0)}")
    print(f"High-valence nodes at z=2.5: {len(high_val_z25)}")


if __name__ == '__main__':
    main()
