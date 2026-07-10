#!/usr/bin/env python3
"""
Analyze the hex core structure of T1_9/T1_9_ru_gridGmsh.msh.

Key finding: hex elements have phys_tag=0 (unassigned) but geom_tag=2..6
which identifies the 5 volume regions. The physical names R_0..R_5 are
defined but no elements reference them.

Goals:
1. Parse all type-5 (hex) elements from the .msh file
2. Group them by geometry tag (the actual volume assignment)
3. For each region, compute the bounding box and element count
4. Analyze connectivity: do hex elements form regular structured blocks?
5. Look for O-grid pattern: is there a central block surrounded by peripheral blocks?
6. Extract a single z-layer of hex elements and visualize their 2D arrangement
7. Find singular nodes (valence > 4 in hex connectivity graph)
"""
import numpy as np
from collections import defaultdict, Counter
import json
import os
import time

MESH_PATH = '/home/t1dde/Duty/projects/domain_partition/domain_partition_3D/T1_9/T1_9_ru_gridGmsh.msh'
EVIDENCE_DIR = '/home/t1dde/Duty/projects/domain_partition/domain_partition_3D/.omo/evidence'
OUTPUT_DIR = '/home/t1dde/Duty/projects/domain_partition/domain_partition_3D/output'

os.makedirs(EVIDENCE_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)


def parse_msh_fast(path):
    """Fast parser: extract nodes and hex elements only."""
    print(f"Parsing {path}...")
    t0 = time.time()

    physical_names = {}
    nodes = {}
    hex_elems = []

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
                if line and line[0].isdigit():
                    parts = line.split()
                    if len(parts) >= 5:
                        try:
                            elem_type = int(parts[1])
                            if elem_type == 5:
                                elem_id = int(parts[0])
                                n_tags = int(parts[2])
                                tags = [int(t) for t in parts[3:3+n_tags]]
                                node_list = [int(n) for n in parts[3+n_tags:]]
                                phys_tag = tags[0] if len(tags) > 0 else 0
                                geom_tag = tags[1] if len(tags) > 1 else 0
                                hex_elems.append((elem_id, phys_tag, geom_tag, node_list))
                        except (ValueError, IndexError):
                            pass

    print(f"  Parsed in {time.time()-t0:.2f}s")
    print(f"  Nodes: {len(nodes)}")
    print(f"  Hex elements: {len(hex_elems)}")
    return physical_names, nodes, hex_elems


def main():
    physical_names, nodes, hex_elems = parse_msh_fast(MESH_PATH)

    print("\n=== Physical Names (3D volumes) ===")
    for tag, name in sorted(physical_names.get(3, {}).items()):
        print(f"  {tag}: {name}")

    print("\n=== Tag Distribution ===")
    phys_tags = Counter(e[1] for e in hex_elems)
    geom_tags = Counter(e[2] for e in hex_elems)
    print(f"  phys_tag: {dict(phys_tags)}")
    print(f"  geom_tag: {dict(geom_tags)}")

    print("\nConverting nodes to numpy array...")
    t0 = time.time()
    max_tag = max(nodes.keys())
    node_coords = np.zeros((max_tag + 1, 3), dtype=np.float64)
    for tag, (x, y, z) in nodes.items():
        node_coords[tag] = (x, y, z)
    print(f"  Done in {time.time()-t0:.2f}s")

    # Group hex elements by geom_tag (the actual volume assignment)
    print("\n=== Hex Elements by Geometry Tag (Volume Region) ===")
    hex_by_geom = defaultdict(list)
    for elem_id, phys_tag, geom_tag, node_list in hex_elems:
        hex_by_geom[geom_tag].append((elem_id, node_list))

    region_stats = {}
    for geom_tag in sorted(hex_by_geom.keys()):
        elems = hex_by_geom[geom_tag]
        n_elems = len(elems)

        centroids = np.zeros((n_elems, 3), dtype=np.float64)
        for i, (_, node_list) in enumerate(elems):
            coords = node_coords[node_list]
            centroids[i] = coords.mean(axis=0)

        bbox_min = centroids.min(axis=0)
        bbox_max = centroids.max(axis=0)
        bbox_size = bbox_max - bbox_min

        centroid = centroids.mean(axis=0)
        r_centroid = np.sqrt(centroid[0]**2 + centroid[1]**2)

        elem_r = np.sqrt(centroids[:, 0]**2 + centroids[:, 1]**2)
        elem_z = centroids[:, 2]
        elem_theta = np.degrees(np.arctan2(centroids[:, 1], centroids[:, 0]))

        region_stats[geom_tag] = {
            'geom_tag': geom_tag,
            'n_elements': n_elems,
            'bbox_min': bbox_min.tolist(),
            'bbox_max': bbox_max.tolist(),
            'bbox_size': bbox_size.tolist(),
            'centroid': centroid.tolist(),
            'r_centroid': float(r_centroid),
            'r_range': [float(elem_r.min()), float(elem_r.max())],
            'z_range': [float(elem_z.min()), float(elem_z.max())],
            'theta_range': [float(elem_theta.min()), float(elem_theta.max())],
        }

        print(f"\n  Geom tag={geom_tag} (R_{geom_tag-2}):")
        print(f"    Elements: {n_elems}")
        print(f"    BBox min: ({bbox_min[0]:.4f}, {bbox_min[1]:.4f}, {bbox_min[2]:.4f})")
        print(f"    BBox max: ({bbox_max[0]:.4f}, {bbox_max[1]:.4f}, {bbox_max[2]:.4f})")
        print(f"    BBox size: ({bbox_size[0]:.4f}, {bbox_size[1]:.4f}, {bbox_size[2]:.4f})")
        print(f"    Centroid: ({centroid[0]:.4f}, {centroid[1]:.4f}, {centroid[2]:.4f})")
        print(f"    r_centroid: {r_centroid:.4f}")
        print(f"    r_range: [{elem_r.min():.4f}, {elem_r.max():.4f}]")
        print(f"    z_range: [{elem_z.min():.4f}, {elem_z.max():.4f}]")
        print(f"    theta_range: [{elem_theta.min():.2f}, {elem_theta.max():.2f}]")

    with open(os.path.join(EVIDENCE_DIR, 'task-hex-core-region-stats.json'), 'w') as f:
        json.dump(region_stats, f, indent=2)

    # Analyze z-coordinates
    print("\n=== Z-Coordinate Analysis ===")
    all_hex_z_mins = []
    all_hex_z_maxs = []
    all_hex_z_centers = []
    for elem_id, phys_tag, geom_tag, node_list in hex_elems:
        coords = node_coords[node_list]
        all_hex_z_mins.append(coords[:, 2].min())
        all_hex_z_maxs.append(coords[:, 2].max())
        all_hex_z_centers.append(coords[:, 2].mean())

    all_hex_z_mins = np.array(all_hex_z_mins)
    all_hex_z_maxs = np.array(all_hex_z_maxs)
    all_hex_z_centers = np.array(all_hex_z_centers)
    all_hex_z_thickness = all_hex_z_maxs - all_hex_z_mins

    print(f"  Z mins: [{all_hex_z_mins.min():.6f}, {all_hex_z_mins.max():.6f}]")
    print(f"  Z maxs: [{all_hex_z_maxs.min():.6f}, {all_hex_z_maxs.max():.6f}]")
    print(f"  Z centers: [{all_hex_z_centers.min():.6f}, {all_hex_z_centers.max():.6f}]")
    print(f"  Z thickness: [{all_hex_z_thickness.min():.6f}, {all_hex_z_thickness.max():.6f}], mean={all_hex_z_thickness.mean():.6f}")

    # Check for discrete z-layers by looking at unique z-mins
    unique_z_mins = np.unique(np.round(all_hex_z_mins, 3))
    print(f"  Unique z-mins (3 decimals): {len(unique_z_mins)}")

    # Histogram of z-centers
    hist, bin_edges = np.histogram(all_hex_z_centers, bins=40)
    print("\n  Z-center histogram (top 10 bins):")
    sorted_bins = sorted(zip(hist, bin_edges[:-1], bin_edges[1:]), key=lambda x: -x[0])[:10]
    for cnt, lo, hi in sorted_bins:
        print(f"    [{lo:.4f}, {hi:.4f}): {cnt} hex")

    # Analyze connectivity: face-sharing between hex elements
    print("\n=== Hex Connectivity Analysis ===")
    t0 = time.time()

    hex_faces_idx = [
        (0, 1, 2, 3),  # bottom
        (4, 5, 6, 7),  # top
        (0, 1, 5, 4),  # front
        (1, 2, 6, 5),  # right
        (2, 3, 7, 6),  # back
        (3, 0, 4, 7),  # left
    ]

    face_to_elems = defaultdict(list)
    for elem_id, phys_tag, geom_tag, node_list in hex_elems:
        for face_idx in range(6):
            face = tuple(sorted(node_list[i] for i in hex_faces_idx[face_idx]))
            face_to_elems[face].append((elem_id, face_idx, geom_tag))

    print(f"  Face map built in {time.time()-t0:.2f}s")

    shared_faces = 0
    boundary_faces = 0
    inter_region_faces = 0
    non_manifold = 0
    for face, elems in face_to_elems.items():
        if len(elems) == 1:
            boundary_faces += 1
        elif len(elems) == 2:
            shared_faces += 1
            if elems[0][2] != elems[1][2]:
                inter_region_faces += 1
        else:
            non_manifold += 1

    total_faces = shared_faces + boundary_faces
    print(f"  Total unique faces: {total_faces}")
    print(f"  Shared faces (2 hex): {shared_faces}")
    print(f"  Boundary faces (1 hex): {boundary_faces}")
    print(f"  Inter-region shared faces: {inter_region_faces}")
    print(f"  Non-manifold faces (>2 hex): {non_manifold}")
    print(f"  Manifold ratio: {shared_faces / total_faces:.4f}")

    # Build node valence in hex connectivity graph
    print("\n=== Node Valence in Hex Connectivity ===")
    t0 = time.time()
    node_hex_count = defaultdict(int)
    for elem_id, phys_tag, geom_tag, node_list in hex_elems:
        for n in node_list:
            node_hex_count[n] += 1
    print(f"  Valence computed in {time.time()-t0:.2f}s")

    valence_dist = Counter(node_hex_count.values())
    print(f"  Total unique nodes in hex: {len(node_hex_count)}")
    print(f"  Valence distribution:")
    for v in sorted(valence_dist.keys()):
        print(f"    Valence {v}: {valence_dist[v]} nodes")

    # Find high-valence nodes (potential singularities)
    # In a regular hex grid, max valence is 8 (8 hexes meet at a corner)
    # Valence > 8 indicates T-junctions or multi-block meeting points
    high_valence = [(n, v) for n, v in node_hex_count.items() if v > 8]
    high_valence.sort(key=lambda x: -x[1])
    print(f"\n  High-valence nodes (valence > 8): {len(high_valence)}")
    for n, v in high_valence[:20]:
        x, y, z = node_coords[n]
        r = np.sqrt(x**2 + y**2)
        print(f"    Node {n}: valence={v}, pos=({x:.4f}, {y:.4f}, {z:.4f}), r={r:.4f}")

    # O-grid pattern detection
    print("\n=== O-Grid Pattern Detection ===")
    sorted_regions = sorted(region_stats.items(), key=lambda x: x[1]['r_centroid'])
    print("  Regions sorted by r_centroid (inner to outer):")
    for tag, stats in sorted_regions:
        print(f"    R_{tag-2} (geom_tag={tag}): r_centroid={stats['r_centroid']:.4f}, "
              f"r_range=[{stats['r_range'][0]:.4f}, {stats['r_range'][1]:.4f}], "
              f"n={stats['n_elements']}")

    # Check if regions form an O-grid pattern
    # O-grid: central block surrounded by peripheral blocks
    # Check radial ordering and coverage
    if len(sorted_regions) >= 3:
        print("\n  Checking O-grid pattern...")
        for i, (tag, stats) in enumerate(sorted_regions):
            print(f"    R_{tag-2}: r=[{stats['r_range'][0]:.3f}, {stats['r_range'][1]:.3f}], "
                  f"theta=[{stats['theta_range'][0]:.1f}, {stats['theta_range'][1]:.1f}], "
                  f"z=[{stats['z_range'][0]:.3f}, {stats['z_range'][1]:.3f}]")

    # Extract a single z-layer
    print("\n=== Z-Layer Extraction ===")
    # Pick a z-value that has many hex elements
    # Use the z-center with the most elements
    hist_max_idx = np.argmax(hist)
    target_z = (bin_edges[hist_max_idx] + bin_edges[hist_max_idx + 1]) / 2
    print(f"  Target z (most populated bin): {target_z:.4f}")

    # Find hex elements whose z-center is close to target_z
    z_tol = 0.02
    layer_hex_elems = []
    for elem_id, phys_tag, geom_tag, node_list in hex_elems:
        coords = node_coords[node_list]
        z_center = coords[:, 2].mean()
        if abs(z_center - target_z) < z_tol:
            layer_hex_elems.append((elem_id, geom_tag, node_list))

    print(f"  Hex elements near z={target_z:.4f} (tol={z_tol}): {len(layer_hex_elems)}")

    if len(layer_hex_elems) > 0:
        # Extract bottom faces
        layer_faces = []
        for elem_id, geom_tag, node_list in layer_hex_elems:
            face = tuple(sorted(node_list[i] for i in hex_faces_idx[0]))
            face_coords = node_coords[list(face)]
            if len(face_coords) == 4:
                layer_faces.append({
                    'elem_id': elem_id,
                    'geom_tag': geom_tag,
                    'face_nodes': face,
                    'coords': face_coords,
                })

        print(f"  Bottom faces extracted: {len(layer_faces)}")

        layer_nodes = set()
        for f in layer_faces:
            layer_nodes.update(f['face_nodes'])
        print(f"  Unique nodes in layer: {len(layer_nodes)}")

        layer_face_to_elems = defaultdict(list)
        for f in layer_faces:
            layer_face_to_elems[f['face_nodes']].append(f['elem_id'])

        layer_shared = sum(1 for elems in layer_face_to_elems.values() if len(elems) == 2)
        layer_boundary = sum(1 for elems in layer_face_to_elems.values() if len(elems) == 1)
        print(f"  Layer shared faces: {layer_shared}")
        print(f"  Layer boundary faces: {layer_boundary}")

        layer_node_face_count = defaultdict(int)
        for f in layer_faces:
            for n in f['face_nodes']:
                layer_node_face_count[n] += 1

        layer_valence_dist = Counter(layer_node_face_count.values())
        print(f"\n  Layer node valence distribution:")
        for v in sorted(layer_valence_dist.keys()):
            print(f"    Valence {v}: {layer_valence_dist[v]} nodes")

        layer_singular = [(n, v) for n, v in layer_node_face_count.items() if v != 4]
        print(f"\n  Layer singular nodes (valence != 4): {len(layer_singular)}")
        for n, v in sorted(layer_singular, key=lambda x: -x[1])[:20]:
            x, y, z = node_coords[n]
            r = np.sqrt(x**2 + y**2)
            print(f"    Node {n}: valence={v}, pos=({x:.4f}, {y:.4f}, {z:.4f}), r={r:.4f}")

        # Save layer data
        layer_data = {
            'z_center': float(target_z),
            'n_hex': len(layer_hex_elems),
            'n_faces': len(layer_faces),
            'n_nodes': len(layer_nodes),
            'n_shared_faces': layer_shared,
            'n_boundary_faces': layer_boundary,
            'faces': [
                {
                    'elem_id': f['elem_id'],
                    'geom_tag': f['geom_tag'],
                    'nodes': list(f['face_nodes']),
                }
                for f in layer_faces
            ],
            'nodes': {str(n): node_coords[n].tolist() for n in layer_nodes},
        }
        with open(os.path.join(EVIDENCE_DIR, 'task-hex-core-zlayer.json'), 'w') as f:
            json.dump(layer_data, f, indent=2)

        # Generate visualization
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt

            fig, axes = plt.subplots(1, 2, figsize=(16, 8))

            ax = axes[0]
            region_colors = {2: 'red', 3: 'blue', 4: 'green', 5: 'orange', 6: 'purple'}
            for f in layer_faces:
                coords = f['coords']
                poly = plt.Polygon(coords[:, :2], alpha=0.3,
                                   facecolor=region_colors.get(f['geom_tag'], 'gray'),
                                   edgecolor='black', linewidth=0.3)
                ax.add_patch(poly)

            ax.set_aspect('equal')
            ax.set_title(f'Hex Bottom Faces at z={target_z:.3f}\n'
                         f'(colored by geometry tag)')
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.grid(True, alpha=0.3)

            ax = axes[1]
            for f in layer_faces:
                coords = f['coords']
                poly = plt.Polygon(coords[:, :2], alpha=0.2,
                                   facecolor='lightblue',
                                   edgecolor='gray', linewidth=0.2)
                ax.add_patch(poly)

            for n, v in layer_singular:
                x, y, z = node_coords[n]
                color = 'red' if v > 4 else 'orange'
                size = 50 + abs(v - 4) * 20
                ax.scatter(x, y, c=color, s=size, edgecolors='black', linewidths=0.5, zorder=5)
                ax.annotate(f'v={v}', (x, y), fontsize=6, ha='center', va='center')

            ax.set_aspect('equal')
            ax.set_title(f'Singular Nodes at z={target_z:.3f}\n'
                         f'({len(layer_singular)} nodes with valence != 4)')
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.grid(True, alpha=0.3)

            plt.tight_layout()
            plt.savefig(os.path.join(EVIDENCE_DIR, 'task-hex-core-zlayer.png'), dpi=150)
            plt.close()
            print(f"\n  Visualization saved to {EVIDENCE_DIR}/task-hex-core-zlayer.png")
        except ImportError:
            print("\n  matplotlib not available, skipping visualization")

    # Summary
    print("\n" + "="*60)
    print("=== SUMMARY ===")
    print("="*60)
    print(f"Total hex elements: {len(hex_elems)}")
    print(f"Volume regions (geom_tag): {len(hex_by_geom)}")
    for tag, stats in sorted(region_stats.items()):
        print(f"  R_{tag-2} (geom_tag={tag}): {stats['n_elements']} elements, "
              f"r=[{stats['r_range'][0]:.3f}, {stats['r_range'][1]:.3f}], "
              f"z=[{stats['z_range'][0]:.3f}, {stats['z_range'][1]:.3f}], "
              f"theta=[{stats['theta_range'][0]:.1f}, {stats['theta_range'][1]:.1f}]")
    print(f"\nZ range: [{all_hex_z_mins.min():.3f}, {all_hex_z_maxs.max():.3f}]")
    print(f"Z thickness range: [{all_hex_z_thickness.min():.4f}, {all_hex_z_thickness.max():.4f}]")
    print(f"Manifold ratio: {shared_faces / total_faces:.4f}")
    print(f"High-valence nodes (>8): {len(high_valence)}")
    if len(layer_hex_elems) > 0:
        print(f"Layer singular nodes (valence != 4): {len(layer_singular)}")


if __name__ == '__main__':
    main()
