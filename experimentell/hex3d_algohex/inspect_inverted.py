#!/usr/bin/env python3
"""Inspect inverted cells in Canada hex meshes.

Parses hex meshes (OVM and VTK formats), computes per-cell scaled Jacobian,
identifies inverted cells (min sJ <= 0), and creates visualization files
with inverted cells, their neighbors, and markers.

For each inverted cell reports:
- mesh name, cell id
- centre xyz, bbox, volume
- min scaled Jacobian
- boundary vs interior (count of boundary faces)
- 1-ring neighbours

Writes:
- inverted_zoom.vtk: inverted cells + 1-2 rings of neighbours + markers
- inverted_zoom.json: metadata and explanation
"""

import sys
import json
import numpy as np
from pathlib import Path
from collections import defaultdict

# Add the module directory to path
REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "experimentell" / "hex3d_algohex"))

import ovm_io


def read_hex_mesh_vtk(path):
    """Read a hex mesh from VTK unstructured grid format.
    
    Returns (points, hexes) where hexes is (n_cells, 8) array.
    """
    with open(path) as f:
        lines = f.readlines()
    
    i = 0
    # Find POINTS section
    while i < len(lines) and not lines[i].startswith("POINTS"):
        i += 1
    if i >= len(lines):
        raise ValueError(f"No POINTS section in {path}")
    
    n_pts = int(lines[i].split()[1])
    i += 1
    points = np.zeros((n_pts, 3))
    for k in range(n_pts):
        parts = lines[i].split()
        points[k] = [float(x) for x in parts[:3]]
        i += 1
    
    # Find CELLS section
    while i < len(lines) and not lines[i].startswith("CELLS"):
        i += 1
    if i >= len(lines):
        raise ValueError(f"No CELLS section in {path}")
    
    n_cells = int(lines[i].split()[1])
    i += 1
    hexes = []
    for k in range(n_cells):
        parts = [int(x) for x in lines[i].split()]
        if parts[0] == 8:  # VTK_HEXAHEDRON
            hexes.append(parts[1:9])
        i += 1
    
    return points, np.array(hexes, dtype=np.int64)


def read_hex_mesh(path):
    """Read hex mesh from OVM or VTK format."""
    path = Path(path)
    if path.suffix == ".ovm":
        points, edges, faces, polys = ovm_io.read_ovm(path)
        hexes, skipped = ovm_io.ovm_to_cells(points, edges, faces, polys)
        return points, hexes
    elif path.suffix == ".vtk":
        return read_hex_mesh_vtk(path)
    else:
        raise ValueError(f"Unknown format: {path.suffix}")


def compute_scaled_jacobians(points, hexes):
    """Compute min scaled Jacobian for each cell.
    
    Uses the same convention as ovm_io.scaled_jacobian and HexBlockValidator.
    Returns array of shape (n_cells,).
    """
    idx = [(0, 1, 3, 4), (1, 2, 0, 5), (2, 3, 1, 6), (3, 0, 2, 7),
           (4, 7, 5, 0), (5, 4, 6, 1), (6, 5, 7, 2), (7, 6, 4, 3)]
    
    Q = points[hexes]  # (n_cells, 8, 3)
    sj_all = np.empty((len(hexes), 8))
    
    for k, (a, b, c, d) in enumerate(idx):
        u = Q[:, b] - Q[:, a]
        v = Q[:, c] - Q[:, a]
        w = Q[:, d] - Q[:, a]
        n = np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1) * np.linalg.norm(w, axis=1)
        denom = np.where(n > 0, n, 1.0)
        sj_all[:, k] = np.einsum('ij,ij->i', u, np.cross(v, w)) / denom
    
    return sj_all.min(axis=1)


def build_cell_topology(hexes):
    """Build face-to-cell and edge-to-cell maps.
    
    Returns (f2h, e2h) where:
    - f2h: dict mapping frozenset(face_verts) -> list of cell indices
    - e2h: dict mapping (min_v, max_v) -> list of cell indices
    """
    f2h = defaultdict(list)
    e2h = defaultdict(list)
    
    # Hex faces in VTK order
    HF = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
          (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
    
    for ci, cell in enumerate(hexes):
        for face in HF:
            verts = tuple(int(cell[i]) for i in face)
            f2h[frozenset(verts)].append(ci)
        
        # Edges
        edges = [(0, 1), (1, 2), (2, 3), (3, 0),
                 (4, 5), (5, 6), (6, 7), (7, 4),
                 (0, 4), (1, 5), (2, 6), (3, 7)]
        for a, b in edges:
            u, v = int(cell[a]), int(cell[b])
            e2h[(min(u, v), max(u, v))].append(ci)
    
    return f2h, e2h


def compute_cell_neighbours(f2h, n_cells):
    """Compute 1-ring cell neighbours (cells sharing a face)."""
    neighbours = defaultdict(set)
    for fk, hs in f2h.items():
        if len(hs) == 2:
            a, b = hs
            neighbours[a].add(b)
            neighbours[b].add(a)
    return neighbours


def compute_cell_boundary_faces(hexes, f2h):
    """Count boundary faces per cell (faces with only one owner)."""
    bnd_faces = {fk for fk, hs in f2h.items() if len(hs) == 1}
    
    HF = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
          (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
    
    bnd_count = np.zeros(len(hexes), dtype=int)
    for ci, cell in enumerate(hexes):
        for face in HF:
            verts = tuple(int(cell[i]) for i in face)
            if frozenset(verts) in bnd_faces:
                bnd_count[ci] += 1
    
    return bnd_count


def compute_cell_volume(points, cell):
    """Compute signed volume of a hexahedron."""
    return ovm_io._hex_volume(points[cell])


def compute_cell_center(points, cell):
    """Compute cell centroid."""
    return points[cell].mean(axis=0)


def compute_cell_bbox(points, cell):
    """Compute cell bounding box."""
    pts = points[cell]
    return pts.min(axis=0), pts.max(axis=0)


def get_rings_around(cells, neighbours, n_rings=1):
    """Get cells within n_rings of the given cells."""
    result = set(cells)
    frontier = set(cells)
    for _ in range(n_rings):
        new_frontier = set()
        for c in frontier:
            new_frontier.update(neighbours[c])
        new_frontier -= result
        result.update(new_frontier)
        frontier = new_frontier
    return sorted(result)


def write_zoom_vtk(path, points, hexes, inverted_ids, ring_ids, cell_data=None):
    """Write a VTK file with inverted cells, their rings, and markers.
    
    Markers are small octahedra at the center of each inverted cell,
    scaled to ~2x the cell size for visibility.
    """
    # Collect all cell indices
    all_cells = sorted(set(inverted_ids) | set(ring_ids))
    
    # Compact points to only those used
    used_pts = np.unique(np.concatenate([hexes[c] for c in all_cells]))
    pt_remap = {int(v): i for i, v in enumerate(used_pts)}
    P = points[used_pts]
    H = np.array([[pt_remap[int(v)] for v in hexes[c]] for c in all_cells])
    
    # Map cell indices
    cell_remap = {old: new for new, old in enumerate(all_cells)}
    inv_local = [cell_remap[c] for c in inverted_ids]
    
    # Create markers (octahedra) at inverted cell centers
    markers = []
    marker_pts = []
    pt_offset = len(P)
    
    for cid in inverted_ids:
        center = compute_cell_center(points, hexes[cid])
        # Estimate cell size from bbox
        bmin, bmax = compute_cell_bbox(points, hexes[cid])
        size = np.linalg.norm(bmax - bmin)
        radius = size * 1.0  # marker radius = cell size
        
        # Octahedron: 6 vertices along axes
        oct_pts = [
            center + [radius, 0, 0],
            center + [-radius, 0, 0],
            center + [0, radius, 0],
            center + [0, -radius, 0],
            center + [0, 0, radius],
            center + [0, 0, -radius],
        ]
        
        # Octahedron faces (8 triangles)
        # Vertices: 0=+x, 1=-x, 2=+y, 3=-y, 4=+z, 5=-z
        oct_faces = [
            [0, 2, 4], [2, 1, 4], [1, 3, 4], [3, 0, 4],
            [0, 4, 2], [2, 4, 1], [1, 4, 3], [3, 4, 0],
            [0, 2, 5], [2, 1, 5], [1, 3, 5], [3, 0, 5],
            [0, 5, 2], [2, 5, 1], [1, 5, 3], [3, 5, 0],
        ]
        
        for pt in oct_pts:
            marker_pts.append(pt)
        
        for face in oct_faces:
            markers.append([pt_offset + f for f in face])
        
        pt_offset += 6
    
    # Write VTK
    with open(path, 'w') as f:
        f.write("# vtk DataFile Version 2.0\n")
        f.write("Inverted cells with 1-2 ring neighbors and markers\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        
        # Points (hex mesh + marker points)
        all_pts = np.vstack([P, marker_pts]) if marker_pts else P
        f.write(f"POINTS {len(all_pts)} float\n")
        for pt in all_pts:
            f.write(f"{pt[0]} {pt[1]} {pt[2]}\n")
        
        # Cells
        n_hex = len(H)
        n_tri = len(markers)
        n_total = n_hex + n_tri
        f.write(f"CELLS {n_total} {9 * n_hex + 4 * n_tri}\n")
        
        for h in H:
            f.write(f"8 {' '.join(str(int(v)) for v in h)}\n")
        
        for m in markers:
            f.write(f"3 {m[0]} {m[1]} {m[2]}\n")
        
        # Cell types
        f.write(f"CELL_TYPES {n_total}\n")
        for _ in range(n_hex):
            f.write("12\n")  # VTK_HEXAHEDRON
        for _ in range(n_tri):
            f.write("5\n")  # VTK_TRIANGLE
        
        # Cell data
        if cell_data or inv_local:
            f.write(f"CELL_DATA {n_total}\n")
            
            # Inverted flag
            f.write("SCALARS is_inverted int 1\n")
            f.write("LOOKUP_TABLE default\n")
            inv_set = set(inv_local)
            for i in range(n_hex):
                f.write(f"{1 if i in inv_set else 0}\n")
            for _ in range(n_tri):
                f.write("1\n")  # markers are all "inverted"
            
            # Original cell ID
            f.write("SCALARS original_cell_id int 1\n")
            f.write("LOOKUP_TABLE default\n")
            for c in all_cells:
                f.write(f"{c}\n")
            for _ in range(n_tri):
                f.write("-1\n")  # markers have no original ID


def write_scaled_vtk(path, points, hexes, inverted_ids, scale=1.5):
    """Write inverted cells with coordinates scaled about their centers."""
    # For each inverted cell, scale its coordinates
    new_points = points.copy()
    new_hexes = []
    
    for cid in inverted_ids:
        center = compute_cell_center(points, hexes[cid])
        cell = hexes[cid]
        
        # Scale each vertex about the center
        for v in cell:
            new_points[v] = center + scale * (points[v] - center)
        
        new_hexes.append(cell)
    
    # Write as separate VTK
    H = np.array(new_hexes)
    used_pts = np.unique(H)
    pt_remap = {int(v): i for i, v in enumerate(used_pts)}
    P = new_points[used_pts]
    H_remapped = np.array([[pt_remap[int(v)] for v in h] for h in H])
    
    with open(path, 'w') as f:
        f.write("# vtk DataFile Version 2.0\n")
        f.write(f"Scaled inverted cells (scale={scale}x about center)\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {len(P)} float\n")
        for pt in P:
            f.write(f"{pt[0]} {pt[1]} {pt[2]}\n")
        f.write(f"CELLS {len(H_remapped)} {9 * len(H_remapped)}\n")
        for h in H_remapped:
            f.write(f"8 {' '.join(str(int(v)) for v in h)}\n")
        f.write(f"CELL_TYPES {len(H_remapped)}\n")
        for _ in H_remapped:
            f.write("12\n")


def check_existing_inverted_vtk(path):
    """Check what the existing inverted_cells.vtk contains."""
    print(f"\n=== Checking existing {path} ===")
    if not Path(path).exists():
        print(f"File not found: {path}")
        return
    
    points, hexes = read_hex_mesh_vtk(path)
    print(f"Contains {len(hexes)} cells")
    
    # Check for duplicates
    cell_set = set(tuple(h) for h in hexes)
    print(f"Unique cells: {len(cell_set)} (duplicates: {len(hexes) - len(cell_set)})")
    
    # Check point duplication
    pt_set = set(tuple(p) for p in points)
    print(f"Unique points: {len(pt_set)} (duplicates: {len(points) - len(pt_set)})")
    
    return points, hexes


def main():
    import argparse
    
    ap = argparse.ArgumentParser(description="Inspect inverted cells in Canada hex meshes")
    ap.add_argument("--check-existing", action="store_true",
                    help="Check existing inverted_cells.vtk")
    args = ap.parse_args()
    
    # Define input files
    ovm_path = Path("/work/canada/sweep/n100000/hex_n100000.ovm")
    repaired_vtk = Path("/work/canada/sweep/inspect/repaired_n100000.vtk")
    blocks_vtk = Path("/work/canada/sweep/struct_beam/blocks.vtk")
    existing_inv = Path("/work/canada/sweep/inspect/inverted_cells.vtk")
    out_dir = Path("/work/canada/sweep/inspect")
    
    # Check existing file if requested
    if args.check_existing:
        check_existing_inverted_vtk(existing_inv)
        return
    
    # Collect inverted cells from all meshes
    all_inverted = []
    
    # 1. Original OVM mesh
    print("\n=== Processing original OVM mesh ===")
    if ovm_path.exists():
        P, H = read_hex_mesh(ovm_path)
        sj = compute_scaled_jacobians(P, H)
        inv_ids = np.where(sj <= 0)[0]
        print(f"Original mesh: {len(H)} cells, {len(inv_ids)} inverted, min sJ {sj.min():.4f}")
        
        for cid in inv_ids:
            center = compute_cell_center(P, H[cid])
            bmin, bmax = compute_cell_bbox(P, H[cid])
            vol = compute_cell_volume(P, H[cid])
            all_inverted.append({
                "mesh": "hex_n100000.ovm",
                "cell_id": int(cid),
                "center": center.tolist(),
                "bbox_min": bmin.tolist(),
                "bbox_max": bmax.tolist(),
                "volume": float(vol),
                "min_sj": float(sj[cid]),
            })
    
    # 2. Repaired mesh
    print("\n=== Processing repaired mesh ===")
    if repaired_vtk.exists():
        P, H = read_hex_mesh(repaired_vtk)
        sj = compute_scaled_jacobians(P, H)
        inv_ids = np.where(sj <= 0)[0]
        print(f"Repaired mesh: {len(H)} cells, {len(inv_ids)} inverted, min sJ {sj.min():.4f}")
        
        f2h, e2h = build_cell_topology(H)
        neighbours = compute_cell_neighbours(f2h, len(H))
        bnd_count = compute_cell_boundary_faces(H, f2h)
        
        for cid in inv_ids:
            center = compute_cell_center(P, H[cid])
            bmin, bmax = compute_cell_bbox(P, H[cid])
            vol = compute_cell_volume(P, H[cid])
            ring1 = get_rings_around([cid], neighbours, 1)
            ring2 = get_rings_around([cid], neighbours, 2)
            
            all_inverted.append({
                "mesh": "repaired_n100000.vtk",
                "cell_id": int(cid),
                "center": center.tolist(),
                "bbox_min": bmin.tolist(),
                "bbox_max": bmax.tolist(),
                "volume": float(vol),
                "min_sj": float(sj[cid]),
                "n_boundary_faces": int(bnd_count[cid]),
                "is_interior": bnd_count[cid] == 0,
                "ring1_neighbours": sorted(ring1),
                "ring2_neighbours": sorted(ring2),
                "points": P,
                "hexes": H,
            })
    
    # 3. Blocks mesh
    print("\n=== Processing blocks mesh ===")
    if blocks_vtk.exists():
        P, H = read_hex_mesh(blocks_vtk)
        sj = compute_scaled_jacobians(P, H)
        inv_ids = np.where(sj <= 0)[0]
        print(f"Blocks mesh: {len(H)} cells, {len(inv_ids)} inverted, min sJ {sj.min():.4f}")
        
        f2h, e2h = build_cell_topology(H)
        neighbours = compute_cell_neighbours(f2h, len(H))
        bnd_count = compute_cell_boundary_faces(H, f2h)
        
        for cid in inv_ids:
            center = compute_cell_center(P, H[cid])
            bmin, bmax = compute_cell_bbox(P, H[cid])
            vol = compute_cell_volume(P, H[cid])
            ring1 = get_rings_around([cid], neighbours, 1)
            ring2 = get_rings_around([cid], neighbours, 2)
            
            all_inverted.append({
                "mesh": "blocks.vtk",
                "cell_id": int(cid),
                "center": center.tolist(),
                "bbox_min": bmin.tolist(),
                "bbox_max": bmax.tolist(),
                "volume": float(vol),
                "min_sj": float(sj[cid]),
                "n_boundary_faces": int(bnd_count[cid]),
                "is_interior": bnd_count[cid] == 0,
                "ring1_neighbours": sorted(ring1),
                "ring2_neighbours": sorted(ring2),
                "points": P,
                "hexes": H,
            })
    
    # Print summary
    print(f"\n=== Summary ===")
    print(f"Total inverted cells found: {len(all_inverted)}")
    for inv in all_inverted:
        print(f"\nMesh: {inv['mesh']}, Cell ID: {inv['cell_id']}")
        print(f"  Center: {inv['center']}")
        print(f"  BBox: {inv['bbox_min']} -> {inv['bbox_max']}")
        print(f"  Volume: {inv['volume']:.6f}")
        print(f"  Min sJ: {inv['min_sj']:.6f}")
        if 'is_interior' in inv:
            print(f"  Interior: {inv['is_interior']} ({inv.get('n_boundary_faces', '?')} boundary faces)")
        if 'ring1_neighbours' in inv:
            print(f"  Ring-1 neighbours: {len(inv['ring1_neighbours'])} cells")
            print(f"  Ring-2 neighbours: {len(inv['ring2_neighbours'])} cells")
    
    # Write visualization files for the first inverted cell with full data
    inv_with_data = [x for x in all_inverted if 'points' in x]
    if inv_with_data:
        inv = inv_with_data[0]
        P = inv['points']
        H = inv['hexes']
        cid = inv['cell_id']
        
        f2h, e2h = build_cell_topology(H)
        neighbours = compute_cell_neighbours(f2h, len(H))
        
        # Get 1-2 rings
        ring1 = inv['ring1_neighbours']
        ring2 = inv['ring2_neighbours']
        
        # Write inverted_zoom.vtk
        zoom_path = out_dir / "inverted_zoom.vtk"
        print(f"\n=== Writing {zoom_path} ===")
        write_zoom_vtk(zoom_path, P, H, [cid], ring2)
        print(f"Wrote {zoom_path}")
        
        # Write scaled version if cell is tiny
        vol = abs(inv['volume'])
        if vol < 1e-6:
            scaled_path = out_dir / "inverted_zoom_scaled.vtk"
            print(f"Cell volume is tiny ({vol:.2e}), writing scaled version")
            write_scaled_vtk(scaled_path, P, H, [cid], scale=1.5)
            print(f"Wrote {scaled_path}")
        
        # Write metadata JSON
        json_path = out_dir / "inverted_zoom.json"
        print(f"\n=== Writing {json_path} ===")
        
        # Check if cell is invisible in surface views
        is_interior = inv.get('is_interior', False)
        is_degenerate = vol < 1e-6
        
        boundary_desc = 'has no boundary faces' if is_interior else f"has {inv.get('n_boundary_faces', 0)} boundary faces"
        interior_desc = 'it is completely enclosed by other cells and does not touch the domain boundary, ' if is_interior else ''
        degenerate_desc = 'its near-zero volume makes it collapse to a point/line/plane, ' if is_degenerate else ''
        
        explanation = (
            f"This inverted cell is {'interior' if is_interior else 'on the boundary'} "
            f"({boundary_desc}). "
            f"It {'is' if is_degenerate else 'is not'} geometrically degenerate (volume = {vol:.2e}). "
            f"The cell is invisible in normal surface views because "
            f"{interior_desc}"
            f"{degenerate_desc}"
            f"surface-only exports only show boundary faces, not interior cells or their quality. "
            f"To visualize it, use the inverted_zoom.vtk file which shows the cell, its 2-ring neighbors, "
            f"and an octahedral marker at its center."
        )
        
        metadata = {
            "mesh": inv['mesh'],
            "cell_id": inv['cell_id'],
            "center": inv['center'],
            "bbox_min": inv['bbox_min'],
            "bbox_max": inv['bbox_max'],
            "volume": inv['volume'],
            "min_sj": inv['min_sj'],
            "is_interior": bool(is_interior),
            "n_boundary_faces": int(inv.get('n_boundary_faces', 0)),
            "ring1_neighbours": [int(x) for x in inv['ring1_neighbours']],
            "ring2_neighbours": [int(x) for x in inv['ring2_neighbours']],
            "explanation": explanation,
        }
        
        with open(json_path, 'w') as f:
            json.dump(metadata, f, indent=2)
        print(f"Wrote {json_path}")
    else:
        print("\nNo inverted cells with full mesh data found, skipping visualization files")


if __name__ == "__main__":
    main()
