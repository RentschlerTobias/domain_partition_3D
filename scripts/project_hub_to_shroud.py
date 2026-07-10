#!/usr/bin/env python3
"""
Project Hub singularities to Shroud along blade curves.

This script:
1. Loads existing Hub singularities from hub_master_ta.json
2. Extracts blade curves from .msh file (geometry tags 5,6)
3. For each Hub singularity, finds the nearest point on the blade curve at z=0
4. Follows the blade curve from hub to shroud (z=0 to z=2.5)
5. Projects the singularity position to z=2.5
6. Validates against actual Shroud geometry
"""

import json
import numpy as np
from pathlib import Path

def load_hub_singularities():
    """Load singularities from existing Hub partition."""
    sing_file = Path('output/tmesh_hub/master/hub_master_ta.json')
    if not sing_file.exists():
        print(f"Warning: {sing_file} not found")
        return []
    
    with open(sing_file, 'r') as f:
        data = json.load(f)
    
    singularities = []
    for sing in data.get('singularities', []):
        singularities.append({
            'id': sing['id'],
            'position': sing['position_st'],
            'index': sing.get('index', -1)
        })
    
    return singularities

def extract_blade_curves(msh_file):
    """Extract blade surface curves from .msh file."""
    nodes = {}
    blade_nodes = set()
    
    with open(msh_file, 'r') as f:
        line = f.readline()
        while line:
            line = line.strip()
            
            if line == '$Nodes':
                n_nodes = int(f.readline().strip())
                for _ in range(n_nodes):
                    parts = f.readline().strip().split()
                    node_id = int(parts[0])
                    coords = [float(parts[1]), float(parts[2]), float(parts[3])]
                    nodes[node_id] = coords
            
            elif line == '$Elements':
                n_elements = int(f.readline().strip())
                for _ in range(n_elements):
                    parts = f.readline().strip().split()
                    elem_type = int(parts[1])
                    n_tags = int(parts[2])
                    tags = [int(p) for p in parts[3:3+n_tags]]
                    node_ids = [int(p) for p in parts[3+n_tags:]]
                    
                    geom_tag = tags[1] if len(tags) > 1 else 0
                    
                    # Collect blade surface nodes (geom tags 5,6)
                    if geom_tag in [5, 6] and elem_type in [2, 3]:
                        blade_nodes.update(node_ids)
            
            line = f.readline()
    
    # Organize blade nodes by z-level
    blade_by_z = {}
    for nid in blade_nodes:
        x, y, z = nodes[nid]
        z_rounded = round(z, 3)
        if z_rounded not in blade_by_z:
            blade_by_z[z_rounded] = []
        blade_by_z[z_rounded].append([x, y, z])
    
    return nodes, blade_by_z

def project_singularity_to_shroud(singularity, blade_by_z):
    """Project a Hub singularity to Shroud along blade curves."""
    sx, sy = singularity['position']
    
    # Find nearest blade point at z=0
    hub_blade = blade_by_z.get(0.0, [])
    if not hub_blade:
        print(f"Warning: No blade nodes at z=0 for {singularity['id']}")
        return None
    
    hub_blade = np.array(hub_blade)
    distances = np.sqrt((hub_blade[:,0] - sx)**2 + (hub_blade[:,1] - sy)**2)
    nearest_idx = np.argmin(distances)
    nearest_hub = hub_blade[nearest_idx]
    
    print(f"  {singularity['id']} at ({sx:.4f}, {sy:.4f})")
    print(f"    Nearest blade point at z=0: ({nearest_hub[0]:.4f}, {nearest_hub[1]:.4f}, {nearest_hub[2]:.4f})")
    print(f"    Distance to blade: {distances[nearest_idx]:.4f}")
    
    # Follow blade curve to z=2.5
    # Strategy: Find corresponding point at z=2.5 based on parametric position
    shroud_blade = blade_by_z.get(2.5, [])
    if not shroud_blade:
        print(f"    Warning: No blade nodes at z=2.5")
        return None
    
    shroud_blade = np.array(shroud_blade)
    
    # Compute parametric position (normalized arc length along blade at z=0)
    hub_sorted = hub_blade[np.argsort(np.arctan2(hub_blade[:,1], hub_blade[:,0]))]
    
    # Find the angular position of the nearest hub point
    nearest_angle = np.arctan2(nearest_hub[1], nearest_hub[0])
    
    # Find corresponding point at z=2.5 with similar angular position
    shroud_angles = np.arctan2(shroud_blade[:,1], shroud_blade[:,0])
    angle_diff = np.abs(shroud_angles - nearest_angle)
    # Handle angle wrapping
    angle_diff = np.minimum(angle_diff, 2*np.pi - angle_diff)
    
    shroud_idx = np.argmin(angle_diff)
    nearest_shroud = shroud_blade[shroud_idx]
    
    print(f"    Corresponding point at z=2.5: ({nearest_shroud[0]:.4f}, {nearest_shroud[1]:.4f}, {nearest_shroud[2]:.4f})")
    print(f"    Angular difference: {angle_diff[shroud_idx]:.4f} rad ({np.degrees(angle_diff[shroud_idx]):.2f}°)")
    
    return {
        'singularity': singularity,
        'hub_blade_point': nearest_hub.tolist(),
        'shroud_blade_point': nearest_shroud.tolist(),
        'projected_position': [nearest_shroud[0], nearest_shroud[1]],
        'distance_to_blade': float(distances[nearest_idx]),
        'angular_diff': float(angle_diff[shroud_idx])
    }

def main():
    print("Loading Hub singularities...")
    singularities = load_hub_singularities()
    print(f"Found {len(singularities)} singularities")
    
    print("\nExtracting blade curves from mesh...")
    msh_file = 'T1_9/T1_9_ru_gridGmsh.msh'
    nodes, blade_by_z = extract_blade_curves(msh_file)
    
    print(f"Found blade nodes at {len(blade_by_z)} z-levels")
    for z in sorted(blade_by_z.keys())[:5]:
        print(f"  z={z}: {len(blade_by_z[z])} nodes")
    
    print(f"\nProjecting singularities to Shroud...")
    print("=" * 60)
    
    projections = []
    for sing in singularities:
        proj = project_singularity_to_shroud(sing, blade_by_z)
        if proj:
            projections.append(proj)
        print()
    
    # Summary
    print("=" * 60)
    print("\nProjection Summary:")
    print(f"Total singularities: {len(singularities)}")
    print(f"Successfully projected: {len(projections)}")
    
    # Compute statistics
    if projections:
        distances = [p['distance_to_blade'] for p in projections]
        angles = [p['angular_diff'] for p in projections]
        
        print(f"\nDistance to blade at Hub:")
        print(f"  Min: {min(distances):.4f}")
        print(f"  Max: {max(distances):.4f}")
        print(f"  Mean: {np.mean(distances):.4f}")
        
        print(f"\nAngular difference Hub→Shroud:")
        print(f"  Min: {np.degrees(min(angles)):.2f}°")
        print(f"  Max: {np.degrees(max(angles)):.2f}°")
        print(f"  Mean: {np.degrees(np.mean(angles)):.2f}°")
    
    # Save results
    output = {
        'metadata': {
            'n_singularities': len(singularities),
            'n_projected': len(projections),
            'mesh_file': msh_file
        },
        'projections': projections
    }
    
    with open('hub_to_shroud_projection.json', 'w') as f:
        json.dump(output, f, indent=2, default=lambda x: float(x) if isinstance(x, np.floating) else x)
    
    print(f"\nSaved to hub_to_shroud_projection.json")
    
    # Check if projection is feasible
    print("\n" + "=" * 60)
    print("FEASIBILITY ASSESSMENT:")
    print("=" * 60)
    
    if len(projections) == len(singularities):
        print("✓ All singularities can be projected to Shroud")
    else:
        print("✗ Some singularities could not be projected")
    
    if projections and max(distances) < 0.1:
        print("✓ All singularities are close to blade surface (< 0.1)")
        print("  → Projection along blade curves is FEASIBLE")
    elif projections:
        print(f"⚠ Some singularities are far from blade (max: {max(distances):.4f})")
        print("  → May need interpolation or different projection method")
    
    if projections and np.degrees(max(angles)) < 5:
        print("✓ Angular difference Hub→Shroud is small (< 5°)")
        print("  → Blade twist is manageable")
    elif projections:
        print(f"⚠ Blade twist is significant (max: {np.degrees(max(angles)):.2f}°)")
        print("  → Need to account for twist in projection")

if __name__ == '__main__':
    main()
