#!/usr/bin/env python3
"""
Generate 28 block instances (7 blocks × 4 rotations) and create interactive 3D plot.

This script:
1. Computes all 7 base blocks from dtOO XML parameters
2. Applies 4-fold rotational symmetry (0°, 90°, 180°, 270°)
3. Generates 28 block instances with unique labels
4. Creates an interactive Plotly 3D visualization with:
   - Block wireframes
   - Corner labels
   - Center labels
   - Color coding by block type
"""

import json
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path

def compute_base_blocks(params):
    """Compute the 7 base blocks in parametric space."""
    p00 = params['cV_ru_divideInternalMeshBlock_0_0']
    p01 = params['cV_ru_divideInternalMeshBlock_0_1']
    p10 = params['cV_ru_divideInternalMeshBlock_1_0']
    p11 = params['cV_ru_divideInternalMeshBlock_1_1']
    base = 0.495  # for 'ru' context
    
    blocks = {}
    
    # Trailing Edge blocks
    blocks['te_0'] = {
        'corners': [
            [base + p00, 0, 0], [1.0, 0, 0],
            [base + p01, 0, 1], [1.0, 0, 1],
            [base + p00, 1, 0], [1.0, 1, 0],
            [base + p01, 1, 1], [1.0, 1, 1]
        ],
        'location': 'trailing_edge',
        'side': 'pressure'
    }
    
    blocks['te_1'] = {
        'corners': [
            [0, 0, 0], [p10 - base, 0, 0],
            [0, 0, 1], [p11 - base, 0, 1],
            [0, 1, 0], [p10 - base, 1, 0],
            [0, 1, 1], [p11 - base, 1, 1]
        ],
        'location': 'trailing_edge',
        'side': 'suction'
    }
    
    # Main blocks
    blocks['0_0'] = {
        'corners': [
            [p10 - 0.43, 0, 0], [p00 - 0.025, 0, 0],
            [p11 - 0.43, 0, 1], [p01 - 0.025, 0, 1],
            [p10 - 0.43, 1, 0], [p00 - 0.025, 1, 0],
            [p11 - 0.43, 1, 1], [p01 - 0.025, 1, 1]
        ],
        'location': 'main',
        'side': 'suction'
    }
    
    blocks['0_1'] = {
        'corners': [
            [p00 - 0.025, 0, 0], [0.025 + p10, 0, 0],
            [p01 - 0.025, 0, 1], [0.025 + p11, 0, 1],
            [p00 - 0.025, 1, 0], [0.025 + p10, 1, 0],
            [p01 - 0.025, 1, 1], [0.025 + p11, 1, 1]
        ],
        'location': 'main',
        'side': 'center'
    }
    
    blocks['0_2'] = {
        'corners': [
            [0.085 + p10, 0, 0], [0.43 + p00, 0, 0],
            [0.085 + p11, 0, 1], [0.43 + p01, 0, 1],
            [0.085 + p10, 1, 0], [0.43 + p00, 1, 0],
            [0.085 + p11, 1, 1], [0.43 + p01, 1, 1]
        ],
        'location': 'main',
        'side': 'pressure'
    }
    
    # Leading Edge blocks
    blocks['le_0'] = {
        'corners': [
            [p00 - 0.085, 0, 0], [0.5, 0, 0],
            [p01 - 0.085, 0, 1], [0.5, 0, 1],
            [p00 - 0.085, 1, 0], [0.5, 1, 0],
            [p01 - 0.085, 1, 1], [0.5, 1, 1]
        ],
        'location': 'leading_edge',
        'side': 'suction'
    }
    
    blocks['le_1'] = {
        'corners': [
            [0.5, 0, 0], [0.085 + p10, 0, 0],
            [0.5, 0, 1], [0.085 + p11, 0, 1],
            [0.5, 1, 0], [0.085 + p10, 1, 0],
            [0.5, 1, 1], [0.085 + p11, 1, 1]
        ],
        'location': 'leading_edge',
        'side': 'pressure'
    }
    
    return blocks

def map_to_physical(parametric_coords, r_hub=0.5, r_shroud=1.9, z_length=2.5):
    """Map parametric [u,v,w] to physical [x,y,z]."""
    u, v, w = parametric_coords
    
    # Radial: linear from hub to shroud
    r = r_hub + u * (r_shroud - r_hub)
    
    # Angular: v * 90 degrees (4 blades = 4 passages)
    theta = v * (np.pi / 2)
    
    # Axial
    z = w * z_length
    
    # Cartesian
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    
    return [x, y, z]

def rotate_90deg(point, n_rotations=1):
    """Rotate point by n_rotations * 90 degrees around Z."""
    angle = n_rotations * np.pi / 2
    x, y, z = point
    return [
        x * np.cos(angle) - y * np.sin(angle),
        x * np.sin(angle) + y * np.cos(angle),
        z
    ]

def generate_all_instances(base_blocks, params):
    """Generate all 28 instances (7 blocks × 4 rotations)."""
    r_hub = params.get('cV_r_hub', 0.5)
    r_shroud = params.get('cV_r_shroud', 1.9)
    z_length = params.get('cV_l_ru', 2.5)
    
    all_instances = []
    
    for block_name, block_data in base_blocks.items():
        for rot in range(4):
            instance_name = f"{block_name}_rot{rot*90}"
            
            # Map corners to physical and rotate
            physical_corners = []
            for p_corner in block_data['corners']:
                phys = map_to_physical(p_corner, r_hub, r_shroud, z_length)
                rotated = rotate_90deg(phys, rot)
                physical_corners.append(rotated)
            
            # Compute center
            center = np.mean(physical_corners, axis=0).tolist()
            
            all_instances.append({
                'name': instance_name,
                'base_block': block_name,
                'rotation': rot * 90,
                'corners': physical_corners,
                'center': center,
                'location': block_data['location'],
                'side': block_data['side']
            })
    
    return all_instances

def create_3d_plot(instances):
    """Create interactive Plotly 3D visualization."""
    fig = go.Figure()
    
    # Color scheme by location
    colors = {
        'leading_edge': '#FF6B6B',  # Red
        'main': '#4ECDC4',           # Teal
        'trailing_edge': '#45B7D1'   # Blue
    }
    
    # Edge pairs for hexahedron wireframe
    edges = [
        (0,1), (1,3), (3,2), (2,0),  # Bottom face
        (4,5), (5,7), (7,6), (6,4),  # Top face
        (0,4), (1,5), (2,6), (3,7)   # Vertical edges
    ]
    
    for inst in instances:
        corners = np.array(inst['corners'])
        color = colors.get(inst['location'], '#999999')
        
        # Draw wireframe
        for edge in edges:
            i, j = edge
            fig.add_trace(go.Scatter3d(
                x=[corners[i,0], corners[j,0]],
                y=[corners[i,1], corners[j,1]],
                z=[corners[i,2], corners[j,2]],
                mode='lines',
                line=dict(color=color, width=1),
                showlegend=False,
                hoverinfo='skip'
            ))
        
        # Draw corner points
        fig.add_trace(go.Scatter3d(
            x=corners[:,0],
            y=corners[:,1],
            z=corners[:,2],
            mode='markers',
            marker=dict(size=3, color=color),
            name=inst['name'],
            text=[f"{inst['name']}_c{i}" for i in range(8)],
            hovertemplate='%{text}<br>x: %{x:.3f}<br>y: %{y:.3f}<br>z: %{z:.3f}<extra></extra>'
        ))
        
        # Label center
        center = inst['center']
        fig.add_trace(go.Scatter3d(
            x=[center[0]],
            y=[center[1]],
            z=[center[2]],
            mode='markers+text',
            marker=dict(size=5, color='black'),
            text=[inst['name']],
            textposition='top center',
            textfont=dict(size=8),
            showlegend=False,
            hoverinfo='skip'
        ))
    
    # Layout
    fig.update_layout(
        title='T1_9 Turbine: 28 Hexa Block Instances (7 blocks × 4 rotations)',
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z',
            aspectmode='data'
        ),
        width=1200,
        height=800,
        showlegend=False
    )
    
    return fig

def main():
    # Parameters from T2_7461 (same parametric structure for T1_9)
    params = {
        'cV_ru_divideInternalMeshBlock_0_0': 0.48500001430511475,
        'cV_ru_divideInternalMeshBlock_0_1': 0.48500001430511475,
        'cV_ru_divideInternalMeshBlock_1_0': 0.5149999856948853,
        'cV_ru_divideInternalMeshBlock_1_1': 0.5149999856948853,
        'cV_r_hub': 0.5,
        'cV_r_shroud': 1.9,  # T1_9 uses 1.9 instead of 2.0
        'cV_l_ru': 2.5,
        'cV_ru_nBlades': 4
    }
    
    print("Computing base blocks...")
    base_blocks = compute_base_blocks(params)
    
    print(f"Found {len(base_blocks)} base block types")
    for name, data in base_blocks.items():
        print(f"  {name}: {data['location']} - {data['side']}")
    
    print("\nGenerating 28 instances (7 blocks × 4 rotations)...")
    instances = generate_all_instances(base_blocks, params)
    
    # Save to JSON
    output = {
        'metadata': {
            'n_base_blocks': len(base_blocks),
            'n_rotations': 4,
            'total_instances': len(instances),
            'r_hub': params['cV_r_hub'],
            'r_shroud': params['cV_r_shroud'],
            'z_length': params['cV_l_ru']
        },
        'instances': instances
    }
    
    with open('block_instances_28.json', 'w') as f:
        json.dump(output, f, indent=2, default=lambda x: float(x) if isinstance(x, np.floating) else x)
    
    print(f"Saved to block_instances_28.json")
    
    # Create 3D plot
    print("\nCreating interactive 3D plot...")
    fig = create_3d_plot(instances)
    
    # Save as HTML
    html_path = Path('block_visualization_28.html')
    fig.write_html(html_path)
    print(f"Saved interactive plot to {html_path}")
    
    # Also create a summary of corner candidates for singularities
    print("\nExtracting corner candidates for singularities...")
    corner_candidates = []
    
    for inst in instances:
        corners = inst['corners']
        # The 8 corners of each block
        for i, corner in enumerate(corners):
            corner_candidates.append({
                'instance': inst['name'],
                'corner_id': i,
                'position': corner,
                'location': inst['location'],
                'side': inst['side']
            })
    
    # Group by unique positions (within tolerance)
    unique_positions = []
    tolerance = 0.01
    
    for candidate in corner_candidates:
        pos = np.array(candidate['position'])
        found = False
        
        for unique in unique_positions:
            if np.linalg.norm(pos - np.array(unique['position'])) < tolerance:
                unique['instances'].append(candidate['instance'])
                unique['corner_ids'].append(candidate['corner_id'])
                found = True
                break
        
        if not found:
            unique_positions.append({
                'position': pos.tolist(),
                'instances': [candidate['instance']],
                'corner_ids': [candidate['corner_id']],
                'location': candidate['location'],
                'side': candidate['side']
            })
    
    print(f"Found {len(unique_positions)} unique corner positions (from {len(corner_candidates)} total corners)")
    
    # Save corner candidates
    with open('corner_candidates.json', 'w') as f:
        json.dump({
            'n_unique_corners': len(unique_positions),
            'corners': unique_positions
        }, f, indent=2)
    
    print(f"Saved corner candidates to corner_candidates.json")
    print("\nDone! Open 'block_visualization_28.html' in a browser to explore the blocks.")
    print("The file 'corner_candidates.json' contains all unique corner positions.")

if __name__ == '__main__':
    main()
