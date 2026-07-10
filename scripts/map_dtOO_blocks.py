#!/usr/bin/env python3
"""
Compute explicit block corners for dtOO T2_7461 turbine mesh.

This script:
1. Reads parameter values from T2_7461.xml
2. Computes block corners in parametric [u,v,w] space
3. Maps them to physical [x,y,z] coordinates using the channel geometry
4. Outputs JSON with both parametric and physical coordinates
"""

import json
import numpy as np
from pathlib import Path

def compute_block_corners(params, context='ru'):
    """Compute all block corners for a given context (ru, ruWithRoundingBlade, etc.)."""
    
    # Get parameters for this context
    p00 = params[f'cV_{context}_divideInternalMeshBlock_0_0']
    p01 = params[f'cV_{context}_divideInternalMeshBlock_0_1']
    p10 = params[f'cV_{context}_divideInternalMeshBlock_1_0']
    p11 = params[f'cV_{context}_divideInternalMeshBlock_1_1']
    
    # Base offset depends on context
    if context == 'ru':
        base = 0.495
    elif context in ['ruWithRoundingBlade', 'ruWithRoundingRounding']:
        base = 0.43
    else:
        base = 0.495
    
    blocks = {}
    
    # Trailing Edge blocks
    blocks[f'{context}_te_0'] = {
        'parametric': [
            [base + p00, 0, 0], [1.0, 0, 0],
            [base + p01, 0, 1], [1.0, 0, 1],
            [base + p00, 1, 0], [1.0, 1, 0],
            [base + p01, 1, 1], [1.0, 1, 1]
        ],
        'location': 'trailing_edge',
        'side': 'pressure'  # x > blade
    }
    
    blocks[f'{context}_te_1'] = {
        'parametric': [
            [0, 0, 0], [p10 - base, 0, 0],
            [0, 0, 1], [p11 - base, 0, 1],
            [0, 1, 0], [p10 - base, 1, 0],
            [0, 1, 1], [p11 - base, 1, 1]
        ],
        'location': 'trailing_edge',
        'side': 'suction'  # x < blade
    }
    
    # Main blocks around blade
    if context == 'ru':
        blocks[f'{context}_0_0'] = {
            'parametric': [
                [p10 - 0.43, 0, 0], [p00 - 0.025, 0, 0],
                [p11 - 0.43, 0, 1], [p01 - 0.025, 0, 1],
                [p10 - 0.43, 1, 0], [p00 - 0.025, 1, 0],
                [p11 - 0.43, 1, 1], [p01 - 0.025, 1, 1]
            ],
            'location': 'main',
            'side': 'suction'
        }
        
        blocks[f'{context}_0_1'] = {
            'parametric': [
                [p00 - 0.025, 0, 0], [0.025 + p10, 0, 0],
                [p01 - 0.025, 0, 1], [0.025 + p11, 0, 1],
                [p00 - 0.025, 1, 0], [0.025 + p10, 1, 0],
                [p01 - 0.025, 1, 1], [0.025 + p11, 1, 1]
            ],
            'location': 'main',
            'side': 'center'
        }
        
        blocks[f'{context}_0_2'] = {
            'parametric': [
                [0.085 + p10, 0, 0], [0.43 + p00, 0, 0],
                [0.085 + p11, 0, 1], [0.43 + p01, 0, 1],
                [0.085 + p10, 1, 0], [0.43 + p00, 1, 0],
                [0.085 + p11, 1, 1], [0.43 + p01, 1, 1]
            ],
            'location': 'main',
            'side': 'pressure'
        }
        
        # Leading Edge blocks
        blocks[f'{context}_le_0'] = {
            'parametric': [
                [p00 - 0.085, 0, 0], [0.5, 0, 0],
                [p01 - 0.085, 0, 1], [0.5, 0, 1],
                [p00 - 0.085, 1, 0], [0.5, 1, 0],
                [p01 - 0.085, 1, 1], [0.5, 1, 1]
            ],
            'location': 'leading_edge',
            'side': 'suction'
        }
        
        blocks[f'{context}_le_1'] = {
            'parametric': [
                [0.5, 0, 0], [0.085 + p10, 0, 0],
                [0.5, 0, 1], [0.085 + p11, 0, 1],
                [0.5, 1, 0], [0.085 + p10, 1, 0],
                [0.5, 1, 1], [0.085 + p11, 1, 1]
            ],
            'location': 'leading_edge',
            'side': 'pressure'
        }
    
    return blocks

def map_to_physical(parametric_coords, r_hub=0.5, r_shroud=2.0, z_length=2.5):
    """
    Map parametric [u,v,w] coordinates to physical [x,y,z].
    
    The dtOO channel geometry is a surface of revolution:
    - u maps to radial position (r_hub to r_shroud)
    - v maps to angular position (0 to 2*pi/nBlades)
    - w maps to axial position (0 to z_length)
    """
    u, v, w = parametric_coords
    
    # Radial position: linear interpolation from hub to shroud
    r = r_hub + u * (r_shroud - r_hub)
    
    # Angular position: v * (2*pi / nBlades)
    # For a 4-blade turbine, each blade passage is 90 degrees
    theta = v * (2 * np.pi / 4)  # 4 blades
    
    # Axial position
    z = w * z_length
    
    # Convert to Cartesian
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    
    return [x, y, z]

def main():
    # Parameters from T2_7461.xml
    params = {
        'cV_ru_divideInternalMeshBlock_0_0': 0.48500001430511475,
        'cV_ru_divideInternalMeshBlock_0_1': 0.48500001430511475,
        'cV_ru_divideInternalMeshBlock_1_0': 0.5149999856948853,
        'cV_ru_divideInternalMeshBlock_1_1': 0.5149999856948853,
        'cV_ruWithRoundingBlade_divideInternalMeshBlock_0_0': 0.48500001430511475,
        'cV_ruWithRoundingBlade_divideInternalMeshBlock_0_1': 0.48500001430511475,
        'cV_ruWithRoundingBlade_divideInternalMeshBlock_1_0': 0.5149999856948853,
        'cV_ruWithRoundingBlade_divideInternalMeshBlock_1_1': 0.5149999856948853,
        'cV_r_hub': 0.5,
        'cV_r_shroud': 2.0,
        'cV_l_ru': 2.5,
        'cV_ru_nBlades': 4
    }
    
    # Compute blocks for ru context (no rounding)
    blocks = compute_block_corners(params, 'ru')
    
    # Map to physical coordinates
    r_hub = params['cV_r_hub']
    r_shroud = params['cV_r_shroud']
    z_length = params['cV_l_ru']
    
    print("Block corners in parametric [u,v,w] and physical [x,y,z] space:")
    print("=" * 80)
    
    for block_name, block_data in blocks.items():
        print(f"\n{block_name}:")
        print(f"  Location: {block_data['location']}, Side: {block_data['side']}")
        print(f"  Corners:")
        
        physical_corners = []
        for i, p_coord in enumerate(block_data['parametric']):
            phys = map_to_physical(p_coord, r_hub, r_shroud, z_length)
            physical_corners.append(phys)
            print(f"    {i}: param=({p_coord[0]:.4f}, {p_coord[1]:.4f}, {p_coord[2]:.4f}) -> "
                  f"phys=({phys[0]:.4f}, {phys[1]:.4f}, {phys[2]:.4f})")
        
        block_data['physical'] = physical_corners
        
        # Compute center
        center = np.mean(physical_corners, axis=0)
        print(f"  Center: ({center[0]:.4f}, {center[1]:.4f}, {center[2]:.4f})")
    
    # Save to JSON
    output = {
        'metadata': {
            'source': 'dtOO T2_7461',
            'context': 'ru (no rounding)',
            'r_hub': r_hub,
            'r_shroud': r_shroud,
            'z_length': z_length,
            'n_blades': 4
        },
        'blocks': blocks
    }
    
    output_path = Path('dtOO_blocks_mapped.json')
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2, default=lambda x: float(x) if isinstance(x, np.float64) else x)
    
    print(f"\n\nSaved to {output_path}")
    print(f"Total blocks: {len(blocks)}")

if __name__ == '__main__':
    main()
