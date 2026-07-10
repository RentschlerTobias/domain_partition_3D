#!/usr/bin/env python3
"""
Parse dtOO XML block definitions and extract explicit hexa block corner coordinates.

This script reads:
1. T2_7461.xml (or similar state file) for parameter values
2. ru_bladeRunnerMeshBlock.xml for block corner definitions

And outputs:
- JSON with all block corners in parametric space [u,v,w]
- The blocks can then be mapped to physical space using the channel geometry
"""

import xml.etree.ElementTree as ET
import json
import re
from pathlib import Path

def parse_parameter_file(filepath):
    """Extract all cV_* parameters from state file."""
    tree = ET.parse(filepath)
    root = tree.getroot()
    
    params = {}
    for const in root.iter('constValue'):
        label = const.get('label', '')
        value = const.get('value', '')
        if label.startswith('cV_'):
            # Evaluate simple expressions
            try:
                # Handle expressions like "2*(#cV_r_shroud#-#cV_r_hub#)*#cV_Cmu#^(-1/4)"
                # First, substitute known params
                val = str(value)
                # Find all #cV_*# references
                refs = re.findall(r'#(cV_[^#]+)#', val)
                for ref in refs:
                    if ref in params:
                        val = val.replace(f'#{ref}#', str(params[ref]))
                
                # Evaluate if it's a simple number
                if val.replace('.', '').replace('-', '').isdigit():
                    params[label] = float(val)
                else:
                    # Try to evaluate expression
                    try:
                        # Replace ^ with ** for Python
                        val = val.replace('^', '**')
                        # Replace _pi with math.pi
                        val = val.replace('_pi', str(3.141592653589793))
                        params[label] = eval(val)
                    except:
                        params[label] = val  # Keep as string if can't evaluate
            except:
                params[label] = value
    
    return params

def parse_block_file(filepath, params):
    """Parse block definitions and compute corner coordinates."""
    tree = ET.parse(filepath)
    root = tree.getroot()
    
    blocks = {}
    
    # Find all vec3dTriLinearThreeD builders (these define hexa blocks)
    for func in root.iter('function'):
        func_label = func.get('label', '')
        if 'meshBlock_closed' in func_label:
            builder = func.find('builder')
            if builder is not None and builder.get('name') == 'vec3dTriLinearThreeD':
                # Extract 8 corner points
                corners = []
                for point in builder.iter('Point_3'):
                    x_expr = point.get('x', '0')
                    y_expr = point.get('y', '0')
                    z_expr = point.get('z', '0')
                    
                    # Substitute parameters
                    x_val = evaluate_expression(x_expr, params)
                    y_val = evaluate_expression(y_expr, params)
                    z_val = evaluate_expression(z_expr, params)
                    
                    corners.append([x_val, y_val, z_val])
                
                if len(corners) == 8:
                    blocks[func_label] = {
                        'corners': corners,
                        'type': 'triLinearHexa'
                    }
    
    return blocks

def evaluate_expression(expr, params):
    """Evaluate an expression with parameter substitution."""
    # Substitute parameters: #cV_*# -> value
    val = str(expr)
    refs = re.findall(r'#(cV_[^#]+)#', val)
    for ref in refs:
        if ref in params:
            val = val.replace(f'#{ref}#', str(params[ref]))
    
    # Try to evaluate
    try:
        return float(val)
    except:
        # Return as string if can't evaluate
        return val

def main():
    repo_path = Path('/tmp/dtOO_repo/test/tistos')
    
    # Parse parameter file
    print("Parsing T2_7461.xml for parameters...")
    params = parse_parameter_file(repo_path / 'T2_7461.xml')
    print(f"Found {len(params)} parameters")
    
    # Show key parameters
    key_params = [
        'cV_ru_divideInternalMeshBlock_0_0',
        'cV_ru_divideInternalMeshBlock_0_1',
        'cV_ru_divideInternalMeshBlock_1_0',
        'cV_ru_divideInternalMeshBlock_1_1',
        'cV_r_meshBlockThickness',
        'cV_r_hub',
        'cV_r_shroud',
        'cV_l_ru'
    ]
    
    print("\nKey block parameters:")
    for kp in key_params:
        if kp in params:
            print(f"  {kp} = {params[kp]}")
    
    # Parse block definitions
    print("\nParsing ru_bladeRunnerMeshBlock.xml for blocks...")
    blocks = parse_block_file(repo_path / 'xml/ru_bladeRunnerMeshBlock.xml', params)
    print(f"Found {len(blocks)} block definitions")
    
    # Print block details
    print("\nBlock details:")
    for block_name, block_data in blocks.items():
        corners = block_data['corners']
        print(f"\n{block_name}:")
        print(f"  Type: {block_data['type']}")
        print(f"  Corners (u,v,w in parametric space):")
        for i, c in enumerate(corners):
            vals = [f"{v:.4f}" if isinstance(v, (int, float)) else str(v) for v in c]
            print(f"    {i}: ({vals[0]}, {vals[1]}, {vals[2]})")
        
        # Compute block dimensions
        numeric_corners = [c for c in corners if all(isinstance(v, (int, float)) for v in c)]
        if numeric_corners:
            u_range = [min(c[0] for c in numeric_corners), max(c[0] for c in numeric_corners)]
            v_range = [min(c[1] for c in numeric_corners), max(c[1] for c in numeric_corners)]
            w_range = [min(c[2] for c in numeric_corners), max(c[2] for c in numeric_corners)]
            print(f"  Dimensions: u={u_range[1]-u_range[0]:.4f}, v={v_range[1]-v_range[0]:.4f}, w={w_range[1]-w_range[0]:.4f}")
        else:
            print("  Dimensions: (contains unevaluated expressions)")
    
    # Save to JSON
    output = {
        'parameters': params,
        'blocks': blocks,
        'metadata': {
            'source': 'dtOO T2_7461',
            'parameter_file': 'T2_7461.xml',
            'block_file': 'ru_bladeRunnerMeshBlock.xml',
            'note': 'Coordinates are in parametric space [u,v,w] before mapping to physical space'
        }
    }
    
    output_path = Path('dtOO_blocks.json')
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"\n\nSaved to {output_path}")
    print(f"Total blocks: {len(blocks)}")
    print("\nNote: These are parametric coordinates. To get physical coordinates,")
    print("they must be mapped through the channel geometry (rM2dTo3d_ru_channel).")

if __name__ == '__main__':
    main()
