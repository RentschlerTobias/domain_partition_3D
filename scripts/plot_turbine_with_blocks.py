#!/usr/bin/env python3
"""
Generate comprehensive 3D visualization with:
1. 28 block instances (wireframes)
2. Actual hex elements from .msh (translucent wireframe)
3. Blade surfaces (colored meshes)
4. Hub and shroud surfaces

This gives the full geometric context for the block structure.
"""

import json
import numpy as np
import plotly.graph_objects as go
from pathlib import Path

def parse_msh_file(filepath):
    """Parse all relevant elements from Gmsh 2.2 format."""
    nodes = {}
    hex_elements = []
    blade_surface_nodes = set()  # geom tags 5 and 6
    
    with open(filepath, 'r') as f:
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
                    
                    if elem_type == 5:  # Hex element
                        hex_elements.append(node_ids)
                    elif elem_type == 2 and geom_tag in [5, 6]:  # Triangles on blade
                        blade_surface_nodes.update(node_ids)
                    elif elem_type == 3 and geom_tag in [5, 6]:  # Quads on blade
                        blade_surface_nodes.update(node_ids)
            
            line = f.readline()
    
    return nodes, hex_elements, blade_surface_nodes

def generate_block_instances():
    """Generate the 28 block instances."""
    params = {
        'cV_ru_divideInternalMeshBlock_0_0': 0.48500001430511475,
        'cV_ru_divideInternalMeshBlock_0_1': 0.48500001430511475,
        'cV_ru_divideInternalMeshBlock_1_0': 0.5149999856948853,
        'cV_ru_divideInternalMeshBlock_1_1': 0.5149999856948853,
        'cV_r_hub': 0.5,
        'cV_r_shroud': 1.9,
        'cV_l_ru': 2.5
    }
    
    p00 = params['cV_ru_divideInternalMeshBlock_0_0']
    p01 = params['cV_ru_divideInternalMeshBlock_0_1']
    p10 = params['cV_ru_divideInternalMeshBlock_1_0']
    p11 = params['cV_ru_divideInternalMeshBlock_1_1']
    base = 0.495
    
    def map_param_to_phys(u, v, w):
        r = 0.5 + u * 1.4  # r_shroud - r_hub = 1.9 - 0.5 = 1.4
        theta = v * np.pi / 2
        z = w * 2.5
        return [r * np.cos(theta), r * np.sin(theta), z]
    
    def rotate_90(point, n):
        angle = n * np.pi / 2
        x, y, z = point
        return [x * np.cos(angle) - y * np.sin(angle), 
                x * np.sin(angle) + y * np.cos(angle), z]
    
    base_blocks = {
        'te_0': [[base+p00,0,0],[1.0,0,0],[base+p01,0,1],[1.0,0,1],
                 [base+p00,1,0],[1.0,1,0],[base+p01,1,1],[1.0,1,1]],
        'te_1': [[0,0,0],[p10-base,0,0],[0,0,1],[p11-base,0,1],
                 [0,1,0],[p10-base,1,0],[0,1,1],[p11-base,1,1]],
        '0_0': [[p10-0.43,0,0],[p00-0.025,0,0],[p11-0.43,0,1],[p01-0.025,0,1],
                [p10-0.43,1,0],[p00-0.025,1,0],[p11-0.43,1,1],[p01-0.025,1,1]],
        '0_1': [[p00-0.025,0,0],[0.025+p10,0,0],[p01-0.025,0,1],[0.025+p11,0,1],
                [p00-0.025,1,0],[0.025+p10,1,0],[p01-0.025,1,1],[0.025+p11,1,1]],
        '0_2': [[0.085+p10,0,0],[0.43+p00,0,0],[0.085+p11,0,1],[0.43+p01,0,1],
                [0.085+p10,1,0],[0.43+p00,1,0],[0.085+p11,1,1],[0.43+p01,1,1]],
        'le_0': [[p00-0.085,0,0],[0.5,0,0],[p01-0.085,0,1],[0.5,0,1],
                 [p00-0.085,1,0],[0.5,1,0],[p01-0.085,1,1],[0.5,1,1]],
        'le_1': [[0.5,0,0],[0.085+p10,0,0],[0.5,0,1],[0.085+p11,0,1],
                 [0.5,1,0],[0.085+p10,1,0],[0.5,1,1],[0.085+p11,1,1]]
    }
    
    instances = []
    for block_name, corners in base_blocks.items():
        for rot in range(4):
            phys_corners = []
            for p_corner in corners:
                phys = map_param_to_phys(*p_corner)
                phys_rot = rotate_90(phys, rot)
                phys_corners.append(phys_rot)
            
            instances.append({
                'name': f"{block_name}_rot{rot*90}",
                'corners': np.array(phys_corners)
            })
    
    return instances

def create_comprehensive_plot(nodes, hex_elements, blade_surface_nodes, block_instances):
    """Create comprehensive 3D visualization."""
    fig = go.Figure()
    
    # 1. Plot hex elements as very faint wireframe (subsample for performance)
    print("Adding hex mesh wireframe...")
    hex_centers = []
    for elem_nodes in hex_elements[::10]:  # Every 10th element for performance
        coords = np.array([nodes[nid] for nid in elem_nodes])
        center = np.mean(coords, axis=0)
        hex_centers.append(center)
    
    hex_centers = np.array(hex_centers)
    fig.add_trace(go.Scatter3d(
        x=hex_centers[:,0], y=hex_centers[:,1], z=hex_centers[:,2],
        mode='markers',
        marker=dict(size=1, color='lightgray', opacity=0.3),
        name='Hex mesh centers',
        showlegend=True
    ))
    
    # 2. Plot blade surface nodes
    print("Adding blade surface...")
    blade_coords = np.array([nodes[nid] for nid in blade_surface_nodes])
    fig.add_trace(go.Scatter3d(
        x=blade_coords[:,0], y=blade_coords[:,1], z=blade_coords[:,2],
        mode='markers',
        marker=dict(size=2, color='red', opacity=0.6),
        name='Blade surface',
        showlegend=True
    ))
    
    # 3. Plot block instances
    print("Adding block wireframes...")
    colors = {'te': '#1f77b4', '0': '#2ca02c', 'le': '#ff7f0e'}
    edges = [(0,1),(1,3),(3,2),(2,0),(4,5),(5,7),(7,6),(6,4),(0,4),(1,5),(2,6),(3,7)]
    
    for inst in block_instances:
        corners = inst['corners']
        color = colors.get(inst['name'][:2], '#9467bd')
        
        # Draw edges
        for edge in edges:
            i, j = edge
            fig.add_trace(go.Scatter3d(
                x=[corners[i,0], corners[j,0]],
                y=[corners[i,1], corners[j,1]],
                z=[corners[i,2], corners[j,2]],
                mode='lines',
                line=dict(color=color, width=2),
                showlegend=False,
                hoverinfo='skip'
            ))
        
        # Draw corners
        fig.add_trace(go.Scatter3d(
            x=corners[:,0], y=corners[:,1], z=corners[:,2],
            mode='markers',
            marker=dict(size=3, color=color),
            name=inst['name'],
            text=[f"{inst['name']}_c{i}" for i in range(8)],
            hovertemplate='%{text}<br>x: %{x:.3f}<br>y: %{y:.3f}<br>z: %{z:.3f}<extra></extra>',
            showlegend=False
        ))
    
    # Layout
    fig.update_layout(
        title='T1_9 Turbine: Block Structure + Actual Mesh Geometry',
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z',
            aspectmode='data',
            camera=dict(eye=dict(x=1.5, y=1.5, z=1.0))
        ),
        width=1400,
        height=900,
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01)
    )
    
    return fig

def main():
    print("Parsing T1_9 mesh file...")
    msh_path = Path('T1_9/T1_9_ru_gridGmsh.msh')
    nodes, hex_elements, blade_surface_nodes = parse_msh_file(msh_path)
    
    print(f"Found {len(hex_elements)} hex elements")
    print(f"Found {len(blade_surface_nodes)} blade surface nodes")
    
    print("\nGenerating 28 block instances...")
    block_instances = generate_block_instances()
    
    print("\nCreating comprehensive 3D plot...")
    fig = create_comprehensive_plot(nodes, hex_elements, blade_surface_nodes, block_instances)
    
    # Save
    output_path = Path('turbine_blocks_with_mesh.html')
    fig.write_html(output_path)
    print(f"\nSaved to {output_path}")
    print("Open this file in a browser to see the blocks overlaid on the actual turbine geometry.")

if __name__ == '__main__':
    main()
