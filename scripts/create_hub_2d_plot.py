#!/usr/bin/env python3
"""
Extract 2D boundary edges on hub surface (z ~ 1.06-1.3, lower hex layer)
and create 2D plot with labeled corners.
"""

import numpy as np
import json
from collections import defaultdict

def parse_hex_regions(msh_file):
    nodes = {}
    hex_elements = []
    
    with open(msh_file, 'r') as f:
        lines = f.readlines()
    
    in_nodes = False
    for line in lines:
        line = line.strip()
        if line == '$Nodes':
            in_nodes = True
            continue
        if line == '$EndNodes':
            break
        if in_nodes:
            if line.isdigit():
                continue
            parts = line.split()
            tag = int(parts[0])
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            nodes[tag] = np.array([x, y, z])
    
    in_elements = False
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line == '$Elements':
            in_elements = True
            continue
        if line == '$EndElements':
            break
        if in_elements:
            if line.isdigit():
                continue
            parts = line.split()
            elem_type = int(parts[1])
            if elem_type == 5:
                tag = int(parts[0])
                num_tags = int(parts[2])
                geom_tag = int(parts[4]) if num_tags >= 2 else 0
                nodes_idx = [int(p) for p in parts[3 + num_tags:]]
                hex_elements.append({
                    'tag': tag,
                    'geom_tag': geom_tag,
                    'nodes': nodes_idx
                })
    
    return nodes, hex_elements

def find_lower_boundary_edges(nodes, hex_elements, z_max=1.3):
    """Find boundary edges on lower hex layer (near hub)."""
    
    # Build face-to-regions mapping
    face_to_regions = defaultdict(set)
    face_info = {}
    
    for elem in hex_elements:
        n = elem['nodes']
        faces = [
            frozenset([n[0], n[1], n[2], n[3]]),
            frozenset([n[4], n[5], n[6], n[7]]),
            frozenset([n[0], n[1], n[5], n[4]]),
            frozenset([n[2], n[3], n[7], n[6]]),
            frozenset([n[0], n[3], n[7], n[4]]),
            frozenset([n[1], n[2], n[6], n[5]]),
        ]
        for face in faces:
            face_to_regions[face].add(elem['geom_tag'])
            if face not in face_info:
                # Compute average z
                zs = [nodes[node][2] for node in face]
                face_info[face] = {
                    'avg_z': np.mean(zs),
                    'min_z': min(zs),
                    'max_z': max(zs)
                }
    
    # Find inter-region faces that are near hub (low z)
    boundary_edges = set()
    for face, regions in face_to_regions.items():
        if len(regions) >= 2:
            info = face_info[face]
            # Keep only faces with average z < z_max (lower half near hub)
            if info['avg_z'] < z_max:
                face_nodes = list(face)
                for i in range(len(face_nodes)):
                    n1, n2 = face_nodes[i], face_nodes[(i+1) % len(face_nodes)]
                    boundary_edges.add(tuple(sorted([n1, n2])))
    
    return boundary_edges

def find_block_corners_2d(nodes, hex_elements, z_max=1.3):
    """Find 2D corners by projecting block BB corners to hub plane."""
    
    regions = defaultdict(list)
    for elem in hex_elements:
        regions[elem['geom_tag']].append(elem)
    
    corners = []
    
    for geom_tag in sorted(regions.keys()):
        elems = regions[geom_tag]
        region_nodes = set()
        for elem in elems:
            region_nodes.update(elem['nodes'])
        
        # Filter nodes near hub (z < z_max)
        lower_nodes = [n for n in region_nodes if nodes[n][2] < z_max]
        
        if len(lower_nodes) < 4:
            continue
        
        coords = np.array([nodes[n] for n in lower_nodes])
        
        # Find bounding box in 2D (x,y)
        bb_min = coords.min(axis=0)
        bb_max = coords.max(axis=0)
        
        # 4 corners of 2D BB (bottom-left, bottom-right, top-left, top-right)
        bb_corners_2d = [
            (bb_min[0], bb_min[1]),
            (bb_max[0], bb_min[1]),
            (bb_min[0], bb_max[1]),
            (bb_max[0], bb_max[1])
        ]
        
        region_corners = []
        used_nodes = set()
        
        for target_x, target_y in bb_corners_2d:
            best_dist = float('inf')
            best_node = None
            
            for n in lower_nodes:
                if n in used_nodes:
                    continue
                x, y = nodes[n][0], nodes[n][1]
                dist = np.sqrt((x - target_x)**2 + (y - target_y)**2)
                if dist < best_dist:
                    best_dist = dist
                    best_node = n
            
            if best_node is not None:
                used_nodes.add(best_node)
                pos = nodes[best_node]
                region_corners.append({
                    'node': best_node,
                    'region': geom_tag,
                    'pos': [float(pos[0]), float(pos[1]), float(pos[2])]
                })
        
        corners.extend(region_corners)
    
    return corners

def create_2d_hub_plot(nodes, edges, corners, output_file):
    """Create 2D HTML plot with Canvas."""
    
    # Convert edges to coordinate pairs
    edge_coords = []
    for n1, n2 in edges:
        p1 = nodes[n1]
        p2 = nodes[n2]
        edge_coords.append({
            'x1': float(p1[0]), 'y1': float(p1[1]),
            'x2': float(p2[0]), 'y2': float(p2[1])
        })
    
    # Get all points for bounding box
    all_x = []
    all_y = []
    for ec in edge_coords:
        all_x.extend([ec['x1'], ec['x2']])
        all_y.extend([ec['y1'], ec['y2']])
    
    margin = 0.1
    min_x, max_x = min(all_x) - margin, max(all_x) + margin
    min_y, max_y = min(all_y) - margin, max(all_y) + margin
    
    width = max_x - min_x
    height = max_y - min_y
    
    # Scale factor for canvas
    canvas_size = 800
    scale = canvas_size / max(width, height)
    
    # Convert to canvas coordinates
    def to_canvas(x, y):
        cx = (x - min_x) * scale + 50
        cy = canvas_size - (y - min_y) * scale - 50
        return cx, cy
    
    # Corner data
    region_colors = {2: '#4488ff', 3: '#44ff88', 4: '#ff8844', 5: '#ff44ff', 6: '#ffff44'}
    
    corners_js = []
    for i, corner in enumerate(corners):
        color = region_colors.get(corner['region'], '#ffffff')
        cx, cy = to_canvas(corner['pos'][0], corner['pos'][1])
        corners_js.append({
            'id': i,
            'node': corner['node'],
            'region': corner['region'],
            'x': corner['pos'][0],
            'y': corner['pos'][1],
            'z': corner['pos'][2],
            'cx': cx,
            'cy': cy,
            'color': color
        })
    
    # Edge data in canvas coordinates
    edges_js = []
    for ec in edge_coords:
        x1, y1 = to_canvas(ec['x1'], ec['y1'])
        x2, y2 = to_canvas(ec['x2'], ec['y2'])
        edges_js.append({'x1': x1, 'y1': y1, 'x2': x2, 'y2': y2})
    
    html = f'''<!DOCTYPE html>
<html>
<head>
    <title>Hub Surface Block Boundaries (2D)</title>
    <style>
        body {{
            margin: 0;
            font-family: 'Segoe UI', Arial, sans-serif;
            background: #0a0a1a;
            color: white;
            display: flex;
            flex-direction: column;
            align-items: center;
            padding: 20px;
        }}
        h2 {{
            margin: 0 0 10px 0;
            color: #44ff88;
            font-size: 20px;
        }}
        #info {{
            margin-bottom: 15px;
            font-size: 13px;
            text-align: center;
            opacity: 0.8;
        }}
        #canvas-container {{
            position: relative;
            border: 2px solid rgba(255,255,255,0.2);
            border-radius: 8px;
            background: #111122;
        }}
        canvas {{
            display: block;
            border-radius: 8px;
        }}
        #legend {{
            display: flex;
            gap: 15px;
            margin-top: 15px;
            font-size: 12px;
        }}
        .legend-item {{
            display: flex;
            align-items: center;
            gap: 5px;
        }}
        .color-box {{
            width: 15px;
            height: 15px;
            border-radius: 3px;
            border: 1px solid white;
        }}
        #corner-list {{
            position: absolute;
            top: 10px;
            right: 10px;
            background: rgba(0,0,0,0.8);
            padding: 10px;
            border-radius: 6px;
            font-size: 10px;
            max-height: 90%;
            overflow-y: auto;
            max-width: 150px;
        }}
        .corner-item {{
            margin: 2px 0;
            padding: 2px 4px;
            cursor: pointer;
            border-left: 2px solid;
            border-radius: 2px;
        }}
        .corner-item:hover {{
            background: rgba(255,255,255,0.2);
        }}
    </style>
</head>
<body>
    <h2>Hub Surface Block Boundaries (2D Projection)</h2>
    <div id="info">
        Lower hex layer boundaries (z &lt; 1.3) | {len(edges)} edges | {len(corners)} corners
    </div>
    
    <div id="canvas-container">
        <canvas id="plot" width="900" height="900"></canvas>
        <div id="corner-list"></div>
    </div>
    
    <div id="legend">
        <div class="legend-item"><div class="color-box" style="background:#4488ff"></div>Region 2</div>
        <div class="legend-item"><div class="color-box" style="background:#44ff88"></div>Region 3</div>
        <div class="legend-item"><div class="color-box" style="background:#ff8844"></div>Region 4</div>
        <div class="legend-item"><div class="color-box" style="background:#ff44ff"></div>Region 5</div>
        <div class="legend-item"><div class="color-box" style="background:#ffff44"></div>Region 6</div>
    </div>

    <script>
        const edges = {json.dumps(edges_js, indent=2)};
        const corners = {json.dumps(corners_js, indent=2)};
        
        const canvas = document.getElementById('plot');
        const ctx = canvas.getContext('2d');
        
        // Clear
        ctx.fillStyle = '#111122';
        ctx.fillRect(0, 0, canvas.width, canvas.height);
        
        // Draw edges
        ctx.strokeStyle = '#44ffaa';
        ctx.lineWidth = 1.5;
        ctx.globalAlpha = 0.7;
        
        edges.forEach(e => {{
            ctx.beginPath();
            ctx.moveTo(e.x1, e.y1);
            ctx.lineTo(e.x2, e.y2);
            ctx.stroke();
        }});
        
        ctx.globalAlpha = 1.0;
        
        // Draw corner points and labels
        corners.forEach(c => {{
            // Point
            ctx.fillStyle = c.color;
            ctx.beginPath();
            ctx.arc(c.cx, c.cy, 6, 0, Math.PI * 2);
            ctx.fill();
            
            // White border
            ctx.strokeStyle = 'white';
            ctx.lineWidth = 2;
            ctx.stroke();
            
            // Label background
            ctx.fillStyle = c.color + 'dd';
            ctx.fillRect(c.cx + 8, c.cy - 12, 28, 18);
            
            // Label text
            ctx.fillStyle = 'white';
            ctx.font = 'bold 12px Arial';
            ctx.textAlign = 'left';
            ctx.textBaseline = 'middle';
            ctx.fillText(c.id.toString(), c.cx + 12, c.cy - 3);
        }});
        
        // Build corner list
        const listDiv = document.getElementById('corner-list');
        corners.forEach(c => {{
            const div = document.createElement('div');
            div.className = 'corner-item';
            div.style.borderLeftColor = c.color;
            div.innerHTML = `
                <b>${{c.id}}</b> (R${{c.region}})<br>
                <span style="opacity:0.7">${{c.x.toFixed(2)}}, ${{c.y.toFixed(2)}}</span>
            `;
            div.addEventListener('mouseenter', () => {{
                div.style.background = 'rgba(255,255,255,0.3)';
            }});
            div.addEventListener('mouseleave', () => {{
                div.style.background = 'transparent';
            }});
            listDiv.appendChild(div);
        }});
    </script>
</body>
</html>'''
    
    with open(output_file, 'w') as f:
        f.write(html)
    
    print(f"2D Hub plot written: {output_file}")
    print(f"  Edges: {len(edges)}")
    print(f"  Corners: {len(corners)}")
    print(f"  Bounds: X=[{min_x:.3f}, {max_x:.3f}], Y=[{min_y:.3f}, {max_y:.3f}]")

def main():
    msh_file = 'T1_9/T1_9_ru_gridGmsh.msh'
    
    print("Parsing hex mesh...")
    nodes, hex_elements = parse_hex_regions(msh_file)
    print(f"Loaded {len(nodes)} nodes, {len(hex_elements)} hex elements")
    
    print("\nFinding lower boundary edges (z < 1.3)...")
    edges = find_lower_boundary_edges(nodes, hex_elements, z_max=1.3)
    
    print("Finding 2D block corners...")
    corners = find_block_corners_2d(nodes, hex_elements, z_max=1.3)
    
    print(f"\nFound {len(edges)} boundary edges near hub")
    print(f"Found {len(corners)} 2D corners")
    
    print("\nCorner list:")
    for i, corner in enumerate(corners):
        pos = corner['pos']
        print(f"  Corner {i:2d}: Node {corner['node']:6d} | Region {corner['region']} | "
              f"({pos[0]:7.3f}, {pos[1]:7.3f}, {pos[2]:7.3f})")
    
    with open('hub_2d_corners.json', 'w') as f:
        json.dump({
            'edges': len(edges),
            'corners': [{
                'id': i,
                'node': c['node'],
                'region': c['region'],
                'position': c['pos']
            } for i, c in enumerate(corners)]
        }, f, indent=2)
    
    create_2d_hub_plot(nodes, edges, corners, 'hub_2d_boundaries.html')
    
    print("\nDone! Open 'hub_2d_boundaries.html' in browser.")

if __name__ == '__main__':
    main()
