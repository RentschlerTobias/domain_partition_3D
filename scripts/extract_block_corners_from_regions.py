#!/usr/bin/env python3
"""
Extract block corners from hex mesh regions.
A corner is a node that belongs to elements in 3 or more different regions.
"""

import numpy as np
import json
from collections import defaultdict

def parse_hex_regions(msh_file):
    """Parse hex elements and their region assignments from .msh file."""
    nodes = {}
    hex_elements = []
    
    with open(msh_file, 'r') as f:
        lines = f.readlines()
    
    # Parse nodes
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
    
    # Parse hex elements
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
            if elem_type == 5:  # Hex
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

def find_region_corners(nodes, hex_elements):
    """
    Find corners where 3 or more regions meet.
    Also find edges where 2 regions meet.
    """
    # For each node, find which regions it belongs to
    node_regions = defaultdict(set)
    
    for elem in hex_elements:
        for node in elem['nodes']:
            node_regions[node].add(elem['geom_tag'])
    
    # Nodes in 3+ regions are corners
    corners = []
    for node, regions in node_regions.items():
        if len(regions) >= 3:
            corners.append({
                'node_id': node,
                'position': nodes[node].tolist(),
                'num_regions': len(regions),
                'regions': sorted(regions)
            })
    
    # Sort by number of regions (descending), then by position
    corners.sort(key=lambda x: (-x['num_regions'], x['position'][2], x['position'][0]))
    
    return corners

def create_wireframe_html(nodes, hex_elements, corners, output_file='block_corners_wireframe.html'):
    """Create HTML wireframe with all hex edges and numbered corner labels."""
    
    # Build unique edges from all hex elements
    edges = set()
    for elem in hex_elements:
        n = elem['nodes']
        # 12 edges of hexahedron
        edge_list = [
            (n[0], n[1]), (n[1], n[2]), (n[2], n[3]), (n[3], n[0]),  # bottom face
            (n[4], n[5]), (n[5], n[6]), (n[6], n[7]), (n[7], n[4]),  # top face
            (n[0], n[4]), (n[1], n[5]), (n[2], n[6]), (n[3], n[7])   # vertical edges
        ]
        for e in edge_list:
            edges.add(tuple(sorted(e)))
    
    # Convert to coordinate pairs
    edge_coords = []
    for n1, n2 in edges:
        p1 = nodes[n1]
        p2 = nodes[n2]
        edge_coords.append({
            'x1': float(p1[0]), 'y1': float(p1[1]), 'z1': float(p1[2]),
            'x2': float(p2[0]), 'y2': float(p2[1]), 'z2': float(p2[2])
        })
    
    # Limit edges for performance (sample every Nth edge)
    if len(edge_coords) > 10000:
        step = len(edge_coords) // 10000 + 1
        edge_coords = edge_coords[::step]
    
    # Corner data
    corners_js = []
    for i, corner in enumerate(corners):
        corners_js.append({
            'id': i,
            'node_id': corner['node_id'],
            'x': corner['position'][0],
            'y': corner['position'][1],
            'z': corner['position'][2],
            'num_regions': corner['num_regions'],
            'regions': corner['regions']
        })
    
    # Bounding box
    all_points = np.array(list(nodes.values()))
    bb_min = all_points.min(axis=0)
    bb_max = all_points.max(axis=0)
    center = (bb_min + bb_max) / 2
    size = np.linalg.norm(bb_max - bb_min)
    
    html_content = f'''<!DOCTYPE html>
<html>
<head>
    <title>T1_9 Hex Block Corners</title>
    <style>
        body {{ margin: 0; overflow: hidden; font-family: 'Segoe UI', Arial, sans-serif; }}
        #info {{
            position: absolute; top: 10px; left: 10px;
            background: rgba(0,0,0,0.8); color: white;
            padding: 15px; border-radius: 8px;
            max-width: 400px; font-size: 14px;
            border: 1px solid rgba(255,255,255,0.2);
        }}
        #corner-list {{
            position: absolute; top: 10px; right: 10px;
            background: rgba(0,0,0,0.8); color: white;
            padding: 15px; border-radius: 8px;
            max-width: 300px; max-height: 90vh;
            overflow-y: auto; font-size: 11px;
            border: 1px solid rgba(255,255,255,0.2);
        }}
        .corner-item {{ 
            margin: 3px 0; padding: 4px 8px; cursor: pointer;
            border-radius: 4px; border-left: 3px solid #ffaa00;
            transition: background 0.2s;
        }}
        .corner-item:hover {{ background: rgba(255,255,255,0.15); }}
        .corner-item.highlight {{ 
            background: rgba(255,255,0,0.25) !important;
            border-left-color: #ffff00;
        }}
        h3, h4 {{ margin: 0 0 10px 0; }}
        .region-badge {{
            display: inline-block; padding: 1px 5px;
            border-radius: 3px; font-size: 10px; margin: 0 2px;
            background: rgba(255,255,255,0.2);
        }}
    </style>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
</head>
<body>
    <div id="info">
        <h3>T1_9 Hex Block Corners</h3>
        <p><b>Corners Found:</b> {len(corners)}</p>
        <p><b>Wireframe Edges:</b> {len(edge_coords)}</p>
        <p>Red spheres = 4-region corners</p>
        <p>Orange spheres = 3-region corners</p>
        <p style="font-size:12px; opacity:0.8; margin-top:10px;">
            Left: Rotate | Right: Pan | Scroll: Zoom
        </p>
    </div>
    
    <div id="corner-list">
        <h4>Corners ({len(corners)} total)</h4>
        <div id="corner-items"></div>
    </div>

    <script>
        // Data
        const corners = {json.dumps(corners_js, indent=2)};
        const edges = {json.dumps(edge_coords, indent=2)};
        
        // Scene
        const scene = new THREE.Scene();
        scene.background = new THREE.Color(0x0a0a1a);
        
        const camera = new THREE.PerspectiveCamera(60, window.innerWidth / window.innerHeight, 0.01, 100);
        camera.position.set({center[0] + size*0.6}, {center[1] + size*0.4}, {center[2] + size*0.9});
        
        const renderer = new THREE.WebGLRenderer({{ antialias: true }});
        renderer.setSize(window.innerWidth, window.innerHeight);
        renderer.setPixelRatio(window.devicePixelRatio);
        document.body.appendChild(renderer.domElement);
        
        const controls = new THREE.OrbitControls(camera, renderer.domElement);
        controls.target.set({center[0]}, {center[1]}, {center[2]});
        controls.enableDamping = true;
        controls.dampingFactor = 0.05;
        controls.update();
        
        // Lighting
        const ambientLight = new THREE.AmbientLight(0x404080, 0.5);
        scene.add(ambientLight);
        const dirLight = new THREE.DirectionalLight(0xffffff, 0.8);
        dirLight.position.set(1, 1, 2);
        scene.add(dirLight);
        const dirLight2 = new THREE.DirectionalLight(0xffffff, 0.3);
        dirLight2.position.set(-1, -1, -1);
        scene.add(dirLight2);
        
        // Wireframe
        const lineGeometry = new THREE.BufferGeometry();
        const positions = new Float32Array(edges.length * 6);
        for (let i = 0; i < edges.length; i++) {{
            const e = edges[i];
            positions[i*6] = e.x1;   positions[i*6+1] = e.y1;   positions[i*6+2] = e.z1;
            positions[i*6+3] = e.x2; positions[i*6+4] = e.y2; positions[i*6+5] = e.z2;
        }}
        lineGeometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
        const lineMaterial = new THREE.LineBasicMaterial({{ 
            color: 0x44ff88, 
            opacity: 0.25, 
            transparent: true 
        }});
        const lineMesh = new THREE.LineSegments(lineGeometry, lineMaterial);
        scene.add(lineMesh);
        
        // Corner spheres and labels
        const cornerGroup = new THREE.Group();
        const sphereMeshes = [];
        const labelSprites = [];
        
        corners.forEach((corner, idx) => {{
            const isMulti = corner.num_regions >= 4;
            const radius = isMulti ? 0.012 * {size} : 0.008 * {size};
            
            // Sphere
            const geometry = new THREE.SphereGeometry(radius, 32, 32);
            const material = new THREE.MeshPhongMaterial({{
                color: isMulti ? 0xff3333 : 0xff8800,
                emissive: isMulti ? 0xff0000 : 0xff4400,
                emissiveIntensity: 0.4,
                shininess: 100
            }});
            const sphere = new THREE.Mesh(geometry, material);
            sphere.position.set(corner.x, corner.y, corner.z);
            sphere.userData = {{ cornerId: idx, originalScale: radius }};
            cornerGroup.add(sphere);
            sphereMeshes.push(sphere);
            
            // Label canvas
            const canvas = document.createElement('canvas');
            const ctx = canvas.getContext('2d');
            canvas.width = 96;
            canvas.height = 48;
            
            // Background
            ctx.fillStyle = isMulti ? 'rgba(255,50,50,0.85)' : 'rgba(255,150,0,0.85)';
            ctx.beginPath();
            ctx.roundRect(0, 0, 96, 48, 8);
            ctx.fill();
            
            // Border
            ctx.strokeStyle = 'white';
            ctx.lineWidth = 2;
            ctx.stroke();
            
            // Text
            ctx.fillStyle = 'white';
            ctx.font = 'bold 22px Arial';
            ctx.textAlign = 'center';
            ctx.textBaseline = 'middle';
            ctx.fillText(idx.toString(), 48, 24);
            
            const texture = new THREE.CanvasTexture(canvas);
            const spriteMat = new THREE.SpriteMaterial({{ map: texture, depthTest: false }});
            const sprite = new THREE.Sprite(spriteMat);
            const labelOffset = radius * 2.5;
            sprite.position.set(
                corner.x + labelOffset, 
                corner.y + labelOffset, 
                corner.z + labelOffset
            );
            sprite.scale.set(radius * 4, radius * 2, 1);
            sprite.userData = {{ cornerId: idx }};
            cornerGroup.add(sprite);
            labelSprites.push(sprite);
        }});
        
        scene.add(cornerGroup);
        
        // Build corner list
        const listDiv = document.getElementById('corner-items');
        corners.forEach((corner, idx) => {{
            const div = document.createElement('div');
            div.className = 'corner-item';
            div.id = 'corner-item-' + idx;
            
            const regionsHtml = corner.regions.map(r => 
                '<span class="region-badge">R' + r + '</span>'
            ).join('');
            
            div.innerHTML = `
                <b>Corner ${{idx}}</b> (Node ${{corner.node_id}})<br>
                Regions: ${{regionsHtml}}<br>
                <span style="opacity:0.7">
                    (${{corner.x.toFixed(2)}}, ${{corner.y.toFixed(2)}}, ${{corner.z.toFixed(2)}})
                </span>
            `;
            
            div.addEventListener('mouseenter', () => {{
                div.classList.add('highlight');
                const s = sphereMeshes[idx];
                if (s) {{
                    s.material.emissiveIntensity = 1.0;
                    s.scale.setScalar(1.5);
                }}
            }});
            
            div.addEventListener('mouseleave', () => {{
                div.classList.remove('highlight');
                const s = sphereMeshes[idx];
                if (s) {{
                    s.material.emissiveIntensity = 0.4;
                    s.scale.setScalar(1.0);
                }}
            }});
            
            listDiv.appendChild(div);
        }});
        
        // Animation
        function animate() {{
            requestAnimationFrame(animate);
            controls.update();
            
            // Make labels always face camera
            labelSprites.forEach(sprite => {{
                sprite.lookAt(camera.position);
            }});
            
            renderer.render(scene, camera);
        }}
        
        window.addEventListener('resize', () => {{
            camera.aspect = window.innerWidth / window.innerHeight;
            camera.updateProjectionMatrix();
            renderer.setSize(window.innerWidth, window.innerHeight);
        }});
        
        animate();
    </script>
</body>
</html>'''
    
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    print(f"\nHTML wireframe written: {output_file}")
    print(f"  Corners: {len(corners)}")
    print(f"  Wireframe edges: {len(edge_coords)}")
    print(f"  Scene size: {size:.3f}")

def main():
    msh_file = 'T1_9/T1_9_ru_gridGmsh.msh'
    
    print("Parsing hex mesh regions...")
    nodes, hex_elements = parse_hex_regions(msh_file)
    print(f"Loaded {len(nodes)} nodes, {len(hex_elements)} hex elements")
    
    print("\nFinding block corners (nodes in 3+ regions)...")
    corners = find_region_corners(nodes, hex_elements)
    
    print(f"\nFound {len(corners)} corners:")
    for i, corner in enumerate(corners):
        print(f"  Corner {i}: Node {corner['node_id']}, "
              f"{corner['num_regions']} regions {corner['regions']}, "
              f"pos=({corner['position'][0]:.3f}, {corner['position'][1]:.3f}, {corner['position'][2]:.3f})")
    
    # Save to JSON
    with open('block_corners.json', 'w') as f:
        json.dump(corners, f, indent=2)
    print("\nCorners saved to block_corners.json")
    
    # Create HTML wireframe
    create_wireframe_html(nodes, hex_elements, corners, 'block_corners_wireframe.html')
    
    print("\nDone! Open 'block_corners_wireframe.html' in a browser to view.")

if __name__ == '__main__':
    main()
