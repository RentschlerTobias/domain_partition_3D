#!/usr/bin/env python3
"""
Extract actual 8 corners from each hex region and create labeled wireframe.
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

def find_region_corners(nodes, hex_elements):
    """Find 8 corners for each hex region using boundary face analysis."""
    regions = defaultdict(list)
    for elem in hex_elements:
        regions[elem['geom_tag']].append(elem)
    
    all_corners = []
    
    for geom_tag in sorted(regions.keys()):
        elems = regions[geom_tag]
        
        # Build face connectivity within region
        face_to_elems = defaultdict(list)
        for elem in elems:
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
                face_to_elems[face].append(elem['tag'])
        
        # Count boundary faces for each node
        node_boundary_faces = defaultdict(int)
        for face, face_elems in face_to_elems.items():
            if len(face_elems) == 1:  # boundary face
                for node in face:
                    node_boundary_faces[node] += 1
        
        # Corner nodes have 3 boundary faces (in a structured hex block)
        # But due to curved geometry, some may have more
        corner_candidates = [(n, c) for n, c in node_boundary_faces.items() if c >= 3]
        corner_candidates.sort(key=lambda x: (-x[1], x[0]))
        
        # Take top 8 as corners
        region_corners = []
        for i, (node, count) in enumerate(corner_candidates[:8]):
            pos = nodes[node]
            region_corners.append({
                'global_id': len(all_corners) + i,
                'node_id': node,
                'region': geom_tag,
                'position': pos.tolist(),
                'boundary_faces': count,
                'corner_index': i
            })
        
        all_corners.extend(region_corners)
        print(f"Region {geom_tag}: {len(region_corners)} corners found")
    
    return all_corners

def create_labeled_wireframe(nodes, hex_elements, corners, output_file):
    """Create HTML wireframe with numbered corner labels."""
    
    # Build unique edges from all hex elements
    edges = set()
    for elem in hex_elements:
        n = elem['nodes']
        edge_list = [
            (n[0], n[1]), (n[1], n[2]), (n[2], n[3]), (n[3], n[0]),
            (n[4], n[5]), (n[5], n[6]), (n[6], n[7]), (n[7], n[4]),
            (n[0], n[4]), (n[1], n[5]), (n[2], n[6]), (n[3], n[7])
        ]
        for e in edge_list:
            edges.add(tuple(sorted(e)))
    
    # Sample edges for performance
    edge_list = list(edges)
    if len(edge_list) > 15000:
        import random
        random.seed(42)
        edge_list = random.sample(edge_list, 15000)
    
    edge_coords = []
    for n1, n2 in edge_list:
        p1 = nodes[n1]
        p2 = nodes[n2]
        edge_coords.append({
            'x1': float(p1[0]), 'y1': float(p1[1]), 'z1': float(p1[2]),
            'x2': float(p2[0]), 'y2': float(p2[1]), 'z2': float(p2[2])
        })
    
    # Corner data
    corners_js = []
    region_colors = {2: '#4488ff', 3: '#44ff88', 4: '#ff8844', 5: '#ff44ff', 6: '#ffff44'}
    
    for corner in corners:
        color = region_colors.get(corner['region'], '#ffffff')
        corners_js.append({
            'id': corner['global_id'],
            'node_id': corner['node_id'],
            'region': corner['region'],
            'x': corner['position'][0],
            'y': corner['position'][1],
            'z': corner['position'][2],
            'color': color,
            'boundary_faces': corner['boundary_faces']
        })
    
    # Bounding box
    all_points = np.array(list(nodes.values()))
    bb_min = all_points.min(axis=0)
    bb_max = all_points.max(axis=0)
    center = (bb_min + bb_max) / 2
    size = np.linalg.norm(bb_max - bb_min)
    
    html = f'''<!DOCTYPE html>
<html>
<head>
    <title>T1_9 Hex Block Corners - {len(corners)} Points</title>
    <style>
        body {{ margin: 0; overflow: hidden; font-family: 'Segoe UI', Arial, sans-serif; background: #0a0a1a; }}
        #info {{
            position: absolute; top: 10px; left: 10px;
            background: rgba(0,0,0,0.85); color: white;
            padding: 15px; border-radius: 8px;
            max-width: 350px; font-size: 13px;
            border: 1px solid rgba(255,255,255,0.15);
            box-shadow: 0 4px 15px rgba(0,0,0,0.5);
        }}
        #corner-list {{
            position: absolute; top: 10px; right: 10px;
            background: rgba(0,0,0,0.85); color: white;
            padding: 12px; border-radius: 8px;
            max-width: 250px; max-height: 90vh;
            overflow-y: auto; font-size: 11px;
            border: 1px solid rgba(255,255,255,0.15);
            box-shadow: 0 4px 15px rgba(0,0,0,0.5);
        }}
        .corner-item {{ 
            margin: 2px 0; padding: 3px 6px; cursor: pointer;
            border-radius: 4px; border-left: 3px solid;
            transition: all 0.2s;
        }}
        .corner-item:hover {{ 
            background: rgba(255,255,255,0.15) !important;
            transform: translateX(-2px);
        }}
        .corner-item.highlight {{ 
            background: rgba(255,255,255,0.25) !important;
        }}
        h3 {{ margin: 0 0 12px 0; font-size: 16px; color: #44ff88; }}
        h4 {{ margin: 0 0 8px 0; font-size: 13px; color: #ffaa44; }}
        .stat {{ margin: 4px 0; opacity: 0.9; }}
        .region-legend {{ display: flex; flex-wrap: wrap; gap: 8px; margin-top: 10px; }}
        .legend-item {{ display: flex; align-items: center; gap: 4px; font-size: 11px; }}
        .color-box {{ width: 12px; height: 12px; border-radius: 2px; }}
    </style>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
</head>
<body>
    <div id="info">
        <h3>T1_9 Hex Block Corners</h3>
        <div class="stat"><b>Total Corners:</b> {len(corners)}</div>
        <div class="stat"><b>Wireframe Edges:</b> {len(edge_coords)}</div>
        <div class="stat"><b>Hex Elements:</b> {len(hex_elements)}</div>
        <div class="stat" style="margin-top:10px; font-size:11px; opacity:0.7;">
            Red/Green/Orange/Pink/Yellow = Region 2/3/4/5/6
        </div>
        <div class="region-legend">
            <div class="legend-item"><div class="color-box" style="background:#4488ff"></div>Region 2</div>
            <div class="legend-item"><div class="color-box" style="background:#44ff88"></div>Region 3</div>
            <div class="legend-item"><div class="color-box" style="background:#ff8844"></div>Region 4</div>
            <div class="legend-item"><div class="color-box" style="background:#ff44ff"></div>Region 5</div>
            <div class="legend-item"><div class="color-box" style="background:#ffff44"></div>Region 6</div>
        </div>
        <div style="margin-top:12px; font-size:11px; opacity:0.6;">
            Left: Rotate | Right: Pan | Scroll: Zoom
        </div>
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
        scene.fog = new THREE.Fog(0x0a0a1a, {size*0.5}, {size*2});
        
        const camera = new THREE.PerspectiveCamera(60, window.innerWidth / window.innerHeight, 0.01, 100);
        camera.position.set({center[0] + size*0.5}, {center[1] + size*0.3}, {center[2] + size*0.8});
        
        const renderer = new THREE.WebGLRenderer({{ antialias: true }});
        renderer.setSize(window.innerWidth, window.innerHeight);
        renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
        document.body.appendChild(renderer.domElement);
        
        const controls = new THREE.OrbitControls(camera, renderer.domElement);
        controls.target.set({center[0]}, {center[1]}, {center[2]});
        controls.enableDamping = true;
        controls.dampingFactor = 0.05;
        controls.update();
        
        // Lighting
        const ambientLight = new THREE.AmbientLight(0x6060a0, 0.4);
        scene.add(ambientLight);
        const dirLight = new THREE.DirectionalLight(0xffffff, 0.7);
        dirLight.position.set(2, 3, 5);
        scene.add(dirLight);
        const dirLight2 = new THREE.DirectionalLight(0x4444ff, 0.3);
        dirLight2.position.set(-2, -1, -3);
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
            color: 0x44ffaa, 
            opacity: 0.15, 
            transparent: true 
        }});
        const lineMesh = new THREE.LineSegments(lineGeometry, lineMaterial);
        scene.add(lineMesh);
        
        // Corner spheres and labels
        const cornerGroup = new THREE.Group();
        const sphereMeshes = [];
        
        corners.forEach((corner) => {{
            const radius = 0.015 * {size};
            
            // Sphere
            const geometry = new THREE.SphereGeometry(radius, 24, 24);
            const material = new THREE.MeshPhongMaterial({{
                color: corner.color,
                emissive: corner.color,
                emissiveIntensity: 0.3,
                shininess: 80
            }});
            const sphere = new THREE.Mesh(geometry, material);
            sphere.position.set(corner.x, corner.y, corner.z);
            sphere.userData = {{ cornerId: corner.id }};
            cornerGroup.add(sphere);
            sphereMeshes.push(sphere);
            
            // Label
            const canvas = document.createElement('canvas');
            const ctx = canvas.getContext('2d');
            canvas.width = 128;
            canvas.height = 64;
            
            ctx.fillStyle = corner.color + 'dd';
            ctx.beginPath();
            ctx.roundRect(2, 2, 124, 60, 10);
            ctx.fill();
            
            ctx.strokeStyle = 'white';
            ctx.lineWidth = 2;
            ctx.stroke();
            
            ctx.fillStyle = 'white';
            ctx.font = 'bold 28px Arial';
            ctx.textAlign = 'center';
            ctx.textBaseline = 'middle';
            ctx.fillText(corner.id.toString(), 64, 32);
            
            const texture = new THREE.CanvasTexture(canvas);
            const spriteMat = new THREE.SpriteMaterial({{ 
                map: texture, 
                depthTest: false,
                transparent: true
            }});
            const sprite = new THREE.Sprite(spriteMat);
            sprite.position.set(corner.x + radius*2, corner.y + radius*2, corner.z + radius*2);
            sprite.scale.set(radius*5, radius*2.5, 1);
            sprite.userData = {{ cornerId: corner.id }};
            cornerGroup.add(sprite);
        }});
        
        scene.add(cornerGroup);
        
        // Build corner list
        const listDiv = document.getElementById('corner-items');
        corners.forEach((corner) => {{
            const div = document.createElement('div');
            div.className = 'corner-item';
            div.style.borderLeftColor = corner.color;
            div.id = 'corner-' + corner.id;
            
            div.innerHTML = `
                <b>Corner ${{corner.id}}</b> (R${{corner.region}})<br>
                Node ${{corner.node_id}}<br>
                <span style="opacity:0.7; font-size:10px;">
                    (${{corner.x.toFixed(2)}}, ${{corner.y.toFixed(2)}}, ${{corner.z.toFixed(2)}})
                </span>
            `;
            
            div.addEventListener('mouseenter', () => {{
                div.classList.add('highlight');
                const s = sphereMeshes.find(m => m.userData.cornerId === corner.id);
                if (s) {{
                    s.material.emissiveIntensity = 0.9;
                    s.scale.setScalar(1.8);
                }}
            }});
            
            div.addEventListener('mouseleave', () => {{
                div.classList.remove('highlight');
                const s = sphereMeshes.find(m => m.userData.cornerId === corner.id);
                if (s) {{
                    s.material.emissiveIntensity = 0.3;
                    s.scale.setScalar(1.0);
                }}
            }});
            
            listDiv.appendChild(div);
        }});
        
        // Animation
        function animate() {{
            requestAnimationFrame(animate);
            controls.update();
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
        f.write(html)
    
    print(f"\nLabeled wireframe written: {output_file}")
    print(f"  Corners: {len(corners)}")
    print(f"  Wireframe edges: {len(edge_coords)}")

def main():
    msh_file = 'T1_9/T1_9_ru_gridGmsh.msh'
    
    print("Parsing hex mesh regions...")
    nodes, hex_elements = parse_hex_regions(msh_file)
    print(f"Loaded {len(nodes)} nodes, {len(hex_elements)} hex elements")
    
    print("\nFinding block corners...")
    corners = find_region_corners(nodes, hex_elements)
    
    print(f"\nTotal corners found: {len(corners)}")
    print("\nCorner list:")
    for corner in corners:
        pos = corner['position']
        print(f"  Corner {corner['global_id']:2d}: Node {corner['node_id']:6d} | "
              f"Region {corner['region']} | "
              f"({pos[0]:7.3f}, {pos[1]:7.3f}, {pos[2]:7.3f})")
    
    # Save to JSON
    with open('block_corners_40.json', 'w') as f:
        json.dump(corners, f, indent=2)
    print("\nSaved to block_corners_40.json")
    
    # Create HTML
    create_labeled_wireframe(nodes, hex_elements, corners, 'block_corners_labeled.html')
    
    print("\nDone! Open 'block_corners_labeled.html' in browser.")

if __name__ == '__main__':
    main()
