#!/usr/bin/env python3
"""
Extract inter-region boundary faces and label 40 block corners.
Create HTML wireframe showing only boundaries between regions.
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

def find_region_boundaries_and_corners(hex_elements, nodes):
    regions = defaultdict(list)
    for elem in hex_elements:
        regions[elem['geom_tag']].append(elem)
    
    # Find inter-region faces
    face_to_regions = defaultdict(set)
    
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
    
    inter_region_faces = []
    for face, regions_set in face_to_regions.items():
        if len(regions_set) >= 2:
            inter_region_faces.append({
                'nodes': list(face),
                'regions': sorted(regions_set)
            })
    
    # Find 8 BB corners per region
    all_corners = []
    
    for geom_tag in sorted(regions.keys()):
        elems = regions[geom_tag]
        region_nodes = set()
        for elem in elems:
            region_nodes.update(elem['nodes'])
        
        coords = np.array([nodes[n] for n in region_nodes])
        bb_min = coords.min(axis=0)
        bb_max = coords.max(axis=0)
        
        bb_corners = []
        for i in [0, 1]:
            for j in [0, 1]:
                for k in [0, 1]:
                    target = np.array([
                        bb_min[0] if i == 0 else bb_max[0],
                        bb_min[1] if j == 0 else bb_max[1],
                        bb_min[2] if k == 0 else bb_max[2]
                    ])
                    bb_corners.append(target)
        
        region_corners = []
        used_nodes = set()
        
        for ci, target in enumerate(bb_corners):
            best_dist = float('inf')
            best_node = None
            
            for n in region_nodes:
                if n in used_nodes:
                    continue
                dist = np.linalg.norm(nodes[n] - target)
                if dist < best_dist:
                    best_dist = dist
                    best_node = n
            
            if best_node is not None:
                used_nodes.add(best_node)
                region_corners.append({
                    'node': best_node,
                    'region': geom_tag,
                    'bb_corner_index': ci,
                    'pos': nodes[best_node].tolist()
                })
        
        all_corners.extend(region_corners)
    
    return inter_region_faces, all_corners

def create_html_wireframe(nodes, boundary_faces, corners, output_file):
    edges = set()
    for face in boundary_faces:
        n = face['nodes']
        if len(n) == 4:
            for i in range(4):
                n1, n2 = n[i], n[(i+1) % 4]
                edges.add(tuple(sorted([n1, n2])))
    
    edge_coords = []
    for n1, n2 in edges:
        p1 = nodes[n1]
        p2 = nodes[n2]
        edge_coords.append({
            'x1': float(p1[0]), 'y1': float(p1[1]), 'z1': float(p1[2]),
            'x2': float(p2[0]), 'y2': float(p2[1]), 'z2': float(p2[2])
        })
    
    region_colors = {2: '#4488ff', 3: '#44ff88', 4: '#ff8844', 5: '#ff44ff', 6: '#ffff44'}
    
    corners_js = []
    for i, corner in enumerate(corners):
        color = region_colors.get(corner['region'], '#ffffff')
        corners_js.append({
            'id': i,
            'node': corner['node'],
            'region': corner['region'],
            'x': corner['pos'][0],
            'y': corner['pos'][1],
            'z': corner['pos'][2],
            'color': color
        })
    
    all_points = np.array([nodes[n] for n in nodes])
    bb_min = all_points.min(axis=0)
    bb_max = all_points.max(axis=0)
    center = (bb_min + bb_max) / 2
    size = np.linalg.norm(bb_max - bb_min)
    
    html = f'''<!DOCTYPE html>
<html>
<head>
    <title>T1_9 Block Boundaries - {len(corners)} Corners</title>
    <style>
        body {{ margin: 0; overflow: hidden; font-family: 'Segoe UI', Arial, sans-serif; background: #0a0a1a; }}
        #info {{
            position: absolute; top: 10px; left: 10px;
            background: rgba(0,0,0,0.85); color: white;
            padding: 15px; border-radius: 8px;
            max-width: 320px; font-size: 13px;
            border: 1px solid rgba(255,255,255,0.15);
        }}
        #corner-list {{
            position: absolute; top: 10px; right: 10px;
            background: rgba(0,0,0,0.85); color: white;
            padding: 12px; border-radius: 8px;
            max-width: 250px; max-height: 90vh;
            overflow-y: auto; font-size: 11px;
            border: 1px solid rgba(255,255,255,0.15);
        }}
        .corner-item {{ 
            margin: 2px 0; padding: 3px 6px; cursor: pointer;
            border-radius: 4px; border-left: 3px solid;
            transition: all 0.2s;
        }}
        .corner-item:hover {{ background: rgba(255,255,255,0.15); }}
        .corner-item.highlight {{ background: rgba(255,255,0,0.25) !important; }}
        h3 {{ margin: 0 0 10px 0; font-size: 16px; color: #44ff88; }}
        h4 {{ margin: 0 0 8px 0; font-size: 12px; color: #ffaa44; }}
        .stat {{ margin: 3px 0; }}
        .legend {{ display: flex; flex-wrap: wrap; gap: 6px; margin-top: 8px; }}
        .legend-item {{ display: flex; align-items: center; gap: 3px; font-size: 10px; }}
        .color-box {{ width: 10px; height: 10px; border-radius: 2px; }}
    </style>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
</head>
<body>
    <div id="info">
        <h3>Block Boundary Wireframe</h3>
        <div class="stat"><b>Inter-Region Faces:</b> {len(boundary_faces)}</div>
        <div class="stat"><b>Boundary Edges:</b> {len(edges)}</div>
        <div class="stat"><b>Block Corners:</b> {len(corners)}</div>
        <div class="stat" style="margin-top:8px; font-size:11px; opacity:0.7;">
            Shows only boundaries between regions<br>
            (not outer mesh boundaries)
        </div>
        <div class="legend">
            <div class="legend-item"><div class="color-box" style="background:#4488ff"></div>R2</div>
            <div class="legend-item"><div class="color-box" style="background:#44ff88"></div>R3</div>
            <div class="legend-item"><div class="color-box" style="background:#ff8844"></div>R4</div>
            <div class="legend-item"><div class="color-box" style="background:#ff44ff"></div>R5</div>
            <div class="legend-item"><div class="color-box" style="background:#ffff44"></div>R6</div>
        </div>
        <div style="margin-top:10px; font-size:11px; opacity:0.6;">
            Left: Rotate | Right: Pan | Scroll: Zoom
        </div>
    </div>
    
    <div id="corner-list">
        <h4>Corners ({len(corners)} total)</h4>
        <div id="corner-items"></div>
    </div>

    <script>
        const corners = {json.dumps(corners_js, indent=2)};
        const edges = {json.dumps(edge_coords, indent=2)};
        
        const scene = new THREE.Scene();
        scene.background = new THREE.Color(0x0a0a1a);
        
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
        
        const ambientLight = new THREE.AmbientLight(0x6060a0, 0.4);
        scene.add(ambientLight);
        const dirLight = new THREE.DirectionalLight(0xffffff, 0.7);
        dirLight.position.set(2, 3, 5);
        scene.add(dirLight);
        
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
            opacity: 0.5, 
            transparent: true
        }});
        const lineMesh = new THREE.LineSegments(lineGeometry, lineMaterial);
        scene.add(lineMesh);
        
        const cornerGroup = new THREE.Group();
        const sphereMeshes = [];
        
        corners.forEach((corner) => {{
            const radius = 0.018 * {size};
            
            const geometry = new THREE.SphereGeometry(radius, 32, 32);
            const material = new THREE.MeshPhongMaterial({{
                color: corner.color,
                emissive: corner.color,
                emissiveIntensity: 0.5,
                shininess: 100
            }});
            const sphere = new THREE.Mesh(geometry, material);
            sphere.position.set(corner.x, corner.y, corner.z);
            sphere.userData = {{ cornerId: corner.id }};
            cornerGroup.add(sphere);
            sphereMeshes.push(sphere);
            
            const canvas = document.createElement('canvas');
            const ctx = canvas.getContext('2d');
            canvas.width = 96;
            canvas.height = 48;
            
            ctx.fillStyle = corner.color + 'ee';
            ctx.beginPath();
            ctx.roundRect(2, 2, 92, 44, 8);
            ctx.fill();
            
            ctx.strokeStyle = 'white';
            ctx.lineWidth = 2;
            ctx.stroke();
            
            ctx.fillStyle = 'white';
            ctx.font = 'bold 24px Arial';
            ctx.textAlign = 'center';
            ctx.textBaseline = 'middle';
            ctx.fillText(corner.id.toString(), 48, 24);
            
            const texture = new THREE.CanvasTexture(canvas);
            const spriteMat = new THREE.SpriteMaterial({{ 
                map: texture, 
                depthTest: false,
                transparent: true
            }});
            const sprite = new THREE.Sprite(spriteMat);
            sprite.position.set(corner.x + radius*2.5, corner.y + radius*2.5, corner.z + radius*2.5);
            sprite.scale.set(radius*6, radius*3, 1);
            cornerGroup.add(sprite);
        }});
        
        scene.add(cornerGroup);
        
        const listDiv = document.getElementById('corner-items');
        corners.forEach((corner) => {{
            const div = document.createElement('div');
            div.className = 'corner-item';
            div.style.borderLeftColor = corner.color;
            div.id = 'corner-' + corner.id;
            
            div.innerHTML = `
                <b>Corner ${{corner.id}}</b> (R${{corner.region}})<br>
                Node ${{corner.node}}<br>
                <span style="opacity:0.7; font-size:10px;">
                    (${{corner.x.toFixed(2)}}, ${{corner.y.toFixed(2)}}, ${{corner.z.toFixed(2)}})
                </span>
            `;
            
            div.addEventListener('mouseenter', () => {{
                div.classList.add('highlight');
                const s = sphereMeshes.find(m => m.userData.cornerId === corner.id);
                if (s) {{
                    s.material.emissiveIntensity = 1.0;
                    s.scale.setScalar(2.0);
                }}
            }});
            
            div.addEventListener('mouseleave', () => {{
                div.classList.remove('highlight');
                const s = sphereMeshes.find(m => m.userData.cornerId === corner.id);
                if (s) {{
                    s.material.emissiveIntensity = 0.5;
                    s.scale.setScalar(1.0);
                }}
            }});
            
            listDiv.appendChild(div);
        }});
        
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
    
    print(f"HTML wireframe written: {output_file}")
    print(f"  Inter-region faces: {len(boundary_faces)}")
    print(f"  Boundary edges: {len(edges)}")
    print(f"  Block corners: {len(corners)}")

def main():
    msh_file = 'T1_9/T1_9_ru_gridGmsh.msh'
    
    print("Parsing hex mesh...")
    nodes, hex_elements = parse_hex_regions(msh_file)
    print(f"Loaded {len(nodes)} nodes, {len(hex_elements)} hex elements")
    
    print("\nFinding inter-region boundaries and corners...")
    boundary_faces, corners = find_region_boundaries_and_corners(hex_elements, nodes)
    
    print(f"\nFound {len(boundary_faces)} inter-region boundary faces")
    print(f"Found {len(corners)} block corners")
    
    print("\nCorner list:")
    for i, corner in enumerate(corners):
        pos = corner['pos']
        print(f"  Corner {i:2d}: Node {corner['node']:6d} | Region {corner['region']} | "
              f"({pos[0]:7.3f}, {pos[1]:7.3f}, {pos[2]:7.3f})")
    
    with open('block_corners_final.json', 'w') as f:
        json.dump({
            'boundary_faces': len(boundary_faces),
            'corners': [{
                'id': i,
                'node': c['node'],
                'region': c['region'],
                'position': c['pos']
            } for i, c in enumerate(corners)]
        }, f, indent=2)
    
    create_html_wireframe(nodes, boundary_faces, corners, 'boundary_wireframe_final.html')
    
    print("\nDone! Open 'boundary_wireframe_final.html' in browser.")

if __name__ == '__main__':
    main()
