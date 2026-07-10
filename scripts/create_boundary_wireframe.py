#!/usr/bin/env python3

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

def find_boundary_faces_and_corners(hex_elements):
    """Find boundary faces and corner nodes."""
    
    face_to_elements = defaultdict(list)
    
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
            face_to_elements[face].append(elem)
    
    boundary_faces = []
    for face, elems in face_to_elements.items():
        if len(elems) == 1:
            boundary_faces.append({
                'nodes': list(face),
                'element': elems[0],
                'region': elems[0]['geom_tag']
            })
    
    node_boundary_count = defaultdict(int)
    node_regions = defaultdict(set)
    
    for face in boundary_faces:
        for node in face['nodes']:
            node_boundary_count[node] += 1
            node_regions[node].add(face['region'])
    
    corners = []
    for node, count in node_boundary_count.items():
        if count >= 3:
            corners.append({
                'node_id': node,
                'boundary_faces': count,
                'regions': sorted(node_regions[node]),
                'num_regions': len(node_regions[node])
            })
    
    corners.sort(key=lambda x: (-x['boundary_faces'], x['node_id']))
    
    return boundary_faces, corners

def create_boundary_wireframe(nodes, boundary_faces, corners, output_file):
    
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
    
    corners_js = []
    for i, corner in enumerate(corners):
        pos = nodes[corner['node_id']]
        is_multi = corner['num_regions'] >= 4
        corners_js.append({
            'id': i,
            'node_id': corner['node_id'],
            'num_regions': corner['num_regions'],
            'regions': corner['regions'],
            'x': float(pos[0]),
            'y': float(pos[1]),
            'z': float(pos[2]),
            'color': '#ff3333' if is_multi else '#ff8800'
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
            max-width: 260px; max-height: 90vh;
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
    </style>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
</head>
<body>
    <div id="info">
        <h3>Block Boundary Wireframe</h3>
        <div class="stat"><b>Boundary Faces:</b> {len(boundary_faces)}</div>
        <div class="stat"><b>Boundary Edges:</b> {len(edges)}</div>
        <div class="stat"><b>Corner Nodes:</b> {len(corners)}</div>
        <div class="stat" style="margin-top:8px; font-size:11px;">
            Red = 4+ region corners<br>
            Orange = 3 region corners<br>
            Green lines = region boundaries
        </div>
        <div style="margin-top:10px; font-size:11px; opacity:0.7;">
            Left: Rotate | Right: Pan | Scroll: Zoom
        </div>
    </div>
    
    <div id="corner-list">
        <h4>Corner Nodes ({len(corners)} total)</h4>
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
        
        // Boundary edges - thicker and brighter
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
            opacity: 0.6, 
            transparent: true,
            linewidth: 2
        }});
        const lineMesh = new THREE.LineSegments(lineGeometry, lineMaterial);
        scene.add(lineMesh);
        
        // Corner spheres and labels
        const cornerGroup = new THREE.Group();
        const sphereMeshes = [];
        
        corners.forEach((corner) => {{
            const radius = 0.02 * {size};
            
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
            
            // Label
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
            ctx.font = 'bold 26px Arial';
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
        
        // Build corner list
        const listDiv = document.getElementById('corner-items');
        corners.forEach((corner) => {{
            const div = document.createElement('div');
            div.className = 'corner-item';
            div.style.borderLeftColor = corner.color;
            div.id = 'corner-' + corner.id;
            
            div.innerHTML = `
                <b>Corner ${{corner.id}}</b> (Node ${{corner.node_id}})<br>
                Regions: ${{corner.regions.join(', ')}}<br>
                <span style="opacity:0.7">
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
    
    print(f"Boundary wireframe written: {output_file}")
    print(f"  Boundary faces: {len(boundary_faces)}")
    print(f"  Boundary edges: {len(edges)}")
    print(f"  Corner nodes: {len(corners)}")

def main():
    msh_file = 'T1_9/T1_9_ru_gridGmsh.msh'
    
    print("Parsing hex mesh...")
    nodes, hex_elements = parse_hex_regions(msh_file)
    print(f"Loaded {len(nodes)} nodes, {len(hex_elements)} hex elements")
    
    print("\nFinding boundary faces and corners...")
    boundary_faces, corners = find_boundary_faces_and_corners(hex_elements)
    
    print(f"\nFound {len(boundary_faces)} boundary faces")
    print(f"Found {len(corners)} corner nodes")
    
    print("\nCorner nodes:")
    for i, corner in enumerate(corners):
        pos = nodes[corner['node_id']]
        print(f"  Corner {i:2d}: Node {corner['node_id']:6d} | "
              f"{corner['num_regions']} regions {corner['regions']} | "
              f"({pos[0]:7.3f}, {pos[1]:7.3f}, {pos[2]:7.3f})")
    
    with open('boundary_corners.json', 'w') as f:
        json.dump({
            'boundary_faces': len(boundary_faces),
            'corners': [{
                'id': i,
                'node_id': c['node_id'],
                'num_regions': c['num_regions'],
                'regions': c['regions'],
                'position': nodes[c['node_id']].tolist()
            } for i, c in enumerate(corners)]
        }, f, indent=2)
    
    create_boundary_wireframe(nodes, boundary_faces, corners, 'boundary_wireframe_labeled.html')
    
    print("\nDone! Open 'boundary_wireframe_labeled.html' in browser.")

if __name__ == '__main__':
    main()
