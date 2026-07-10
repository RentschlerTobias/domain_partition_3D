#!/usr/bin/env python3
"""
Create proper turbine passage visualization showing:
1. Blade surface (red, semi-transparent)
2. Hub and shroud surfaces (blue/cyan, transparent)
3. Hex elements in passage (green wireframes)
4. Proper 3D perspective
"""

import json
import numpy as np
from pathlib import Path

def parse_mesh_passage(filepath):
    """Parse single blade passage from .msh file."""
    nodes = {}
    blade_faces = []
    hub_faces = []
    shroud_faces = []
    hex_elements = []
    
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
                    
                    # Blade surfaces (geom tags 5,6)
                    if geom_tag in [5, 6] and elem_type in [2, 3]:
                        face_coords = [nodes[nid] for nid in node_ids]
                        blade_faces.append(face_coords)
                    
                    # Hub surfaces (geom tags 3,4 - based on z=0)
                    elif geom_tag in [3, 4] and elem_type in [2, 3]:
                        face_coords = [nodes[nid] for nid in node_ids]
                        # Verify z=0
                        if all(abs(c[2]) < 0.01 for c in face_coords):
                            hub_faces.append(face_coords)
                    
                    # Shroud surfaces (geom tags 3,4 - based on z=2.5)
                    elif geom_tag in [3, 4] and elem_type in [2, 3]:
                        face_coords = [nodes[nid] for nid in node_ids]
                        # Verify z=2.5
                        if all(abs(c[2] - 2.5) < 0.01 for c in face_coords):
                            shroud_faces.append(face_coords)
                    
                    # Hex elements (sample)
                    elif elem_type == 5:
                        coords = [nodes[nid] for nid in node_ids]
                        hex_elements.append(coords)
            
            line = f.readline()
    
    return nodes, blade_faces, hub_faces, shroud_faces, hex_elements

def create_passage_viz_html(blade_faces, hex_elements):
    """Create HTML with proper passage visualization."""
    
    # Flatten blade surface
    blade_vertices = []
    blade_indices = []
    vertex_map = {}
    
    for face in blade_faces[:500]:  # Limit for performance
        face_indices = []
        for coord in face:
            key = tuple(round(c, 6) for c in coord)
            if key not in vertex_map:
                vertex_map[key] = len(blade_vertices)
                blade_vertices.append(coord)
            face_indices.append(vertex_map[key])
        
        if len(face) == 3:
            blade_indices.extend(face_indices)
        elif len(face) == 4:
            blade_indices.extend([face_indices[0], face_indices[1], face_indices[2]])
            blade_indices.extend([face_indices[0], face_indices[2], face_indices[3]])
    
    # Sample hex wireframes (every 50th)
    hex_samples = hex_elements[::50][:100]
    
    blade_data = {
        'vertices': blade_vertices,
        'indices': blade_indices
    }
    
    html = f'''<!DOCTYPE html>
<html>
<head>
    <title>T1_9 Blade Passage</title>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
    <style>
        body {{ margin: 0; overflow: hidden; background: #111; }}
        #info {{
            position: absolute; top: 10px; left: 10px;
            background: rgba(0,0,0,0.8); color: #fff;
            padding: 15px; border-radius: 8px; font-size: 13px;
            max-width: 320px; font-family: 'Segoe UI', Arial, sans-serif;
            border: 1px solid #333;
        }}
        #info h2 {{ margin: 0 0 12px 0; color: #4CAF50; font-size: 16px; }}
        .stat {{ margin: 4px 0; color: #ccc; }}
        .highlight {{ color: #ff9800; font-weight: bold; }}
        #controls {{
            position: absolute; bottom: 20px; left: 50%;
            transform: translateX(-50%);
            background: rgba(0,0,0,0.7); color: #fff;
            padding: 10px 20px; border-radius: 20px;
            font-family: Arial, sans-serif; font-size: 12px;
        }}
    </style>
</head>
<body>
    <div id="info">
        <h2>T1_9 Single Blade Passage</h2>
        <div class="stat">Mesh: <span class="highlight">44,800 hex</span> elements</div>
        <div class="stat">Passage: <span class="highlight">90°</span> sector (4 blades total)</div>
        <div class="stat">Z-range: <span class="highlight">0.0 to 2.5</span></div>
        <div class="stat">Radius: <span class="highlight">0.5 (hub) to 1.9 (shroud)</span></div>
        <hr style="border-color: #333; margin: 10px 0;">
        <div style="font-size: 11px; color: #888;">
            Red: Blade surface (pressure/suction)<br>
            Cyan: Hex core elements<br>
            Blue: Hub/Shroud reference planes
        </div>
    </div>
    
    <div id="controls">
        Left Click: Rotate | Right Click: Pan | Scroll: Zoom
    </div>

    <script>
        // Scene
        const scene = new THREE.Scene();
        scene.background = new THREE.Color(0x0a0a1a);
        
        const camera = new THREE.PerspectiveCamera(50, window.innerWidth/window.innerHeight, 0.1, 50);
        camera.position.set(3, 2, 4);
        
        const renderer = new THREE.WebGLRenderer({{ antialias: true }});
        renderer.setSize(window.innerWidth, window.innerHeight);
        renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
        document.body.appendChild(renderer.domElement);
        
        // Controls
        const controls = new THREE.OrbitControls(camera, renderer.domElement);
        controls.enableDamping = true;
        controls.dampingFactor = 0.05;
        controls.target.set(0.8, 0, 1.25);
        controls.minDistance = 1;
        controls.maxDistance = 10;
        
        // Lights
        const ambient = new THREE.AmbientLight(0x404040, 1.5);
        scene.add(ambient);
        
        const dir1 = new THREE.DirectionalLight(0xffffff, 1);
        dir1.position.set(5, 5, 5);
        scene.add(dir1);
        
        const dir2 = new THREE.DirectionalLight(0x8888ff, 0.5);
        dir2.position.set(-3, -2, 3);
        scene.add(dir2);
        
        // Blade surface
        const bladeData = {json.dumps(blade_data)};
        
        const bladeGeo = new THREE.BufferGeometry();
        const positions = new Float32Array(bladeData.vertices.length * 3);
        bladeData.vertices.forEach((v, i) => {{
            positions[i*3] = v[0];
            positions[i*3+1] = v[1];
            positions[i*3+2] = v[2];
        }});
        bladeGeo.setAttribute('position', new THREE.BufferAttribute(positions, 3));
        bladeGeo.setIndex(bladeData.indices);
        bladeGeo.computeVertexNormals();
        
        // Two materials for pressure/suction sides
        const bladeMat = new THREE.MeshPhongMaterial({{
            color: 0xe74c3c,
            side: THREE.DoubleSide,
            transparent: true,
            opacity: 0.85,
            shininess: 80,
            flatShading: false
        }});
        
        const blade = new THREE.Mesh(bladeGeo, bladeMat);
        scene.add(blade);
        
        // Blade wireframe
        const bladeWire = new THREE.WireframeGeometry(bladeGeo);
        const bladeLines = new THREE.LineSegments(
            bladeWire,
            new THREE.LineBasicMaterial({{ color: 0xff6b6b, transparent: true, opacity: 0.2 }})
        );
        scene.add(bladeLines);
        
        // Hub and Shroud planes
        const hubGeo = new THREE.RingGeometry(0.4, 2.0, 64);
        const hubMat = new THREE.MeshBasicMaterial({{
            color: 0x3498db, side: THREE.DoubleSide,
            transparent: true, opacity: 0.1
        }});
        const hub = new THREE.Mesh(hubGeo, hubMat);
        hub.rotation.x = -Math.PI / 2;
        scene.add(hub);
        
        const shroud = hub.clone();
        shroud.position.z = 2.5;
        scene.add(shroud);
        
        // Hex wireframes (subsampled)
        const hexData = {json.dumps([{'coords': h, 'edges': [[0,1],[1,2],[2,3],[3,0],[4,5],[5,6],[6,7],[7,4],[0,4],[1,5],[2,6],[3,7]]} for h in hex_samples])};
        
        hexData.forEach(hex => {{
            hex.edges.forEach(edge => {{
                const geo = new THREE.BufferGeometry().setFromPoints([
                    new THREE.Vector3(hex.coords[edge[0]][0], hex.coords[edge[0]][1], hex.coords[edge[0]][2]),
                    new THREE.Vector3(hex.coords[edge[1]][0], hex.coords[edge[1]][1], hex.coords[edge[1]][2])
                ]);
                const mat = new THREE.LineBasicMaterial({{
                    color: 0x00d2ff, transparent: true, opacity: 0.08
                }});
                const line = new THREE.Line(geo, mat);
                scene.add(line);
            }});
        }});
        
        // Coordinate system
        const axes = new THREE.AxesHelper(0.5);
        axes.position.set(-0.5, -0.5, 0);
        scene.add(axes);
        
        // Grid
        const grid = new THREE.GridHelper(4, 20, 0x222244, 0x111122);
        scene.add(grid);
        
        // Animation
        function animate() {{
            requestAnimationFrame(animate);
            controls.update();
            renderer.render(scene, camera);
        }}
        animate();
        
        // Resize
        window.addEventListener('resize', () => {{
            camera.aspect = window.innerWidth / window.innerHeight;
            camera.updateProjectionMatrix();
            renderer.setSize(window.innerWidth, window.innerHeight);
        }});
    </script>
</body>
</html>'''
    
    return html

def main():
    print("Parsing T1_9 blade passage...")
    nodes, blade_faces, hub_faces, shroud_faces, hex_elements = parse_mesh_passage(
        'T1_9/T1_9_ru_gridGmsh.msh'
    )
    
    print(f"Blade surface faces: {len(blade_faces)}")
    print(f"Hub faces: {len(hub_faces)}")
    print(f"Shroud faces: {len(shroud_faces)}")
    print(f"Hex elements: {len(hex_elements)}")
    
    # Check blade extent
    all_blade = np.array([c for face in blade_faces[:100] for c in face])
    print(f"\nBlade surface extent (sampled):")
    print(f"  X: [{all_blade[:,0].min():.3f}, {all_blade[:,0].max():.3f}]")
    print(f"  Y: [{all_blade[:,1].min():.3f}, {all_blade[:,1].max():.3f}]")
    print(f"  Z: [{all_blade[:,2].min():.3f}, {all_blade[:,2].max():.3f}]")
    
    print("\nGenerating visualization...")
    html = create_passage_viz_html(blade_faces, hex_elements)
    
    with open('turbine_passage_viz.html', 'w') as f:
        f.write(html)
    
    print("Saved to turbine_passage_viz.html")
    print("\nThis shows a SINGLE blade passage with:")
    print("  - Red: Blade pressure/suction surface")
    print("  - Cyan: Hex elements (subsampled)")
    print("  - Blue rings: Hub (z=0) and Shroud (z=2.5)")
    print("\nThe mesh contains only one 90° passage, not the full turbine.")

if __name__ == '__main__':
    main()
