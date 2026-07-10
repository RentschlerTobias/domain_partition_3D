#!/usr/bin/env python3
"""
Create a proper turbine visualization showing actual mesh geometry.
This extracts blade surface elements and a subset of hex elements
with correct coordinates from the .msh file.
"""

import json
import numpy as np
from pathlib import Path

def parse_msh_for_visualization(filepath):
    """Parse mesh with focus on visualization elements."""
    nodes = {}
    blade_faces = []  # Triangles and quads on blade surface
    hex_samples = []  # Sample of hex elements
    
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
                count = 0
                for _ in range(n_elements):
                    parts = f.readline().strip().split()
                    elem_type = int(parts[1])
                    n_tags = int(parts[2])
                    tags = [int(p) for p in parts[3:3+n_tags]]
                    node_ids = [int(p) for p in parts[3+n_tags:]]
                    
                    geom_tag = tags[1] if len(tags) > 1 else 0
                    
                    # Collect blade surface elements (geom tags 5,6)
                    if geom_tag in [5, 6]:
                        if elem_type in [2, 3]:  # triangles or quads
                            face_coords = [nodes[nid] for nid in node_ids]
                            blade_faces.append({
                                'type': 'tri' if elem_type == 2 else 'quad',
                                'nodes': node_ids,
                                'coords': face_coords,
                                'geom_tag': geom_tag
                            })
                    
                    # Sample hex elements (every 20th)
                    elif elem_type == 5 and count % 20 == 0:
                        coords = [nodes[nid] for nid in node_ids]
                        hex_samples.append({
                            'nodes': node_ids,
                            'coords': coords,
                            'center': np.mean(coords, axis=0).tolist()
                        })
                        count += 1
                    elif elem_type == 5:
                        count += 1
            
            line = f.readline()
    
    return nodes, blade_faces, hex_samples

def create_html_visualization(nodes, blade_faces, hex_samples):
    """Create standalone HTML with embedded Three.js visualization."""
    
    # Prepare blade surface geometry
    blade_vertices = []
    blade_indices = []
    vertex_map = {}
    
    for face in blade_faces:
        face_indices = []
        for coord in face['coords']:
            key = tuple(round(c, 6) for c in coord)
            if key not in vertex_map:
                vertex_map[key] = len(blade_vertices)
                blade_vertices.append(coord)
            face_indices.append(vertex_map[key])
        
        if face['type'] == 'tri':
            blade_indices.extend(face_indices)
        elif face['type'] == 'quad':
            # Split quad into two triangles
            blade_indices.extend([face_indices[0], face_indices[1], face_indices[2]])
            blade_indices.extend([face_indices[0], face_indices[2], face_indices[3]])
    
    # Sample hex wireframes
    hex_wireframes = []
    for hex_elem in hex_samples[:100]:  # Limit to 100 for performance
        coords = hex_elem['coords']
        # Hex edges
        edges = [
            [0,1], [1,2], [2,3], [3,0],  # Bottom face
            [4,5], [5,6], [6,7], [7,4],  # Top face  
            [0,4], [1,5], [2,6], [3,7]   # Vertical
        ]
        hex_wireframes.append({
            'coords': coords,
            'edges': edges
        })
    
    # Convert to JSON for embedding
    blade_data = {
        'vertices': blade_vertices,
        'indices': blade_indices
    }
    
    html_content = f'''<!DOCTYPE html>
<html>
<head>
    <title>T1_9 Turbine - Actual Mesh Geometry</title>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
    <style>
        body {{ margin: 0; overflow: hidden; font-family: Arial, sans-serif; }}
        #info {{
            position: absolute; top: 10px; left: 10px;
            background: rgba(0,0,0,0.8); color: white;
            padding: 15px; border-radius: 5px;
            max-width: 350px; font-size: 12px;
            z-index: 100;
        }}
        #info h3 {{ margin: 0 0 10px 0; color: #4CAF50; }}
        .legend-item {{ display: flex; align-items: center; margin: 5px 0; }}
        .color-box {{ width: 20px; height: 12px; margin-right: 8px; border: 1px solid white; }}
        .stats {{ margin-top: 10px; padding-top: 10px; border-top: 1px solid #555; }}
    </style>
</head>
<body>
    <div id="info">
        <h3>T1_9 Turbine Mesh Geometry</h3>
        <div class="legend-item">
            <div class="color-box" style="background: #ff6b6b;"></div>
            <span>Blade Surface (Pressure/Suction)</span>
        </div>
        <div class="legend-item">
            <div class="color-box" style="background: #4ecdc4;"></div>
            <span>Blade Surface (Hub/Shroud)</span>
        </div>
        <div class="legend-item">
            <div class="color-box" style="background: #95e1d3;"></div>
            <span>Hex Core (subsampled wireframes)</span>
        </div>
        <div class="stats">
            <strong>Mesh Statistics:</strong><br>
            Nodes: 130,198<br>
            Hex Elements: 44,800<br>
            Blade Triangles: 2,828<br>
            Blade Quads: 1,728<br>
            Hub Radius: 0.5<br>
            Shroud Radius: 1.9
        </div>
        <p><strong>Controls:</strong><br>
        Left click: Rotate | Right click: Pan | Scroll: Zoom</p>
    </div>

    <script>
        // Scene setup
        const scene = new THREE.Scene();
        scene.background = new THREE.Color(0x1a1a2e);
        
        const camera = new THREE.PerspectiveCamera(60, window.innerWidth / window.innerHeight, 0.1, 100);
        camera.position.set(2, 2, 3);
        
        const renderer = new THREE.WebGLRenderer({{ antialias: true }});
        renderer.setSize(window.innerWidth, window.innerHeight);
        renderer.setPixelRatio(window.devicePixelRatio);
        document.body.appendChild(renderer.domElement);
        
        const controls = new THREE.OrbitControls(camera, renderer.domElement);
        controls.enableDamping = true;
        controls.dampingFactor = 0.05;
        controls.target.set(0, 0, 1.25);
        
        // Lighting
        const ambientLight = new THREE.AmbientLight(0x404040, 0.6);
        scene.add(ambientLight);
        const dirLight1 = new THREE.DirectionalLight(0xffffff, 0.8);
        dirLight1.position.set(5, 5, 5);
        scene.add(dirLight1);
        const dirLight2 = new THREE.DirectionalLight(0xffffff, 0.3);
        dirLight2.position.set(-5, -5, 5);
        scene.add(dirLight2);
        
        // Blade surface geometry (embedded from mesh data)
        const bladeData = {json.dumps(blade_data)};
        
        const bladeGeometry = new THREE.BufferGeometry();
        const bladePositions = new Float32Array(bladeData.vertices.length * 3);
        bladeData.vertices.forEach((v, i) => {{
            bladePositions[i*3] = v[0];
            bladePositions[i*3+1] = v[1];
            bladePositions[i*3+2] = v[2];
        }});
        bladeGeometry.setAttribute('position', new THREE.BufferAttribute(bladePositions, 3));
        bladeGeometry.setIndex(bladeData.indices);
        bladeGeometry.computeVertexNormals();
        
        const bladeMaterial = new THREE.MeshPhongMaterial({{
            color: 0xff6b6b,
            side: THREE.DoubleSide,
            transparent: true,
            opacity: 0.9,
            shininess: 100
        }});
        const bladeMesh = new THREE.Mesh(bladeGeometry, bladeMaterial);
        scene.add(bladeMesh);
        
        // Add wireframe overlay for blade
        const bladeWireframe = new THREE.WireframeGeometry(bladeGeometry);
        const bladeLine = new THREE.LineSegments(
            bladeWireframe, 
            new THREE.LineBasicMaterial({{ color: 0xffaaaa, transparent: true, opacity: 0.3 }})
        );
        scene.add(bladeLine);
        
        // Hex elements as wireframes (subsampled)
        const hexData = {json.dumps(hex_wireframes)};
        
        hexData.forEach(hex => {{
            hex.edges.forEach(edge => {{
                const geometry = new THREE.BufferGeometry().setFromPoints([
                    new THREE.Vector3(hex.coords[edge[0]][0], hex.coords[edge[0]][1], hex.coords[edge[0]][2]),
                    new THREE.Vector3(hex.coords[edge[1]][0], hex.coords[edge[1]][1], hex.coords[edge[1]][2])
                ]);
                const material = new THREE.LineBasicMaterial({{
                    color: 0x4ecdc4, 
                    transparent: true, 
                    opacity: 0.15,
                    linewidth: 1
                }});
                const line = new THREE.Line(geometry, material);
                scene.add(line);
            }});
        }});
        
        // Add hub and shroud planes (approximate)
        const hubGeometry = new THREE.RingGeometry(0.3, 0.7, 32);
        const hubMaterial = new THREE.MeshBasicMaterial({{
            color: 0x4444ff, 
            side: THREE.DoubleSide, 
            transparent: true, 
            opacity: 0.1
        }});
        const hubPlane = new THREE.Mesh(hubGeometry, hubMaterial);
        hubPlane.rotation.x = -Math.PI / 2;
        scene.add(hubPlane);
        
        const shroudGeometry = new THREE.RingGeometry(1.7, 2.1, 32);
        const shroudPlane = new THREE.Mesh(shroudGeometry, hubMaterial);
        shroudPlane.rotation.x = -Math.PI / 2;
        shroudPlane.position.y = 2.5;
        scene.add(shroudPlane);
        
        // Coordinate axes
        const axesHelper = new THREE.AxesHelper(1);
        scene.add(axesHelper);
        
        // Grid
        const gridHelper = new THREE.GridHelper(4, 20, 0x444444, 0x222222);
        scene.add(gridHelper);
        
        // Animation
        function animate() {{
            requestAnimationFrame(animate);
            controls.update();
            renderer.render(scene, camera);
        }}
        animate();
        
        // Resize handler
        window.addEventListener('resize', () => {{
            camera.aspect = window.innerWidth / window.innerHeight;
            camera.updateProjectionMatrix();
            renderer.setSize(window.innerWidth, window.innerHeight);
        }});
    </script>
</body>
</html>'''
    
    return html_content

def main():
    print("Parsing mesh file for visualization...")
    msh_path = Path('T1_9/T1_9_ru_gridGmsh.msh')
    nodes, blade_faces, hex_samples = parse_msh_for_visualization(msh_path)
    
    print(f"Found {len(blade_faces)} blade surface faces")
    print(f"Found {len(hex_samples)} hex samples")
    
    # Check coordinate ranges
    blade_vertices = []
    for face in blade_faces:
        blade_vertices.extend(face['coords'])
    
    blade_array = np.array(blade_vertices)
    print(f"\nBlade coordinate ranges:")
    print(f"  X: [{blade_array[:,0].min():.3f}, {blade_array[:,0].max():.3f}]")
    print(f"  Y: [{blade_array[:,1].min():.3f}, {blade_array[:,1].max():.3f}]")
    print(f"  Z: [{blade_array[:,2].min():.3f}, {blade_array[:,2].max():.3f}]")
    
    print("\nGenerating HTML visualization...")
    html_content = create_html_visualization(nodes, blade_faces, hex_samples)
    
    output_path = Path('turbine_actual_geometry.html')
    with open(output_path, 'w') as f:
        f.write(html_content)
    
    print(f"Saved to {output_path}")
    print(f"File size: {output_path.stat().st_size / 1024:.1f} KB")
    print("\nThis file contains the ACTUAL mesh geometry embedded directly.")
    print("Open it in a browser to see the real turbine shape.")

if __name__ == '__main__':
    main()
