#!/usr/bin/env python3
"""
Parse T1_9_hex_boundaries.vtk and extract block corner points.
Create HTML wireframe with numbered labels.
"""

import numpy as np
from collections import defaultdict

def parse_vtk(filename):
    """Parse ASCII VTK unstructured grid file."""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    points = []
    cells = []
    cell_data = {}
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        if line.startswith('POINTS'):
            parts = line.split()
            num_points = int(parts[1])
            i += 1
            for _ in range(num_points):
                coords = [float(x) for x in lines[i].strip().split()]
                points.append(coords)
                i += 1
            continue
        
        elif line.startswith('CELLS'):
            parts = line.split()
            num_cells = int(parts[1])
            i += 1
            for _ in range(num_cells):
                cell_line = [int(x) for x in lines[i].strip().split()]
                n_nodes = cell_line[0]
                nodes = cell_line[1:]
                cells.append(nodes)
                i += 1
            continue
        
        elif line.startswith('CELL_TYPES'):
            num_types = int(line.split()[1])
            i += 1
            for _ in range(num_types):
                i += 1
            continue
        
        elif line.startswith('CELL_DATA'):
            num_data = int(line.split()[1])
            i += 1
            
            while i < len(lines) and lines[i].strip().startswith('SCALARS'):
                data_name = lines[i].split()[1]
                i += 2  # Skip LOOKUP_TABLE line
                
                values = []
                for _ in range(num_data):
                    values.append(int(lines[i].strip()))
                    i += 1
                
                cell_data[data_name] = values
            continue
        
        i += 1
    
    return np.array(points), cells, cell_data

def find_block_corners(points, cells, cell_data):
    """
    Find block corners by analyzing where boundary faces from different regions meet.
    A corner is where 3 or more boundary faces (from different regions) share a node.
    """
    # Build node-to-faces mapping
    node_to_faces = defaultdict(list)
    
    for face_idx, cell in enumerate(cells):
        for node in cell:
            node_to_faces[node].append(face_idx)
    
    # For each node, count how many different regions it touches
    node_region_info = {}
    
    for node, faces in node_to_faces.items():
        regions = set()
        for face_idx in faces:
            if 'region_id' in cell_data and face_idx < len(cell_data['region_id']):
                regions.add(cell_data['region_id'][face_idx])
        
        node_region_info[node] = {
            'num_faces': len(faces),
            'regions': sorted(regions),
            'num_regions': len(regions)
        }
    
    # Corners are nodes that:
    # 1. Belong to 3 or more different regions (corner of block)
    # OR
    # 2. Have high valence (many boundary faces meeting)
    
    corners = []
    for node, info in node_region_info.items():
        if info['num_regions'] >= 3 or info['num_faces'] >= 6:
            corners.append({
                'node_id': node,
                'position': points[node].tolist(),
                'num_faces': info['num_faces'],
                'num_regions': info['num_regions'],
                'regions': info['regions']
            })
    
    # Sort by number of regions (most connected first), then by position
    corners.sort(key=lambda x: (-x['num_regions'], -x['num_faces'], x['position'][2], x['position'][0]))
    
    return corners

def create_html_wireframe(points, cells, corners, output_file='block_corners_wireframe.html'):
    """Create HTML/Three.js wireframe with numbered labels."""
    
    # Prepare corner data for JavaScript
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
    
    # Prepare edges (wireframe lines from cells)
    edges = set()
    for cell in cells:
        if len(cell) == 4:  # Quad face
            # 4 edges per quad
            for i in range(4):
                n1 = cell[i]
                n2 = cell[(i+1) % 4]
                edge = tuple(sorted([n1, n2]))
                edges.add(edge)
    
    # Convert to lines for Three.js
    lines_js = []
    for edge in edges:
        n1, n2 = edge
        p1 = points[n1]
        p2 = points[n2]
        lines_js.append({
            'x1': float(p1[0]), 'y1': float(p1[1]), 'z1': float(p1[2]),
            'x2': float(p2[0]), 'y2': float(p2[1]), 'z2': float(p2[2])
        })
    
    # Compute bounding box for camera
    bb_min = points.min(axis=0)
    bb_max = points.max(axis=0)
    center = (bb_min + bb_max) / 2
    size = np.linalg.norm(bb_max - bb_min)
    
    html = f'''<!DOCTYPE html>
<html>
<head>
    <title>T1_9 Hex Block Corners Wireframe</title>
    <style>
        body {{ margin: 0; overflow: hidden; font-family: Arial, sans-serif; }}
        #info {{
            position: absolute; top: 10px; left: 10px;
            background: rgba(0,0,0,0.7); color: white;
            padding: 15px; border-radius: 5px;
            max-width: 400px; font-size: 14px;
        }}
        #corner-list {{
            position: absolute; top: 10px; right: 10px;
            background: rgba(0,0,0,0.7); color: white;
            padding: 15px; border-radius: 5px;
            max-width: 300px; max-height: 90vh;
            overflow-y: auto; font-size: 12px;
        }}
        .corner-item {{ margin: 2px 0; padding: 2px 5px; cursor: pointer; }}
        .corner-item:hover {{ background: rgba(255,255,255,0.2); }}
        .highlight {{ background: rgba(255,255,0,0.3) !important; }}
    </style>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
</head>
<body>
    <div id="info">
        <h3>T1_9 Hex Block Corners</h3>
        <p><b>Total Corners:</b> {len(corners)}</p>
        <p><b>Lines:</b> {len(lines_js)}</p>
        <p>Hover over corner list to highlight</p>
        <p>Left click: Rotate | Right click: Pan | Scroll: Zoom</p>
    </div>
    
    <div id="corner-list">
        <h4>Corners (sorted by connectivity)</h4>
        <div id="corner-items"></div>
    </div>

    <script>
        // Data
        const corners = {json.dumps(corners_js, indent=2)};
        const lines = {json.dumps(lines_js[:5000], indent=2)}; // Limit lines for performance
        
        // Scene setup
        const scene = new THREE.Scene();
        scene.background = new THREE.Color(0x1a1a2e);
        
        const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 100);
        camera.position.set({center[0] + size*0.8}, {center[1] + size*0.5}, {center[2] + size*1.2});
        
        const renderer = new THREE.WebGLRenderer({{ antialias: true }});
        renderer.setSize(window.innerWidth, window.innerHeight);
        document.body.appendChild(renderer.domElement);
        
        const controls = new THREE.OrbitControls(camera, renderer.domElement);
        controls.target.set({center[0]}, {center[1]}, {center[2]});
        controls.update();
        
        // Lighting
        const ambientLight = new THREE.AmbientLight(0xffffff, 0.6);
        scene.add(ambientLight);
        const directionalLight = new THREE.DirectionalLight(0xffffff, 0.4);
        directionalLight.position.set(1, 1, 1);
        scene.add(directionalLight);
        
        // Wireframe lines
        const lineMaterial = new THREE.LineBasicMaterial({{ 
            color: 0x00ff88, 
            opacity: 0.6, 
            transparent: true 
        }});
        
        const lineGeometry = new THREE.BufferGeometry();
        const positions = new Float32Array(lines.length * 6);
        
        for (let i = 0; i < lines.length; i++) {{
            const l = lines[i];
            positions[i * 6] = l.x1;
            positions[i * 6 + 1] = l.y1;
            positions[i * 6 + 2] = l.z1;
            positions[i * 6 + 3] = l.x2;
            positions[i * 6 + 4] = l.y2;
            positions[i * 6 + 5] = l.z2;
        }}
        
        lineGeometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
        const lineObject = new THREE.LineSegments(lineGeometry, lineMaterial);
        scene.add(lineObject);
        
        // Corner points and labels
        const cornerGroup = new THREE.Group();
        const labelSprites = [];
        
        corners.forEach((corner, idx) => {{
            // Point
            const geometry = new THREE.SphereGeometry(0.015 * {size}, 16, 16);
            const material = new THREE.MeshPhongMaterial({{ 
                color: corner.num_regions >= 4 ? 0xff4444 : 0xffaa00,
                emissive: corner.num_regions >= 4 ? 0xff0000 : 0xff6600,
                emissiveIntensity: 0.5
            }});
            const sphere = new THREE.Mesh(geometry, material);
            sphere.position.set(corner.x, corner.y, corner.z);
            sphere.userData = {{ cornerId: idx }};
            cornerGroup.add(sphere);
            
            // Label
            const canvas = document.createElement('canvas');
            const context = canvas.getContext('2d');
            canvas.width = 128;
            canvas.height = 64;
            context.fillStyle = 'rgba(0, 0, 0, 0.7)';
            context.fillRect(0, 0, 128, 64);
            context.strokeStyle = 'white';
            context.lineWidth = 2;
            context.strokeRect(0, 0, 128, 64);
            context.font = 'bold 24px Arial';
            context.fillStyle = 'white';
            context.textAlign = 'center';
            context.fillText(idx.toString(), 64, 40);
            
            const texture = new THREE.CanvasTexture(canvas);
            const spriteMaterial = new THREE.SpriteMaterial({{ map: texture }});
            const sprite = new THREE.Sprite(spriteMaterial);
            sprite.position.set(corner.x + 0.05 * {size}, corner.y + 0.05 * {size}, corner.z + 0.05 * {size});
            sprite.scale.set(0.15 * {size}, 0.075 * {size}, 1);
            sprite.userData = {{ cornerId: idx }};
            
            cornerGroup.add(sprite);
            labelSprites.push(sprite);
        }});
        
        scene.add(cornerGroup);
        
        // Build corner list HTML
        const cornerItemsDiv = document.getElementById('corner-items');
        corners.forEach((corner, idx) => {{
            const div = document.createElement('div');
            div.className = 'corner-item';
            div.innerHTML = `
                <b>Corner ${{idx}}</b> (Node ${{corner.node_id}})<br>
                Regions: ${{corner.regions.join(', ')}}<br>
                Pos: (${{corner.x.toFixed(3)}}, ${{corner.y.toFixed(3)}}, ${{corner.z.toFixed(3)}})
            `;
            div.addEventListener('mouseenter', () => {{
                div.classList.add('highlight');
                // Highlight 3D point
                cornerGroup.children.forEach(child => {{
                    if (child.userData.cornerId === idx) {{
                        if (child.material.emissive) {{
                            child.material.emissiveIntensity = 1.0;
                        }}
                    }}
                }});
            }});
            div.addEventListener('mouseleave', () => {{
                div.classList.remove('highlight');
                cornerGroup.children.forEach(child => {{
                    if (child.userData.cornerId === idx) {{
                        if (child.material.emissive) {{
                            child.material.emissiveIntensity = 0.5;
                        }}
                    }}
                }});
            }});
            cornerItemsDiv.appendChild(div);
        }});
        
        // Animation loop
        function animate() {{
            requestAnimationFrame(animate);
            controls.update();
            renderer.render(scene, camera);
        }}
        
        // Handle window resize
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
    print(f"  Corners: {len(corners)}")
    print(f"  Lines: {len(lines_js)}")
    print(f"  Bounding box: [{bb_min[0]:.3f}, {bb_min[1]:.3f}, {bb_min[2]:.3f}] -> [{bb_max[0]:.3f}, {bb_max[1]:.3f}, {bb_max[2]:.3f}]")

def main():
    vtk_file = 'T1_9_hex_boundaries.vtk'
    
    print("Parsing VTK boundaries...")
    points, cells, cell_data = parse_vtk(vtk_file)
    print(f"Loaded {len(points)} points, {len(cells)} boundary faces")
    
    print("\nFinding block corners...")
    corners = find_block_corners(points, cells, cell_data)
    
    print(f"\nFound {len(corners)} block corners:")
    for i, corner in enumerate(corners[:20]):  # Print first 20
        print(f"  Corner {i}: Node {corner['node_id']}, "
              f"Regions {corner['regions']}, "
              f"Faces={corner['num_faces']}, "
              f"Pos=({corner['position'][0]:.3f}, {corner['position'][1]:.3f}, {corner['position'][2]:.3f})")
    
    if len(corners) > 20:
        print(f"  ... and {len(corners)-20} more")
    
    # Create HTML
    create_html_wireframe(points, cells, corners, 'block_corners_wireframe.html')
    
    # Save corners to JSON
    import json
    with open('block_corners.json', 'w') as f:
        json.dump(corners, f, indent=2)
    print("\nCorners saved to block_corners.json")

if __name__ == '__main__':
    main()
