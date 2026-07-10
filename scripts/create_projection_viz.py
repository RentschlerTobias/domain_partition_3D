#!/usr/bin/env python3
"""
Create visualization showing Hub singularities, blade curves, and projected Shroud positions.
Demonstrates why direct projection along blade curves is problematic.
"""

import json
import numpy as np
from pathlib import Path

def create_projection_viz_html():
    """Create HTML visualization of the projection problem."""
    
    # Load projection results
    with open('hub_to_shroud_projection.json', 'r') as f:
        proj_data = json.load(f)
    
    projections = proj_data['projections']
    
    # Prepare data for visualization
    hub_points = []
    shroud_points = []
    blade_curves = []
    
    for p in projections:
        sing = p['singularity']
        hub_points.append({
            'id': sing['id'],
            'pos': sing['position'],
            'index': sing['index']
        })
        
        shroud_points.append({
            'id': f"{sing['id']}_projected",
            'pos': p['projected_position'],
            'blade_point': p['shroud_blade_point']
        })
        
        # Blade curve from hub to shroud
        blade_curves.append({
            'hub': p['hub_blade_point'],
            'shroud': p['shroud_blade_point']
        })
    
    html = f'''<!DOCTYPE html>
<html>
<head>
    <title>Hub-to-Shroud Projection Problem</title>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
    <style>
        body {{ margin: 0; overflow: hidden; background: #1a1a2e; }}
        #info {{
            position: absolute; top: 10px; left: 10px;
            background: rgba(0,0,0,0.8); color: white;
            padding: 15px; border-radius: 5px;
            max-width: 400px; font-size: 12px;
            z-index: 100;
        }}
        #info h3 {{ margin: 0 0 10px 0; color: #ff6b6b; }}
        .warning {{ color: #ff9800; font-weight: bold; }}
        .stats {{ margin-top: 10px; padding-top: 10px; border-top: 1px solid #555; }}
        .problem {{ color: #ff6b6b; }}
    </style>
</head>
<body>
    <div id="info">
        <h3>Hub-to-Shroud Projection Problem</h3>
        <p class="warning">Direct projection along blade curves does NOT work well</p>
        
        <div class="stats">
            <strong>Issues:</strong><br>
            - Singularities are far from blade surface (0.26-0.40)<br>
            - Blade twist is 28-32°<br>
            - Multiple singularities project to same Shroud point<br>
            <br>
            <strong>Colors:</strong><br>
            Red dots: Hub singularities<br>
            Green dots: Shroud projected positions<br>
            Blue lines: Blade curves (hub→shroud)<br>
            Yellow lines: Projection paths
        </div>
    </div>

    <script>
        // Scene
        const scene = new THREE.Scene();
        scene.background = new THREE.Color(0x0a0a1a);
        
        const camera = new THREE.PerspectiveCamera(50, window.innerWidth/window.innerHeight, 0.1, 50);
        camera.position.set(2, 2, 3);
        
        const renderer = new THREE.WebGLRenderer({{ antialias: true }});
        renderer.setSize(window.innerWidth, window.innerHeight);
        document.body.appendChild(renderer.domElement);
        
        const controls = new THREE.OrbitControls(camera, renderer.domElement);
        controls.enableDamping = true;
        controls.target.set(0.5, 0.5, 1.25);
        
        // Lights
        const ambient = new THREE.AmbientLight(0x404040, 1);
        scene.add(ambient);
        const dir = new THREE.DirectionalLight(0xffffff, 1);
        dir.position.set(5, 5, 5);
        scene.add(dir);
        
        // Hub plane (z=0)
        const hubGeo = new THREE.PlaneGeometry(2, 2);
        const hubMat = new THREE.MeshBasicMaterial({{ color: 0x222244, side: THREE.DoubleSide, transparent: true, opacity: 0.3 }});
        const hub = new THREE.Mesh(hubGeo, hubMat);
        hub.rotation.x = -Math.PI / 2;
        scene.add(hub);
        
        // Shroud plane (z=2.5)
        const shroud = hub.clone();
        shroud.position.z = 2.5;
        scene.add(shroud);
        
        // Hub singularities (red)
        const hubPoints = {json.dumps(hub_points)};
        hubPoints.forEach(pt => {{
            const geom = new THREE.SphereGeometry(0.04, 8, 8);
            const mat = new THREE.MeshPhongMaterial({{ color: 0xff0000 }});
            const sphere = new THREE.Mesh(geom, mat);
            sphere.position.set(pt.pos[0], pt.pos[1], 0);
            scene.add(sphere);
            
            // Label
            const label = document.createElement('div');
            label.className = 'label';
            label.textContent = pt.id;
            // Add label to scene using canvas texture
        }});
        
        // Shroud projected points (green)
        const shroudPoints = {json.dumps(shroud_points)};
        shroudPoints.forEach(pt => {{
            const geom = new THREE.SphereGeometry(0.04, 8, 8);
            const mat = new THREE.MeshPhongMaterial({{ color: 0x00ff00 }});
            const sphere = new THREE.Mesh(geom, mat);
            sphere.position.set(pt.pos[0], pt.pos[1], 2.5);
            scene.add(sphere);
        }});
        
        // Blade curves (blue)
        const bladeCurves = {json.dumps(blade_curves)};
        bladeCurves.forEach(curve => {{
            const geom = new THREE.BufferGeometry().setFromPoints([
                new THREE.Vector3(curve.hub[0], curve.hub[1], curve.hub[2]),
                new THREE.Vector3(curve.shroud[0], curve.shroud[1], curve.shroud[2])
            ]);
            const mat = new THREE.LineBasicMaterial({{ color: 0x0000ff, linewidth: 2 }});
            const line = new THREE.Line(geom, mat);
            scene.add(line);
        }});
        
        // Projection paths (yellow dashed)
        for (let i = 0; i < hubPoints.length; i++) {{
            const hub = hubPoints[i].pos;
            const shroud = shroudPoints[i].pos;
            
            const geom = new THREE.BufferGeometry().setFromPoints([
                new THREE.Vector3(hub[0], hub[1], 0),
                new THREE.Vector3(shroud[0], shroud[1], 2.5)
            ]);
            const mat = new THREE.LineDashedMaterial({{
                color: 0xffff00,
                dashSize: 0.1,
                gapSize: 0.05,
                linewidth: 2
            }});
            const line = new THREE.Line(geom, mat);
            line.computeLineDistances();
            scene.add(line);
        }}
        
        // Axes
        const axes = new THREE.AxesHelper(1);
        scene.add(axes);
        
        // Grid
        const grid = new THREE.GridHelper(3, 10, 0x444444, 0x222222);
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
    print("Creating Hub-to-Shroud projection visualization...")
    html = create_projection_viz_html()
    
    with open('projection_problem_viz.html', 'w') as f:
        f.write(html)
    
    print("Saved to projection_problem_viz.html")
    print("\nThis visualization shows:")
    print("  - Red dots: Hub singularities")
    print("  - Green dots: Projected Shroud positions")
    print("  - Blue lines: Blade curves (hub→shroud)")
    print("  - Yellow dashed: Projection paths")
    print("\nThe problem is clear: singularities are NOT on the blade surface,")
    print("so projection along blade curves doesn't make geometric sense.")

if __name__ == '__main__':
    main()
