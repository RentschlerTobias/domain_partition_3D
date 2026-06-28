#!/usr/bin/env python3
"""
Process T1_9 mesh in dtOO Docker container:
1. Load T1_9 mesh from /work/T1_9/T1_9_ru_gridGmsh.msh
2. Identify and tag surfaces (hub, shroud, blade, inlet, outlet, periodic)
3. Apply quasi-structured quad algorithm (Algorithm 11) to hub/shroud
4. Export coarse block structure
"""
import sys
import os

# Add gmsh package from workspace
sys.path.insert(0, '/work/gmsh_pkg')

# Import gmsh
import gmsh
import numpy as np

# Set OpenFOAM environment
os.environ['LD_LIBRARY_PATH'] = '/dtOO-install/lib:/usr/lib/openfoam/openfoam2406/platforms/linux64GccDPInt32Opt/lib/sys-openmpi:/usr/lib64/mpi/gcc/openmpi4/lib64:/usr/lib/hpc/gnu7/mpi/openmpi/4.1.6/lib64'

def identify_and_tag_surfaces():
    """
    Identify surfaces by their element count and tag them with physical groups.
    """
    surfaces = gmsh.model.getEntities(2)
    print(f"Found {len(surfaces)} surfaces")
    
    # Get surface info
    surface_info = {}
    for dim, tag in surfaces:
        elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(2, tag)
        total_elements = sum(len(t) for t in elem_tags)
        
        # Get center point
        bbox = gmsh.model.getBoundingBox(dim, tag)
        center = [
            (bbox[0] + bbox[3]) / 2,
            (bbox[1] + bbox[4]) / 2,
            (bbox[2] + bbox[5]) / 2
        ]
        
        surface_info[tag] = {
            'elements': total_elements,
            'center': center,
            'bbox': bbox
        }
        print(f"  Surface {tag}: {total_elements} elements, center={center}")
    
    # Sort by Z-coordinate and element count to identify surfaces
    # Hub surfaces typically have smaller Z range, shroud larger Z
    # Blade is between them
    
    sorted_by_z = sorted(surface_info.items(), key=lambda x: x[1]['center'][2])
    
    print("\nSurfaces sorted by Z-coordinate:")
    for tag, info in sorted_by_z:
        print(f"  Tag {tag}: Z={info['center'][2]:.3f}, elements={info['elements']}")
    
    # Identify surfaces based on characteristics
    # Based on previous analysis, T1_9 has:
    # - hub_0..5 (tags 6-11): quads, boundary layer patches
    # - shroud_0..5 (tags 12-17): quads, boundary layer patches
    # - blade: triangles
    # - inlet, outlet, periodic
    
    # For now, create physical groups based on z-position and element count
    hub_tags = []
    shroud_tags = []
    blade_tags = []
    inlet_tags = []
    outlet_tags = []
    periodic_tags = []
    
    for tag, info in surface_info.items():
        z = info['center'][2]
        elements = info['elements']
        
        # Very low Z (z < 0.5): hub
        if z < 0.5:
            hub_tags.append(tag)
        # Very high Z (z > 2.0): shroud
        elif z > 2.0:
            shroud_tags.append(tag)
        # Middle Z with many elements: blade
        elif elements > 1000:
            blade_tags.append(tag)
        # Other: periodic/inlet/outlet
        else:
            # Check X/Y extent to distinguish
            dx = info['bbox'][3] - info['bbox'][0]
            dy = info['bbox'][4] - info['bbox'][1]
            
            if dx > 1.0 or dy > 1.0:
                periodic_tags.append(tag)
            else:
                # Distinguish inlet/outlet by X position
                if info['center'][0] < 0.5:
                    inlet_tags.append(tag)
                else:
                    outlet_tags.append(tag)
    
    # Create physical groups
    if hub_tags:
        gmsh.model.addPhysicalGroup(2, hub_tags, tag=1)
        gmsh.model.setPhysicalName(2, 1, "RU_HUB")
        print(f"\nRU_HUB: surfaces {hub_tags}")
    
    if shroud_tags:
        gmsh.model.addPhysicalGroup(2, shroud_tags, tag=2)
        gmsh.model.setPhysicalName(2, 2, "RU_SHROUD")
        print(f"RU_SHROUD: surfaces {shroud_tags}")
    
    if blade_tags:
        gmsh.model.addPhysicalGroup(2, blade_tags, tag=3)
        gmsh.model.setPhysicalName(2, 3, "RU_BLADE")
        print(f"RU_BLADE: surfaces {blade_tags}")
    
    if inlet_tags:
        gmsh.model.addPhysicalGroup(2, inlet_tags, tag=4)
        gmsh.model.setPhysicalName(2, 4, "RU_INLET")
        print(f"RU_INLET: surfaces {inlet_tags}")
    
    if outlet_tags:
        gmsh.model.addPhysicalGroup(2, outlet_tags, tag=5)
        gmsh.model.setPhysicalName(2, 5, "RU_OUTLET")
        print(f"RU_OUTLET: surfaces {outlet_tags}")
    
    if periodic_tags:
        gmsh.model.addPhysicalGroup(2, periodic_tags, tag=6)
        gmsh.model.setPhysicalName(2, 6, "RU_PERIODIC")
        print(f"RU_PERIODIC: surfaces {periodic_tags}")
    
    return hub_tags, shroud_tags, blade_tags

def remesh_hub_shroud(hub_tags, shroud_tags, element_size=0.5):
    """
    Apply quasi-structured quad algorithm to hub and shroud surfaces.
    """
    # Set algorithm to quasi-structured quad (11)
    gmsh.option.setNumber("Mesh.Algorithm", 11)
    
    # Set element size
    for tag in hub_tags + shroud_tags:
        gmsh.model.mesh.setSize([(2, tag)], element_size)
    
    # Clear mesh
    gmsh.model.mesh.clear()
    
    # Generate 2D mesh
    gmsh.model.mesh.generate(2)
    
    # Get quad elements
    _, quad_tags, _ = gmsh.model.mesh.getElementsByType(3)
    print(f"\nGenerated {len(quad_tags)} quadrilateral elements")
    
    return quad_tags

def extract_coarse_blocks():
    """
    Extract coarse block structure from the quad mesh.
    """
    # Get all nodes
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    vertices = coords.reshape(-1, 3)
    
    # Get quad elements
    _, quad_tags, _ = gmsh.model.mesh.getElementsByType(3)
    
    print(f"\nExtracting coarse blocks from {len(quad_tags)} quads...")
    
    # For now, just return the info
    # Block extraction would require graph analysis
    return vertices, quad_tags

def main():
    # Initialize Gmsh
    gmsh.initialize()
    
    # Load T1_9 mesh
    mesh_path = "/work/T1_9/T1_9_ru_gridGmsh.msh"
    print(f"Loading mesh from {mesh_path}...")
    gmsh.open(mesh_path)
    
    # Identify and tag surfaces
    print("\n=== Identifying Surfaces ===")
    hub_tags, shroud_tags, blade_tags = identify_and_tag_surfaces()
    
    # Export mesh with physical groups
    output_path = "/work/T1_9_tagged_surfaces.msh"
    print(f"\nExporting tagged mesh to {output_path}...")
    gmsh.write(output_path)
    
    # Try to remesh hub and shroud with available algorithms
    print("\n=== Remeshing Hub/Shroud ===")
    try:
        # Algorithm 11 requires QUADMESHINGTOOLS, try Algorithm 8 instead
        # Algorithm 8: Frontal Quad (quad-dominant)
        gmsh.option.setNumber("Mesh.Algorithm", 8)
        
        # Set element size
        for tag in hub_tags + shroud_tags:
            gmsh.model.mesh.setSize([(2, tag)], 0.5)
        
        # Clear mesh and regenerate
        gmsh.model.mesh.clear()
        gmsh.model.mesh.generate(2)
        
        # Get quad elements
        _, quad_tags, _ = gmsh.model.mesh.getElementsByType(3)
        print(f"Generated {len(quad_tags)} quadrilateral elements")
        
        # Export remeshed surfaces
        output_path = "/work/T1_9_remeshed_quad8.msh"
        print(f"\nExporting remeshed mesh to {output_path}...")
        gmsh.write(output_path)
    except Exception as e:
        print(f"Remeshing failed: {e}")
        print("Skipping remeshing, exporting original mesh with tags...")
    
    # Finalize
    gmsh.finalize()
    print("\nDone!")

if __name__ == '__main__':
    main()
