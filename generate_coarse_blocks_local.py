#!/usr/bin/env python3
"""
Generate coarse block-structured mesh for T1_9 hub/shroud surfaces.
Uses local Gmsh Python API with full Algorithm 11 support.
"""
import sys
sys.path.insert(0, '/root/venv/lib/python3.12/site-packages')

import gmsh
import numpy as np

def identify_surfaces_by_physical_names():
    """Identify surfaces by their physical names in the mesh."""
    # Get all physical groups
    phys_groups = gmsh.model.getPhysicalGroups(2)
    
    hub_tags = []
    shroud_tags = []
    blade_tags = []
    inlet_tags = []
    outlet_tags = []
    periodic_tags = []
    
    print("Physical groups found:")
    for dim, tag in phys_groups:
        name = gmsh.model.getPhysicalName(dim, tag)
        entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
        print(f"  {name} (tag={tag}): entities {entities}")
        
        if 'hub' in name.lower():
            hub_tags.extend(entities)
        elif 'shroud' in name.lower():
            shroud_tags.extend(entities)
        elif 'blade' in name.lower():
            blade_tags.extend(entities)
        elif 'inlet' in name.lower():
            inlet_tags.extend(entities)
        elif 'outlet' in name.lower():
            outlet_tags.extend(entities)
        elif 'period' in name.lower():
            periodic_tags.extend(entities)
    
    return hub_tags, shroud_tags, blade_tags, inlet_tags, outlet_tags, periodic_tags

def get_surface_elements(surface_tags):
    """Get all elements for given surface tags."""
    all_elements = []
    for tag in surface_tags:
        elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(2, tag)
        for etype, tags, nodes in zip(elem_types, elem_tags, elem_node_tags):
            if etype == 3:  # Quads
                all_elements.extend(tags)
    return all_elements

def remesh_surface_with_quads(surface_tags, element_size=0.5):
    """Remesh surfaces using quasi-structured quad algorithm (Algorithm 11)."""
    # Set Algorithm 11: Quasi-structured Quad
    gmsh.option.setNumber("Mesh.Algorithm", 11)
    
    # Set element size for the surfaces
    all_points = []
    for tag in surface_tags:
        # Get boundary points of surface
        boundary = gmsh.model.getBoundary([(2, tag)], recursive=True)
        for bdim, btag in boundary:
            if bdim == 0:
                all_points.append((0, btag))
    
    # Remove duplicates
    all_points = list(set(all_points))
    
    # Set size
    gmsh.model.mesh.setSize(all_points, element_size)
    
    # Clear mesh for these surfaces
    gmsh.model.mesh.clear()
    
    # Generate 2D mesh
    gmsh.model.mesh.generate(2)
    
    # Get quads
    _, quad_tags, _ = gmsh.model.mesh.getElementsByType(3)
    
    return quad_tags

def extract_coarse_block_structure(quad_tags, surface_tags):
    """Extract coarse block structure from quad mesh."""
    # Get nodes
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    vertices = coords.reshape(-1, 3)
    
    # Build adjacency graph
    import networkx as nx
    G = nx.Graph()
    
    # Add nodes
    for i, tag in enumerate(node_tags):
        G.add_node(int(tag), pos=vertices[i])
    
    # Add edges from quads
    _, _, node_tags_quad = gmsh.model.mesh.getElementsByType(3)
    node_tags_quad = node_tags_quad.reshape(-1, 4)
    
    for quad in node_tags_quad:
        for i in range(4):
            G.add_edge(int(quad[i]), int(quad[(i+1)%4]))
    
    # Find coarse blocks (faces with 4 edges)
    # This is a simplified approach - we need to find 4-sided faces
    blocks = []
    
    # For each surface, try to find 4-sided blocks
    for surf_tag in surface_tags:
        # Get elements for this surface
        elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(2, surf_tag)
        
        for etype, tags, nodes in zip(elem_types, elem_tags, elem_node_tags):
            if etype == 3:  # Quads
                # Each quad is a potential block
                nodes = nodes.reshape(-1, 4)
                for quad in nodes:
                    blocks.append(quad)
    
    return blocks

def main():
    # Initialize Gmsh
    gmsh.initialize()
    
    # Load T1_9 mesh
    mesh_path = "/root/repos/block_structured_meshing/T1_9/T1_9_ru_gridGmsh.msh"
    print(f"Loading mesh from {mesh_path}...")
    gmsh.open(mesh_path)
    
    # Identify surfaces by physical names
    print("\n=== Identifying Surfaces ===")
    hub_tags, shroud_tags, blade_tags, inlet_tags, outlet_tags, periodic_tags = identify_surfaces_by_physical_names()
    
    print(f"\nHub surfaces: {hub_tags}")
    print(f"Shroud surfaces: {shroud_tags}")
    print(f"Blade surfaces: {blade_tags}")
    
    # Get current elements on hub and shroud
    hub_elements = get_surface_elements(hub_tags)
    shroud_elements = get_surface_elements(shroud_tags)
    
    print(f"\nCurrent hub elements: {len(hub_elements)} quads")
    print(f"Current shroud elements: {len(shroud_elements)} quads")
    
    # Try to remesh with Algorithm 11
    if hub_tags:
        print(f"\n=== Remeshing Hub (Algorithm 11) ===")
        try:
            hub_quads = remesh_surface_with_quads(hub_tags, element_size=0.5)
            print(f"Generated {len(hub_quads)} quads on hub")
        except Exception as e:
            print(f"Hub remeshing failed: {e}")
    
    if shroud_tags:
        print(f"\n=== Remeshing Shroud (Algorithm 11) ===")
        try:
            shroud_quads = remesh_surface_with_quads(shroud_tags, element_size=0.5)
            print(f"Generated {len(shroud_quads)} quads on shroud")
        except Exception as e:
            print(f"Shroud remeshing failed: {e}")
    
    # Export
    output_path = "/root/repos/block_structured_meshing/T1_9_coarse_blocks.msh"
    print(f"\nExporting to {output_path}...")
    gmsh.write(output_path)
    
    gmsh.finalize()
    print("\nDone!")

if __name__ == '__main__':
    main()
