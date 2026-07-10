#!/usr/bin/env python3
"""
Map dtOO block corners to actual hex mesh regions.
Compare dtOO block definitions with hex element regions.
"""

import numpy as np
import json
from collections import defaultdict

def parse_msh(filename):
    nodes = {}
    hex_elements = []
    
    with open(filename, 'r') as f:
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

def load_dtOO_blocks(filename):
    with open(filename, 'r') as f:
        return json.load(f)

def find_nearest_hex_to_corner(dtOO_corner, hex_elements, nodes, max_dist=0.1):
    """Find hex elements near a dtOO block corner."""
    nearby_elems = []
    
    for elem in hex_elements:
        elem_nodes = [nodes[n] for n in elem['nodes']]
        elem_center = np.mean(elem_nodes, axis=0)
        dist = np.linalg.norm(elem_center - dtOO_corner)
        
        if dist < max_dist:
            nearby_elems.append({
                'elem_tag': elem['tag'],
                'geom_tag': elem['geom_tag'],
                'distance': float(dist),
                'center': elem_center.tolist()
            })
    
    return sorted(nearby_elems, key=lambda x: x['distance'])

def map_dtOO_to_regions(dtOO_blocks, hex_elements, nodes):
    """Map dtOO blocks to hex regions."""
    
    # Group hex by region
    hex_by_region = defaultdict(list)
    for elem in hex_elements:
        hex_by_region[elem['geom_tag']].append(elem)
    
    mapping = {}
    
    for block_name, block_data in dtOO_blocks.items():
        print(f"\nBlock: {block_name}")
        
        if isinstance(block_data, dict) and 'physical' in block_data:
            corners = block_data['physical']
        else:
            corners = block_data
            
        print(f"  Corners: {len(corners)}")
        
        # Check which region contains most corners
        region_votes = defaultdict(int)
        
        for i, corner in enumerate(corners):
            corner = np.array(corner, dtype=float)
            
            # Find nearest hex element
            best_dist = float('inf')
            best_region = None
            
            for geom_tag, elems in hex_by_region.items():
                for elem in elems:
                    elem_nodes = [nodes[n] for n in elem['nodes']]
                    elem_center = np.mean(elem_nodes, axis=0)
                    dist = np.linalg.norm(elem_center - corner)
                    
                    if dist < best_dist:
                        best_dist = dist
                        best_region = geom_tag
            
            if best_region is not None:
                region_votes[best_region] += 1
            
            print(f"    Corner {i}: nearest region={best_region}, dist={best_dist:.4f}")
        
        # Determine best matching region
        if region_votes:
            best_region = max(region_votes, key=region_votes.get)
            confidence = region_votes[best_region] / len(corners)
            print(f"  -> Best match: Region {best_region} ({region_votes[best_region]}/{len(corners)} corners, {confidence:.1%} confidence)")
            
            mapping[block_name] = {
                'region': best_region,
                'confidence': confidence,
                'votes': dict(region_votes)
            }
    
    return mapping

def main():
    # Load dtOO blocks (if they exist)
    dtOO_files = [
        'dtOO_blocks_mapped.json',
        'block_instances_28.json'
    ]
    
    dtOO_blocks = None
    for f in dtOO_files:
        try:
            with open(f, 'r') as file:
                data = json.load(file)
                if 'blocks' in data:
                    dtOO_blocks = data['blocks']
                else:
                    dtOO_blocks = data
                print(f"Loaded dtOO blocks from {f}")
                break
        except:
            continue
    
    if dtOO_blocks is None:
        print("No dtOO block files found. Please run parse_dtOO_blocks.py first.")
        return
    
    # Parse mesh
    print("\nParsing mesh...")
    nodes, hex_elements = parse_msh('T1_9/T1_9_ru_gridGmsh.msh')
    print(f"Loaded {len(nodes)} nodes, {len(hex_elements)} hex elements")
    
    # Map dtOO blocks to regions
    print("\n" + "="*60)
    print("MAPPING dtOO BLOCKS TO HEX REGIONS")
    print("="*60)
    
    mapping = map_dtOO_to_regions(dtOO_blocks, hex_elements, nodes)
    
    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    for block_name, info in mapping.items():
        print(f"{block_name} -> Region {info['region']} (confidence: {info['confidence']:.1%})")
    
    # Save mapping
    with open('dtoo_to_region_mapping.json', 'w') as f:
        json.dump(mapping, f, indent=2)
    print("\nMapping saved to dtoo_to_region_mapping.json")

if __name__ == '__main__':
    main()
