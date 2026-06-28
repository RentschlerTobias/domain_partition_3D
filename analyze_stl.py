import gmsh
import numpy as np

OUT_DIR = "/root/repos/block_structured_meshing"

def analyze_stl(stl_file, name):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    
    gmsh.merge(stl_file)
    
    # get nodes
    tags, coords, _ = gmsh.model.mesh.getNodes()
    coords = np.array(coords).reshape(-1, 3)
    
    # get triangles
    tri_nodes = []
    elem_types = gmsh.model.mesh.getElementTypes()
    for et in elem_types:
        if et == 2:  # 3-node triangle
            etags, enodes = gmsh.model.mesh.getElementsByType(et)
            tri_nodes = np.array(enodes).reshape(-1, 3)
            break
    
    print(f"[{name}] {len(tri_nodes)} triangles, {len(tags)} nodes")
    print(f"[{name}] coord range: x={coords[:,0].min():.3f}-{coords[:,0].max():.3f}, y={coords[:,1].min():.3f}-{coords[:,1].max():.3f}, z={coords[:,2].min():.3f}-{coords[:,2].max():.3f}")
    
    # compute boundary edges
    edge_map = {}
    for tri in tri_nodes:
        for i in range(3):
            n1, n2 = sorted([tri[i], tri[(i+1)%3]])
            edge = (n1, n2)
            if edge in edge_map:
                edge_map[edge] += 1
            else:
                edge_map[edge] = 1
    
    boundary_edges = [e for e, count in edge_map.items() if count == 1]
    print(f"[{name}] boundary edges: {len(boundary_edges)}")
    
    # find connected components of boundary edges
    # build graph
    from collections import defaultdict
    graph = defaultdict(set)
    for n1, n2 in boundary_edges:
        graph[n1].add(n2)
        graph[n2].add(n1)
    
    visited = set()
    loops = []
    
    for start in graph:
        if start in visited:
            continue
        # BFS/DFS to find connected component
        loop = []
        stack = [start]
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            loop.append(node)
            for nb in graph[node]:
                if nb not in visited:
                    stack.append(nb)
        
        # order loop
        if len(loop) > 2:
            ordered = [loop[0]]
            current = loop[0]
            while True:
                neighbors = [n for n in graph[current] if n not in ordered or n == loop[0]]
                if not neighbors:
                    break
                next_n = neighbors[0]
                if next_n == ordered[0] and len(ordered) > 2:
                    break
                if next_n in ordered and next_n != ordered[0]:
                    break
                ordered.append(next_n)
                current = next_n
            
            loops.append(ordered)
    
    print(f"[{name}] boundary loops: {len(loops)}")
    for i, loop in enumerate(loops):
        print(f"  loop {i}: {len(loop)} nodes")
        # compute center
        loop_coords = coords[[np.where(tags == n)[0][0] for n in loop]]
        center = loop_coords.mean(axis=0)
        print(f"    center: {center}")
        # compute radius
        radii = np.linalg.norm(loop_coords - center, axis=1)
        print(f"    radius: min={radii.min():.3f}, max={radii.max():.3f}, avg={radii.mean():.3f}")
    
    gmsh.finalize()

if __name__ == "__main__":
    analyze_stl(f"{OUT_DIR}/T1_9_hub_raw.stl", "hub")
    analyze_stl(f"{OUT_DIR}/T1_9_shroud_raw.stl", "shroud")
