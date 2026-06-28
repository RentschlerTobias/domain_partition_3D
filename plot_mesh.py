import gmsh
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

OUT_DIR = "/root/repos/block_structured_meshing"

def plot_msh(msh_file, name, png_file):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    
    gmsh.open(msh_file)
    
    # get nodes
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    node_coords = np.array(node_coords).reshape(-1, 3)
    node_map = {int(tag): i for i, tag in enumerate(node_tags)}
    
    # get quads
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(2)
    quads = []
    for etype, etags, enodes in zip(elem_types, elem_tags, elem_node_tags):
        if etype == 3:  # 4-node quad
            n_node_per_elem = 4
            for i in range(len(etags)):
                n = enodes[i*n_node_per_elem:(i+1)*n_node_per_elem]
                quads.append([int(x) for x in n])
    
    print(f"[{name}] {len(quads)} quads, {len(node_tags)} nodes")
    print(f"[{name}] node coords: min={node_coords.min(axis=0)}, max={node_coords.max(axis=0)}")
    
    # plot
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # plot quads
    for quad in quads:
        idx = [node_map[x] for x in quad]
        pts = node_coords[idx]
        # close the quad
        pts_closed = np.vstack([pts, pts[0]])
        ax.plot(pts_closed[:, 0], pts_closed[:, 1], pts_closed[:, 2], 'b-', alpha=0.5, linewidth=0.5)
        
        # fill
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        verts = [pts]
        poly3d = Poly3DCollection(verts, alpha=0.1, facecolor='cyan', edgecolor='none')
        ax.add_collection3d(poly3d)
    
    # plot nodes
    ax.scatter(node_coords[:, 0], node_coords[:, 1], node_coords[:, 2], c='red', s=5)
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'{name}: {len(quads)} quads')
    
    plt.savefig(png_file, dpi=150)
    print(f"[{name}] saved plot to {png_file}")
    
    gmsh.finalize()

if __name__ == "__main__":
    plot_msh(f"{OUT_DIR}/T1_9_hub_quad.msh", "hub", f"{OUT_DIR}/T1_9_hub_quad.png")
    plot_msh(f"{OUT_DIR}/T1_9_shroud_quad.msh", "shroud", f"{OUT_DIR}/T1_9_shroud_quad.png")
