import gmsh
import numpy as np
import openmesh as om


class MshLoader:
    def __init__(self, path):
        self.path = path
        self.vertices = None
        self.faces = None
        self._load()

    def _load(self):
        """Internal gmsh extractor."""
        gmsh.initialize()
        gmsh.open(self.path)

        # Extract nodes
        node_tags, coords, _ = gmsh.model.mesh.getNodes()
        self.vertices = np.reshape(coords, (-1, 3))

        # Tag to index mapping
        max_tag = int(np.max(node_tags))
        tag_to_idx = np.full(max_tag + 1, -1, dtype=np.int32)
        tag_to_idx[node_tags.astype(int)] = np.arange(len(node_tags))

        # Extract Quads from all 2D surfaces
        all_faces = []
        surfaces = gmsh.model.getEntities(dim=2)
        for _, s_tag in surfaces:
            _, q_node_tags = gmsh.model.mesh.getElementsByType(3, tag=s_tag)
            if len(q_node_tags) > 0:
                indices = tag_to_idx[q_node_tags.astype(int)]
                all_faces.append(indices.reshape(-1, 4))

        if all_faces:
            self.faces = np.vstack(all_faces)

        gmsh.finalize()

    def get_mesh_data(self):
        """Returns vertices and face connectivity."""
        return self.vertices, self.faces
