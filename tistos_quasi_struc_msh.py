
from tistos.build import build
import sys
sys.path.insert(0, '/mnt/opt.net/src/dtOO/ThirdParty~opensuse-leap~15.6/lib64')
import gmsh
import numpy as np

#
class QuasiStructuredDTOO:
    def __init__(self, element_size=4,output_path="./dtOO_remeshed.msh"):

        self.vertices = None
        self.surfaces = {}  # Dictionary: {surface_tag: faces_array}

        self._get_gmodel(element_size,output_path) 
        # self._extract_data(output_path)

    def _get_gmodel(self, element_size,output_path):

        gmsh.initialize()
        bC, cV, aF, aG, bV, dC, dP = build()
        # bv_n = 2
        # gmodel = bV[bv_n]
        gmsh.model.set_current('ru_gridGmsh')

        gmsh.model.mesh.generate(2)
        # surfaces = gmsh.model.getEntities(2)
        # tags = [s[1] for s in surfaces]
        # gmsh.model.geo.addVolume([gmsh.model.geo.addSurfaceLoop(tags)])
        # gmsh.model.geo.synchronize()
        #
        field_id = gmsh.model.mesh.field.add("MathEval")
        # gmsh.model.mesh.field.setString(
        #     field_id, "F", f"{element_size} + 0.1 * x")
        gmsh.model.mesh.field.setString(field_id, "F", str(element_size))
        gmsh.model.mesh.field.setAsBackgroundMesh(field_id)
        # Set Algorithm 11: Quasi-structured Quad
        gmsh.option.setNumber("Mesh.Algorithm", 11)
        gmsh.model.mesh.generate(2)

        gmsh.write(output_path)
        gmsh.finalize()
    
    def _extract_data(self, path):
        gmsh.initialize()
        gmsh.open(path)

        # 1. Nodes (Shared globally)
        node_tags, coords, _ = gmsh.model.mesh.getNodes()
        self.vertices = coords.reshape(-1, 3)

        tag_to_idx = np.full(int(node_tags.max()) + 1, -1)
        tag_to_idx[node_tags.astype(int)] = np.arange(len(node_tags))

        # 2. Extract faces per Surface
        for _, s_tag in gmsh.model.getEntities(2):
            _, q_tags = gmsh.model.mesh.getElementsByType(3, s_tag)
            if q_tags.size:
                # Store faces for this specific surface tag
                self.surfaces[s_tag] = tag_to_idx[q_tags.astype(
                    int)].reshape(-1, 4)

        gmsh.finalize()


if __name__ == '__main__':
    output_path = './tistos_test.msh'
    extractor = QuasiStructuredDTOO(
        element_size=1.0,
        output_path=output_path,
    )

