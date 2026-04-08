
from msh_extractor import MshExtractor
from streamline_extractor import StreamlineExtractor
from streamlines_to_msh import export_streamlines_to_msh
from streamline_splitter import split_streamlines
from streamline_intersections import *
from plotting_tools import *
import numpy as np


def main():
    input_path = './stl_files/Case09_post/LV_outer.stl'
    # input_path = './stl_files/Case09_post/LV_inner.stl'
    # input_path = './stl_files/Case199_post/LV_outer.stl'
    # input_path = './stl_files/Case199_pre_all.stl'
    output_path = './remeshed_quads.msh'

    extractor = MshExtractor(
        input_path=input_path,          # Path to your STL
        remesh=True,                    # Enable GMSH remeshing
        output_path=output_path,        # Where to save the intermediate MSH
        element_size=2.5,               # Target edge length
        angle=50,                        # Angle for surface classification
        curve_angle=180
    )

    vertices = extractor.vertices  # Numpy array (N, 3)
    surfaces = extractor.surfaces

    streamlines = {}

    for surface_tag in surfaces.keys():
        faces = surfaces[surface_tag]
        streamline_extractor = StreamlineExtractor(vertices, faces)
        streamlines[surface_tag] = streamline_extractor.get_streamlines()

    splitted_streamlines = split_streamlines(streamlines)
    block_structure = get_block_structure_from_streamlines(
        splitted_streamlines)

    for surface_tag in block_structure.keys():
        vertices = block_structure[surface_tag]['vertices']
        edges = block_structure[surface_tag]['edges']

        edge_to_streamline = block_structure[surface_tag]['edge_to_streamline']
        block_structure[surface_tag]['faces'] = detect_quad_faces(
            vertices, edges)

    output_path_figure_surfaces = f"./html_files/Case09_post_LV_outer.html"
    output_path_blocking_structure = f"./figures/blocking_structure_linear.html"

    # Optional for visualisation purpose
    # html graphics only for small geometries recommendated and small n_u,. n_v values

    # export_streamlines_to_msh(streamlines, output_path='./streamlines.msh')
    #
    faces_to_html(block_structure,
                  output_path=output_path_figure_surfaces, n_u=20, n_v=20)
    #
    # block_structure_to_html(
    #     block_structure, output_path=output_path_blocking_structure)
    #
    # surfaces_to_msh(block_structure, n_u=None, n_v=None)

    # Visualise the surfaces induvidually
    # for surface_tag in block_structure.keys():
    #     output_path = f"blocking_{surface_tag}.html"
    #     s = new_streamlines[surface_tag]
    #     # streamlines_to_html(s, output_path=output_path)
    #     block_structure_to_html(s, output_path=output_path)


if __name__ == '__main__':
    main()
