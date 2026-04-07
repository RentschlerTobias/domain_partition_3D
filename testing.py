from msh_extractor import MshExtractor
from streamline_extractor import StreamlineExtractor
from streamlines_to_msh import export_streamlines_to_msh
from streamline_splitter import split_streamlines
from plotting_tools import streamlines_to_html, block_structure_to_html
from streamline_intersections_v2 import get_block_structure_from_streamlines, detect_quad_faces
import numpy as np
# input_path = './stl_files/Case09_post/LV_outer.stl'
input_path = './stl_files/Case09_post/LV_inner.stl'
# input_path = './stl_files/Case199_pre_all.stl'
output_path = './remeshed_quads.msh'

extractor = MshExtractor(
    input_path=input_path,          # Path to your STL
    remesh=True,                    # Enable GMSH remeshing
    output_path=output_path,  # Where to save the intermediate MSH
    element_size=2.5,               # Target edge length
    angle=50,                        # Angle for surface classification
    curve_angle=180
)

# Access the extracted data
vertices = extractor.vertices  # Numpy array (N, 3)
surfaces = extractor.surfaces
streamlines = {}

for surface_tag in surfaces.keys():
    faces = surfaces[surface_tag]
    streamline_extractor = StreamlineExtractor(vertices, faces)
    streamlines[surface_tag] = streamline_extractor.get_streamlines()

output_path = './streamlines_1.msh'
export_streamlines_to_msh(streamlines, output_path=output_path)

new_streamlines = split_streamlines(streamlines)

block_structure = get_block_structure_from_streamlines(new_streamlines)

vertices = block_structure[2]['vertices']
edges = block_structure[2]['edges']

quad_faces = detect_quad_faces(vertices, edges)
len(quad_faces)


type(streamlines)
len(streamlines[2])
len(new_streamlines[2])

data = {}


for surface_tag in new_streamlines.keys():
    streamlines_of_current_surfaces = new_streamlines[surface_tag]
    for current_streamline in streamlines_of_current_surfaces:
        start_point = current_streamline[0, :]
        end_point = current_streamline[-1, :]
        print(f'streamline: \n {current_streamline}\n\n start: {
              start_point} \n end: {end_point}\n')
        key = (start_point, end_point)
        print(key)
        break


type(block_structure)
block_structure.keys()
type(block_structure[2])
block_structure[2].keys()
type(block_structure)

output_path = f"./figures/blocking.html"
block_structure_to_html(block_structure, output_path=output_path)

for surface_tag in block_structure.keys():
    output_path = f"blocking_{surface_tag}.html"
    s = new_streamlines[surface_tag]
    # streamlines_to_html(s, output_path=output_path)
    block_structure_to_html(s, output_path=output_path)
