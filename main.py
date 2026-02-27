

from remesh import mesh_stl
from streamline_extraction import StrealimeExtractor
from streamline_intersections import *

# stl to msh
# triangulated surface to quasi-structured surface
# triangulated volume to quasi-structured surface is possible, gmsh automatically seperates surface into multiple surfaces based on params in remsh.py file
path_st = "./stl_files/Case09_post/LV_outer.stl"
path_msh = "LV_outer.msh"
name_streamlines_msh = 'LV_outer_streamlines_surface_'
mesh_stl(
    path_to_stl=path_st,
    output_path=path_msh,
    element_size=4
)

streamline_extractor = StrealimeExtractor()
extracted_data = streamline_extractor.extract_mesh_from_msh(
    path_msh_file=path_msh)

for surface_tag in extracted_data['streamlines'].keys():

    current_streamlines = extracted_data['streamlines'][surface_tag]['points']
    path = f'{name_streamlines_msh}_{surface_tag}.msh'
    streamline_extractor.export_streamlines_to_msh(
        streamlines=current_streamlines, output_path=path)
    print(f'tag is: {surface_tag}')
