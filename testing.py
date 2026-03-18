from msh_extractor import MshExtractor
from streamline_extractor import StreamlineExtractor
from streamlines_to_msh import export_streamlines_to_msh

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
