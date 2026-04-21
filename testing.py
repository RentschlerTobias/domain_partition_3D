import gmsh
from collections import defaultdict
from mesh_generator2D import *
from naca_airfoil import *

from streamline_extractor import StreamlineExtractor
from streamlines_to_msh import export_streamlines_to_msh
from streamline_splitter import split_streamlines
from streamline_intersections import *
from plotting_tools import *
import numpy as np


from msh_extractor import MshExtractor
from streamline_extractor import StreamlineExtractor


airfoil = NACA_airfoil()
random_lc = 0.04 + 0.02 * np.random.rand()
mesh_gen = MeshGenerator(airfoil, quadMesh=False, lc=random_lc)
mesh = mesh_gen.mesh

vertices = mesh.x
faces = mesh.faces

plot_mesh(vertices, faces, output_file="./quad_mesh.png")

vertices, surfaces = get_mesh_of_airfoil(random_lc, airfoil)

faces = surfaces[1].T


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


export_streamlines_to_msh(streamlines, output_path='./streamlines2D.msh')
