

import torch
from remesh import mesh_stl

from streamline_extractor import StreamlineExtractor
# from streamline_extraction import StreamlineExtractor
from streamline_intersections import StreamlineIntersections
from remeshing import MshExtractor
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
import webbrowser
import numpy as np


def plt_streamlines(streamlines, output_path='./figures/streamlines.png'):

    figsize = (5, 5)
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(projection='3d')

    for streamline in streamlines:
        x = streamline[:, 0]
        y = streamline[:, 1]
        z = streamline[:, 2]

        color = np.random.rand(3)
        ax.plot(x, y, z, color=color)
    plt.savefig(output_path)
#
#
# streamline_extractor = StreamlineExtractor()
# extracted_data = streamline_extractor.extract_mesh_from_msh(
#     path_msh_file="./msh_files/remeshing.msh")
# streamline_intersection = StreamlineIntersections()
# streamlines = streamline_intersection.get_intersections(
#     extracted_data['streamlines'])
#
# for surface_tag in streamlines.keys():
#     streamlines_to_html(
#         streamlines[surface_tag]['streamlines_split'], output_path="./figures/stromlinien.html")
#
#
# type(extracted_data['streamlines'][2])
# streamlines_preProcessed = extracted_data['streamlines'][2]['points']
# streamlines_postProcessed = streamlines[surface_tag]['streamlines_split']
#
#
# len(streamlines_preProcessed)
# len(streamlines_postProcessed)
#
# streamlines_to_html(streamlines_preProcessed,
#                     output_path="./figures/streamlines_preProcessed.html")
