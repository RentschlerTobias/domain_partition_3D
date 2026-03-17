

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
import plotly.graph_objects as go
import webbrowser
import os


def streamlines_to_html(streamlines, output_path="streamlines.html"):
    fig = go.Figure()

    for streamline in streamlines:
        x, y, z = streamline[:, 0], streamline[:, 1], streamline[:, 2]

        # zufällige Farbe für jede Stromlinie
        color = np.random.randint(0, 255, size=3)
        color = f'rgb({color[0]},{color[1]},{color[2]})'

        fig.add_trace(go.Scatter3d(
            x=x, y=y, z=z,
            mode='lines',
            line=dict(color=color, width=4)
        ))

    fig.update_layout(
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z'
        ),
        width=800,
        height=800,
        title="3D Stromlinien"
    )

    # Speichern als HTML
    fig.write_html(output_path)
    print(f"Plot gespeichert als: {output_path}")


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


streamline_extractor = StreamlineExtractor()
extracted_data = streamline_extractor.extract_mesh_from_msh(
    path_msh_file="./msh_files/remeshing.msh")
streamline_intersection = StreamlineIntersections()
streamlines = streamline_intersection.get_intersections(
    extracted_data['streamlines'])

for surface_tag in streamlines.keys():
    streamlines_to_html(
        streamlines[surface_tag]['streamlines_split'], output_path="./figures/stromlinien.html")


type(extracted_data['streamlines'][2])
streamlines_preProcessed = extracted_data['streamlines'][2]['points']
streamlines_postProcessed = streamlines[surface_tag]['streamlines_split']


len(streamlines_preProcessed)
len(streamlines_postProcessed)

streamlines_to_html(streamlines_preProcessed,
                    output_path="./figures/streamlines_preProcessed.html")
