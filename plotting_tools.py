import plotly.graph_objects as go
import webbrowser
import os
import numpy as np


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


def block_structure_to_html(graph_data, output_path="blocking.html"):
    """
    Input: dict -> keys are surface tags, values contain 'vertices', 'edges', 'edge_to_streamline', 'streamlines'

    Output: None (saves HTML file)

    Plots 3D visualization showing:
        - Streamlines as colored curves
        - Vertices as points at start/end locations
        - Edges as lines connecting vertices with matching streamline colors
    """

    fig = go.Figure()

    # Assign unique color to each streamline for consistency between streamline and edge
    streamline_colors = {}  # (surface_tag, streamline_idx) -> color string

    for surface_tag, data in graph_data.items():
        vertices = data['vertices']
        edges = data['edges']
        streamlines = data['streamlines']

        # Generate unique colors for each streamline in this surface
        for sl_idx in range(len(streamlines)):
            color = np.random.randint(0, 255, size=3)
            color_str = f'rgb({color[0]},{color[1]},{color[2]})'
            streamline_colors[(surface_tag, sl_idx)] = color_str

        # Plot streamlines with assigned colors
        for sl_idx, streamline in enumerate(streamlines):
            x, y, z = streamline[:, 0], streamline[:, 1], streamline[:, 2]
            color_str = streamline_colors[(surface_tag, sl_idx)]
            fig.add_trace(go.Scatter3d(
                x=x, y=y, z=z,
                mode='lines',
                line=dict(color=color_str, width=4),
                name=f'{surface_tag} streamline {
                    sl_idx}' if sl_idx == 0 else None
            ))

        # Plot edges with matching colors (use color of first associated streamline)
        for edge_idx, edge in enumerate(edges):
            start_vertex = vertices[edge[0]]
            end_vertex = vertices[edge[1]]

            # Get the first streamline index associated with this edge to determine color
            sl_indices = data['edge_to_streamline'].get(edge_idx, [])
            if sl_indices:
                first_sl_idx = sl_indices[0]
                edge_color = streamline_colors.get(
                    (surface_tag, first_sl_idx), 'gray')
            else:
                edge_color = 'gray'

            fig.add_trace(go.Scatter3d(
                x=[start_vertex[0], end_vertex[0]],
                y=[start_vertex[1], end_vertex[1]],
                z=[start_vertex[2], end_vertex[2]],
                mode='lines',
                line=dict(color=edge_color, width=4),
                name=f'{surface_tag} edge' if edge_idx == 0 else None
            ))

        # Plot vertices as points
        fig.add_trace(go.Scatter3d(
            x=vertices[:, 0], y=vertices[:, 1], z=vertices[:, 2],
            mode='markers',
            marker=dict(size=8, color='black'),
            name=f'{surface_tag} vertices' if surface_tag == list(graph_data.keys())[
                0] else None
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
