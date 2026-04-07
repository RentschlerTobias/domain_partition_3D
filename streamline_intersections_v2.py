import networkx as nx
import numpy as np
from collections import defaultdict


def get_block_structure_from_streamlines(streamlines_of_all_surfaces):
    """
    Input: dict -> keys are surface tags, values are lists of streamlines (numpy arrays of shape [n_points, n_dims])

    Output: dict -> each key is a surface tag with value containing:
        - 'vertices': numpy.ndarray of shape [n_vertices, n_dimensions], unique start/end points for this surface
        - 'edges': numpy.ndarray of shape [n_edges, 2], indices into vertices array for each streamline connection
        - 'edge_to_streamline': dict mapping edge index to list of streamline_index integers

    Builds a graph structure per surface where:
        - Each streamline becomes an edge between two vertices (start and end points)
        - Duplicate start/end points share the same vertex index within that surface
        - Multiple streamlines can connect the same pair of vertices
    """

    result = {}

    for surface_tag, streamlines in streamlines_of_all_surfaces.items():
        # Unique point storage with hashable keys
        point_to_vertex = {}  # Maps tuple(point) -> vertex_index
        vertices_list = []     # List of unique points (as arrays)

        # Edge data structures
        edges_list = []        # Each entry: [start_vertex_idx, end_vertex_idx]
        # Maps edge_idx -> list of streamline indices
        edge_to_streamline = defaultdict(list)

        vertex_counter = 0
        edge_counter = 0

        # Iterate through all streamlines for this surface
        for streamline_idx, current_streamline in enumerate(streamlines):
            start_point = current_streamline[0, :]
            end_point = current_streamline[-1, :]

            # Create hashable keys from points (round to avoid floating point issues)
            start_key = tuple(np.round(start_point, decimals=9))
            end_key = tuple(np.round(end_point, decimals=9))

            # Add start vertex if not exists
            if start_key not in point_to_vertex:
                point_to_vertex[start_key] = vertex_counter
                vertices_list.append(start_point)
                vertex_counter += 1

            # Add end vertex if not exists
            if end_key not in point_to_vertex:
                point_to_vertex[end_key] = vertex_counter
                vertices_list.append(end_point)
                vertex_counter += 1

            # Store edge with vertex indices
            start_idx = point_to_vertex[start_key]
            end_idx = point_to_vertex[end_key]

            edges_list.append([start_idx, end_idx])
            edge_to_streamline[edge_counter].append(streamline_idx)
            edge_counter += 1

        # Convert lists to numpy arrays
        vertices = np.array(
            vertices_list) if vertices_list else np.empty((0, 0))
        edges = np.array(edges_list) if edges_list else np.empty(
            (0, 2)).astype(int)

        result[surface_tag] = {
            'vertices': vertices,
            'edges': edges,
            'edge_to_streamline': dict(edge_to_streamline),
            'streamlines': streamlines
        }

    return result


def detect_quad_faces(vertices: np.ndarray, edges: np.ndarray) -> list[list[int]]:
    """
    Detects quadrilateral faces in a mesh by finding cycles of length 4.

    Input: 
        vertices (numpy.ndarray): Array of shape [n_vertices, n_dimensions].
                                 Unique start/end points for the surface.
        edges (numpy.ndarray): Array of shape [n_edges, 2].
                               Indices into vertices array for connected vertices.

    Output: 
        List[List[int]]: A list where each element is a list of 4 vertex indices forming a quad face.
                         Returns an empty list if no quads are found.
    """
    # 1. Construct the graph using NetworkX
    G = nx.Graph()

    # Add edges directly; nodes will be created automatically based on edge indices
    # Using np.nditer for efficient iteration over the edges array
    G.add_edges_from(edges)

    # 2. Find the cycle basis of the graph
    # The cycle_basis returns a list of cycles (lists of node indices).
    # These form a basis for the cycle space of the graph.
    try:
        cycles = nx.cycle_basis(G)
    except Exception:
        return []

    # 3. Filter for cycles with exactly length 4 (Quadrilaterals)
    quad_faces = [cycle for cycle in cycles if len(cycle) == 4]

    return quad_faces
