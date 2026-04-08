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

        vertex_pair_to_streamline = {}
        for edge_idx, sl_indices in edge_to_streamline.items():
            va, vb = int(edges[edge_idx][0]), int(edges[edge_idx][1])
            sl = streamlines[sl_indices[0]]
            vertex_pair_to_streamline[(va, vb)] = sl
            vertex_pair_to_streamline[(vb, va)] = sl[::-1]

        result[surface_tag] = {
            'vertices': vertices,
            'edges': edges,
            'edge_to_streamline': dict(vertex_pair_to_streamline),
            'streamlines': streamlines
        }

    return result


def detect_quad_faces(vertices: np.ndarray, edges: np.ndarray) -> list[list[int]]:
    """
    Detects quadrilateral faces in a mesh by finding all simple cycles of length 4.

    Uses diagonal-based detection: for each pair of non-adjacent vertices that share 
    exactly 2 common neighbors, those 4 vertices form a quad face.

    Input: 
        vertices (numpy.ndarray): Array of shape [n_vertices, n_dimensions].
                                  Unique start/end points for the surface.
        edges (np.ndarray): Array of shape [n_edges, 2].
                            Indices into vertices array for connected vertices.

    Output: 
        List[List[int]]: A list where each element is a list of 4 vertex indices forming a quad face.
                         Returns an empty list if no quads are found.
    """
    n_vertices = len(vertices)

    # Build adjacency matrix and neighbor lists (vectorized)
    adj_matrix = np.zeros((n_vertices, n_vertices), dtype=np.bool_)
    adj_matrix[edges[:, 0], edges[:, 1]] = True
    adj_matrix[edges[:, 1], edges[:, 0]] = True

    # Build neighbor list for iteration
    neighbors = defaultdict(list)
    for u, v in edges:
        neighbors[u].append(v)
        neighbors[v].append(u)

    quad_faces = []
    seen_quads = set()

    # For each pair of non-adjacent vertices (potential diagonal), check if they share exactly 2 common neighbors
    # These form a quadrilateral face
    for i in range(n_vertices):
        for j in range(i + 1, n_vertices):
            # Skip if directly connected (would be edge, not diagonal)
            if adj_matrix[i, j]:
                continue

            # Find common neighbors of i and j
            common = np.where((adj_matrix[i] & adj_matrix[j]))[0]

            # Exactly 2 common neighbors means a quad face exists with i and j as diagonal
            if len(common) == 2:
                k, l = common[0], common[1]

                # Verify this forms a valid cycle (all edges exist)
                if adj_matrix[i, k] and adj_matrix[k, j] and adj_matrix[j, l] and adj_matrix[l, i]:
                    quad = tuple(sorted([i, j, k, l]))
                    if quad not in seen_quads:
                        seen_quads.add(quad)
                        # Order vertices to form a proper cycle around the face
                        quad_faces.append([i, k, j, l])

    quads = np.array(quad_faces).T

    return quads
