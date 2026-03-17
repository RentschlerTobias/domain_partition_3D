import gmsh
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.colors import to_rgba


class StreamlineExtractor:

    def get_singularities(self, vertices, faces):
        # 1. Erzeuge alle Kanten (wie bisher)
        edges = np.concatenate([
            faces[[0, 1], :],
            faces[[1, 2], :],
            faces[[2, 3], :],
            faces[[3, 0], :]
        ], axis=1)

        # 2. Eindeutige Kanten finden (Wichtig: sortieren für ungerichtete Erkennung)
        sorted_edges = np.sort(edges, axis=0)
        unique_edges, edge_counts = np.unique(
            sorted_edges, axis=1, return_counts=True)

        # 3. Boundary finden
        boundary_edges = unique_edges[:, edge_counts == 1]
        boundary_nodes = np.unique(boundary_edges)

        # 4. ECHTE Valenz berechnen (basierend auf unique_edges)
        real_valences = np.bincount(
            unique_edges.flatten(), minlength=vertices.shape[0])

        # 5. Singularitäts-Maske
        # Wir suchen Punkte, die:
        # - Nicht am Rand liegen
        # - Eine Valenz ungleich 4 haben
        # - Überhaupt Teil der Oberflaeche sind (Valenz > 0)
        singular_mask = (real_valences != 4) & (real_valences > 0)

        # Boundary-Knoten explizit ausschließen
        # (denn ein Randknoten mit Valenz 3 ist keine Singularität)
        is_boundary = np.zeros(vertices.shape[0], dtype=bool)
        is_boundary[boundary_nodes] = True
        singular_mask &= ~is_boundary

        # 6. Ergebnisse extrahieren
        singular_node_indices = np.where(singular_mask)[0]
        singular_nodes = vertices[singular_mask, :]
        singular_valences = real_valences[singular_mask]

        return singular_nodes, singular_node_indices, singular_valences, boundary_nodes, unique_edges

    def build_adjacency(self, vertices, faces, unique_edges):
        """Baut Adjazenzliste: node -> Liste von Nachbarknoten"""
        from collections import defaultdict
        adj = defaultdict(list)
        for i in range(unique_edges.shape[1]):
            u, v = unique_edges[0, i], unique_edges[1, i]
            adj[u].append(v)
            adj[v].append(u)
        return adj

    def trace_separatrix(self, start_node, start_direction_node, adj,
                         singular_set, boundary_set, vertices, max_steps=10000):
        """
        Verfolgt eine Separatrix vom start_node aus in Richtung start_direction_node.

        Wichtig: Bei Quad-Meshes an jedem regulären Knoten (Valenz 4) 
        'geradeaus' weitergehen = gegenüberliegenden Nachbarn nehmen.
        """
        path = [start_node, start_direction_node]
        prev = start_node
        curr = start_direction_node

        for _ in range(max_steps):
            # Stopp-Bedingungen
            if curr in singular_set or curr in boundary_set:
                break

            neighbors = adj[curr]

            # "Geradeaus" = nicht zurück (prev) und nicht 90°-Abbieger
            # Bei Valenz 4: der gegenüberliegende Nachbar ist derjenige,
            # der am meisten in Richtung (curr - prev) liegt
            incoming_dir = vertices[curr] - vertices[prev]
            incoming_dir /= (np.linalg.norm(incoming_dir) + 1e-12)

            best_node = None
            best_dot = -np.inf

            for nb in neighbors:
                if nb == prev:
                    continue
                outgoing_dir = vertices[nb] - vertices[curr]
                outgoing_dir /= (np.linalg.norm(outgoing_dir) + 1e-12)
                dot = np.dot(incoming_dir, outgoing_dir)
                if dot > best_dot:
                    best_dot = dot
                    best_node = nb

            if best_node is None:
                break

            path.append(best_node)
            prev = curr
            curr = best_node

        return path, curr  # path + Endknoten

    def trace_all_separatrices(self, vertices, faces, singular_node_indices,
                               boundary_nodes, unique_edges):
        adj = self.build_adjacency(vertices, faces, unique_edges)
        singular_set = set(singular_node_indices)
        boundary_set = set(boundary_nodes)

        all_separatrices = []

        for s in singular_node_indices:
            neighbors = adj[s]
            # Starte in jede Richtung (jeder Nachbar = eine Separatrix)
            for nb in neighbors:
                path, end_node = self.trace_separatrix(
                    s, nb, adj, singular_set, boundary_set, vertices
                )
                end_type = "singularity" if end_node in singular_set else "boundary"
                all_separatrices.append({
                    "path": path,
                    "start": s,
                    "end": end_node,
                    "end_type": end_type
                })

        return all_separatrices

    def filter_unique_streamlines(self, streamlines, decimals=9):
        """
        Filters duplicate streamlines, handling both directions and closed loops.
        """
        unique_streamlines = []
        seen_hashes = set()

        for sl in streamlines:
            if len(sl) < 2:
                continue

            # Round to handle floating point noise from integration
            sl_norm = np.round(sl, decimals=decimals)

            # Create reversed version
            sl_rev = sl_norm[::-1]

            # Canonical form: lexicographical comparison of the entire point sequence
            # Ensures [A, B, C, A] and [A, C, B, A] are treated identically
            if np.lexsort(sl_rev.T)[0] < np.lexsort(sl_norm.T)[0]:
                target = sl_rev
            else:
                # For arrays of equal start/end, compare the first differing point
                if np.array_repr(sl_rev) < np.array_repr(sl_norm):
                    target = sl_rev
                else:
                    target = sl_norm

            # Hash the byte representation
            sl_hash = target.tobytes()

            if sl_hash not in seen_hashes:
                seen_hashes.add(sl_hash)
                unique_streamlines.append(sl)

        return unique_streamlines

    def extract_mesh_from_msh(self, path_msh_file):

        gmsh.initialize()
        gmsh.open(path_msh_file)

        # 1. Extract All Nodes
        node_tags, coords, what = gmsh.model.mesh.getNodes()
        vertices = np.reshape(coords, (-1, 3))

        # Mapping GMSH-Node-Tag -> Index Vertices-Array
        max_node_tag = int(np.max(node_tags))
        node_tag_to_index = np.full(max_node_tag + 1, -1, dtype=np.int32)
        node_tag_to_index[node_tags.astype(int)] = np.arange(len(node_tags))

        # Get Nodes dimension (0,1,2 = Boundary, 3 = Interior)
        node_dim_array = np.full(max_node_tag + 1, -1, dtype=np.int32)
        for dim in range(4):
            for _, entity_tag in gmsh.model.getEntities(dim):
                tags, _, _ = gmsh.model.mesh.getNodes(dim, entity_tag)
                node_dim_array[tags.astype(int)] = dim

        # Geordnet nach Vertices-Array Index: shape (N_vertices,)
        vertex_dims = node_dim_array[node_tags.astype(int)]

        # 2. Get 2D-Surfaces
        surfaces = gmsh.model.getEntities(dim=2)
        surface_data = {}
        for _, s_tag in surfaces:
            quad_type = 3
            quad_tags, quad_node_tags = gmsh.model.mesh.getElementsByType(
                quad_type, tag=s_tag)

            if len(quad_tags) > 0:
                quad_indices = node_tag_to_index[quad_node_tags.astype(int)]
                quad_faces = quad_indices.reshape(-1, 4)
                quad_faces_tensor = quad_faces.T

                singular_nodes, singular_node_indices, singular_valences, boundary_nodes, unique_edges = self.get_singularities(
                    vertices, quad_faces_tensor)

                separatrices = self.trace_all_separatrices(
                    vertices, quad_faces_tensor, singular_node_indices, boundary_nodes, unique_edges
                )

                separatrices_ids = []

                for sepa in separatrices:
                    ids = np.array(sepa['path'])
                    separatrices_ids.append(ids)
                surface_data[s_tag] = {
                    "quad_faces": quad_faces_tensor,
                    "singular_nodes": singular_nodes,
                    "singular_node_indices": singular_node_indices,
                    "singular_valences": singular_valences,
                    "boundary_nodes": boundary_nodes,
                    "separatrices_data": separatrices,
                    "separatrices_ids": separatrices_ids,
                }

                streamlines = self.get_streamlines(surface_data, vertices)

                mesh_data = {
                    "vertices":   vertices,
                    "vertex_dims": vertex_dims,
                    "surfaces": surface_data,
                    "streamlines": streamlines,
                }
                print(f"Surface {s_tag}: {
                      quad_faces_tensor.shape[1]} Quads extrahiert.")

        return mesh_data

    def get_streamlines(self, surfaces_data, vertices):
        streamlines = {}

        for surface_tag in surfaces_data.keys():
            surface_data = surfaces_data[surface_tag]

            streamlines_surfaces_points = []
            streamlines_surfaces_splines = []

            for separatrix_id in surface_data['separatrices_ids']:

                points = vertices[separatrix_id]
                # tck         = points_to_spline(points)

                streamlines_surfaces_points.append(points)
                # streamlines_surfaces_splines.append(tck)

            streamlines_unique = self.filter_unique_streamlines(
                streamlines_surfaces_points)
            streamlines[surface_tag] = {
                'points': streamlines_unique,
                # 'splines' : streamlines_surfaces_splines,
            }

        return streamlines

    def export_streamlines_to_msh(self, streamlines, output_path="streamlines.msh"):
        """
        streamlines: Liste von numpy arrays, shape (N, 3)
        """
        gmsh.initialize()
        gmsh.model.add("streamlines")

        for line in streamlines:
            point_tags = []
            for pt in line:
                tag = gmsh.model.geo.addPoint(pt[0], pt[1], pt[2])
                point_tags.append(tag)

            # Spline durch alle Punkte der Stromlinie
            gmsh.model.geo.addSpline(point_tags)

        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(1)  # 1D mesh (nur Kurven)
        gmsh.write(output_path)
        gmsh.finalize()
