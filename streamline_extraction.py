import gmsh
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.colors import to_rgba


class StrealimeExtractor:

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

                mesh_data = {
                    "vertices":   vertices,
                    "vertex_dims": vertex_dims,
                    "surfaces": surface_data,
                }
                print(f"Surface {s_tag}: {
                      quad_faces_tensor.shape[1]} Quads extrahiert.")

        return mesh_data

    def export_separatrices_to_msh(self, data, name="separatrices"):
        """
        Exportiert Separatrizen als Gmsh .msh Datei (Format 2.2).
        Jede Separatrix bekommt einen eigenen Physical Tag.
        """
        # Alle Punkte und Linien sammeln
        all_points = []       # Liste von (x, y, z)
        all_elements = []     # Liste von (node1_idx, node2_idx)
        point_map = {}        # node-index -> msh-node-index (1-basiert)

        msh_node_id = 1
        for key in data:
            filename = f"{name}_{key}.msh"
            nodes = data[key]["vertices"]
            separatrices = data[key]["separatrices"]
            for sep in separatrices:
                path = sep['path']
                for ni in path:
                    if ni not in point_map:
                        point_map[ni] = msh_node_id
                        all_points.append(nodes[ni, :3])
                        msh_node_id += 1
                # Liniensegmente
                for k in range(len(path) - 1):
                    all_elements.append(
                        (point_map[path[k]], point_map[path[k+1]]))

            with open(filename, 'w') as f:
                # Header
                f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")

                # Nodes
                f.write("$Nodes\n")
                f.write(f"{len(all_points)}\n")
                for i, (x, y, z) in enumerate(all_points, start=1):
                    f.write(f"{i} {x:.10e} {y:.10e} {z:.10e}\n")
                f.write("$EndNodes\n")

                # Elements (Typ 1 = 2-Knoten-Linie)
                f.write("$Elements\n")
                f.write(f"{len(all_elements)}\n")
                for i, (n1, n2) in enumerate(all_elements, start=1):
                    # Format: elem_id  elem_type  n_tags  tag1  tag2  node1  node2
                    f.write(f"{i} 1 2 1 0 {n1} {n2}\n")
                f.write("$EndElements\n")

            print(f"Gespeichert: {filename}  ({len(all_points)} Knoten, {
                  len(all_elements)} Segmente)")
