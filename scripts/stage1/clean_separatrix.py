#!/usr/bin/env python3
"""
A robust replacement for domain_partition's SeparatrixGenerator.

domain_partition's v1/v2 separatrix generators are buggy and were never
exercised on meshes with real singularities (the NACA data generation filters
to near-zero-singularity fields). We only need them to populate
``mesh.separatrices`` in the format that ``StreamlineGenerator_v2.get_streamlines``
consumes; that RK/Heun integrator itself is sound and is reused unchanged.

Emanation directions are found by walking a small circle around each singularity
and locating the angles where the (4-fold) cross field aligns radially - those
are exactly the separatrix departure directions. Boundary corner (dim-0) nodes
emit separatrices along the interior cross direction.

Interface (matches what StreamlineGenerator_v2 expects):
    mesh.singularities          (F,)   Poincare index per face   [reused logic]
    mesh.singularities_coords   {face: [x,y]}
    mesh.expected_separatrices  {face: count}
    mesh.separatrices           [ {coordinates, vector, singularity_coords} ]
"""

import numpy as np
import torch


def _poincare_indices(mesh):
    u = mesh.u
    ang = torch.atan2(u[:, 1], u[:, 0])
    f = mesh.faces
    a = ang[f[0]]
    b = ang[f[1]]
    c = ang[f[2]]
    def wrap(d):
        return torch.remainder(d + torch.pi, 2 * torch.pi) - torch.pi
    s = wrap(b - a) + wrap(c - b) + wrap(a - c)
    return torch.round(s / (2 * torch.pi))


def _singularity_coords(mesh, face_id):
    face = mesh.faces[:, face_id]
    vecs = mesh.u[face]
    nodes = mesh.x[face, 0:2]
    A = torch.stack([vecs[0] - vecs[2], vecs[1] - vecs[2]], dim=1)
    try:
        co = torch.linalg.solve(A, -vecs[2])
        if (co >= -0.2).all() and co.sum() <= 1.2:
            return nodes[2] + co[0] * (nodes[0] - nodes[2]) + co[1] * (nodes[1] - nodes[2])
    except RuntimeError:
        pass
    return nodes.mean(dim=0)


class CleanSeparatrixGenerator:
    def __init__(self, mesh, corner_merge_tol=0.06):
        self.mesh = mesh
        self.corner_merge_tol = corner_merge_tol
        self._build_locator()
        mesh.singularities = _poincare_indices(mesh)
        self._coords_and_expected()
        self._delete_corner_singularities()
        seps = self._emanate_from_singularities()
        seps += self._emanate_from_corners()
        mesh.separatrices = seps
        self.found_separatrices = seps

    def _delete_corner_singularities(self):
        """Drop singularities that sit right next to a c0 boundary corner.

        Slanted (non-90deg) domain corners get a bisector Dirichlet cross that
        conflicts with the wall-aligned representative of the adjacent boundary
        nodes, forcing a spurious singularity next to each corner. The corner
        itself is the real topological feature (it emits a separatrix), so we
        delete these artifacts (mirrors domain_partition's disabled
        delete_singularities)."""
        mesh = self.mesh
        corner_xy = self.nodes[mesh.x[:, 2].numpy() == 0]
        if len(corner_xy) == 0:
            return
        removed = 0
        for fid in list(mesh.singularities_coords.keys()):
            c = np.array(mesh.singularities_coords[fid])
            if np.min(np.linalg.norm(corner_xy - c, axis=1)) < self.corner_merge_tol:
                mesh.singularities[fid] = 0
                del mesh.singularities_coords[fid]
                mesh.expected_separatrices.pop(fid, None)
                removed += 1
        print(f"[clean_sep] deleted {removed} near-corner singularities, "
              f"{len(mesh.singularities_coords)} remain")

    # --- point location ---------------------------------------------------
    def _build_locator(self):
        self.nodes = self.mesh.x[:, 0:2].numpy().astype(np.float64)
        self.faces = self.mesh.faces.numpy().T  # (F,3)
        self.n2f = self.mesh.nodes_faces_ids
        self.u = self.mesh.u.numpy().astype(np.float64)

    def _find_face(self, p):
        d = np.sum((self.nodes - p) ** 2, axis=1)
        for fi in self.n2f.get(int(np.argmin(d)), []):
            tri = self.faces[fi]
            if self._bary_inside(p, self.nodes[tri]):
                return fi
        return None

    @staticmethod
    def _bary(p, tri):
        v0 = tri[1] - tri[0]
        v1 = tri[2] - tri[0]
        v2 = p - tri[0]
        den = v0[0] * v1[1] - v1[0] * v0[1]
        if abs(den) < 1e-18:
            return None
        a = (v2[0] * v1[1] - v1[0] * v2[1]) / den
        b = (v0[0] * v2[1] - v2[0] * v0[1]) / den
        return np.array([1 - a - b, a, b])

    def _bary_inside(self, p, tri):
        bc = self._bary(p, tri)
        return bc is not None and (bc >= -1e-9).all()

    def _cross_dirs(self, p):
        """Return the 4 unit cross directions of the field at point p (or None)."""
        fi = self._find_face(p)
        if fi is None:
            return None
        tri = self.faces[fi]
        bc = self._bary(p, self.nodes[tri])
        if bc is None:
            return None
        iv = bc @ self.u[tri]
        n = np.linalg.norm(iv)
        if n < 1e-9:
            return None
        base = np.arctan2(iv[1], iv[0]) / 4.0
        return base + np.arange(4) * (np.pi / 2)

    # --- singularity metadata --------------------------------------------
    def _coords_and_expected(self):
        mesh = self.mesh
        mesh.singularities_coords = {}
        mesh.expected_separatrices = {}
        face_ids = torch.where(mesh.singularities != 0)[0]
        for fid in face_ids.tolist():
            mesh.singularities_coords[fid] = _singularity_coords(mesh, fid).tolist()
            idx = int(mesh.singularities[fid].item())
            mesh.expected_separatrices[fid] = 3 if idx == 1 else (5 if idx == -1 else abs(4 * idx - 1))

    def _local_scale(self, fid):
        tri = self.nodes[self.faces[fid]]
        return float(np.mean([np.linalg.norm(tri[i] - tri[(i + 1) % 3]) for i in range(3)]))

    # --- emanation --------------------------------------------------------
    def _emanate_from_singularities(self, n_samples=720):
        seps = []
        mesh = self.mesh
        alphas = np.linspace(0, 2 * np.pi, n_samples, endpoint=False)
        for fid, coords in mesh.singularities_coords.items():
            c = np.array(coords)
            eps = 1.2 * self._local_scale(fid)
            g = np.full(n_samples, np.pi)
            for i, al in enumerate(alphas):
                p = c + eps * np.array([np.cos(al), np.sin(al)])
                cd = self._cross_dirs(p)
                if cd is None:
                    continue
                # angular distance from radial direction al to nearest cross dir
                diff = np.angle(np.exp(1j * (cd - al)))
                g[i] = np.min(np.abs(diff))
            # local minima of g below threshold = separatrix departure angles
            dirs = self._local_minima_angles(alphas, g, thresh=0.12)
            for al in dirs:
                d = np.array([np.cos(al), np.sin(al)])
                seps.append({
                    "coordinates": torch.tensor(c + eps * d, dtype=torch.float),
                    "vector": torch.tensor(d, dtype=torch.float),
                    "singularity_coords": torch.tensor(c, dtype=torch.float),
                    "face_id": fid,
                })
        return seps

    @staticmethod
    def _local_minima_angles(alphas, g, thresh):
        n = len(g)
        mins = []
        for i in range(n):
            if g[i] < thresh and g[i] <= g[(i - 1) % n] and g[i] < g[(i + 1) % n]:
                mins.append(alphas[i])
        # merge angles closer than 10 deg
        merged = []
        for a in mins:
            if all(abs(np.angle(np.exp(1j * (a - m)))) > np.radians(10) for m in merged):
                merged.append(a)
        return merged

    def _emanate_from_corners(self):
        seps = []
        mesh = self.mesh
        corner_ids = torch.where(mesh.x[:, 2] == 0)[0].tolist()
        for nid in corner_ids:
            c = self.nodes[nid]
            # inward direction = mean of edges to face-neighbours minus boundary
            faces = self.n2f.get(nid, [])
            if not faces:
                continue
            nbrs = set()
            for fi in faces:
                nbrs.update(self.faces[fi].tolist())
            nbrs.discard(nid)
            inward = np.mean([self.nodes[j] - c for j in nbrs], axis=0)
            ni = np.linalg.norm(inward)
            if ni < 1e-9:
                continue
            inward /= ni
            eps = 1.2 * np.mean([self._local_scale(fi) for fi in faces])
            cd = self._cross_dirs(c + eps * inward)
            if cd is None:
                continue
            # pick the cross direction best aligned with inward
            dots = np.cos(cd - np.arctan2(inward[1], inward[0]))
            al = cd[int(np.argmax(dots))]
            d = np.array([np.cos(al), np.sin(al)])
            if np.dot(d, inward) < 0:
                d = -d
            seps.append({
                "coordinates": torch.tensor(c + eps * d, dtype=torch.float),
                "vector": torch.tensor(d, dtype=torch.float),
                "singularity_coords": torch.tensor(c, dtype=torch.float),
                "face_id": -1 - nid,  # negative marker = boundary-corner origin
            })
        return seps
