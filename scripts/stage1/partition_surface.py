#!/usr/bin/env python3
"""
Stage 1 / Step C: Run domain_partition's 2D cross-field block partition on the
unwrapped surface mesh.

Pipeline (mirrors domain_partition/data_generator.py:get_mesh):
    FrameField -> detect_singularities -> StreamlineGenerator_v2
               -> StreamlinePostProcessor -> block_mesh

One deviation: StreamlineGenerator_v2 hardcodes the v1 ``SeparatrixGenerator``,
whose vector interpolation divides by a zero-norm field near singularities and
returns ``None`` (the NACA data-gen never hit it because it filters to
near-zero-singularity cases). We monkeypatch in ``SeparatrixGenerator_v2``,
which guards that case. No edits to the domain_partition repo.

Returns the 2D block mesh (block_mesh.x = (N,2) corners, block_mesh.faces =
(4,B) quad blocks) plus the affine (s,t)<-[0,1]^2 transform for back-mapping.
"""

import sys
from pathlib import Path

import numpy as np
import torch
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

sys.path.insert(0, "/root/repos/domain_partition")
sys.path.insert(0, str(Path(__file__).resolve().parent))

# --- swap in our robust separatrix emanation finder ---
# domain_partition's v1/v2 generators are buggy and untested on real
# singularities; only the RK/Heun integrator in StreamlineGenerator_v2 is sound.
import tools.streamline_generator_v2 as _sg2  # noqa: E402
from clean_separatrix import CleanSeparatrixGenerator  # noqa: E402
_sg2.SeparatrixGenerator = CleanSeparatrixGenerator

from tools import FrameField, StreamlineGenerator_v2, StreamlinePostProcessor  # noqa: E402
from tools.singularity_detector import detect_singularities  # noqa: E402
from tools.quad_partition_validator import QuadPartitionValidator  # noqa: E402
from tools.streamline_merging import StreamlineMerging  # noqa: E402
from dp_adapter import build_dp_data  # noqa: E402


# --- enable Xiao 2020 Algorithm 2 Case 2 in StreamlineMerging -----------------
# StreamlineMerging.find_missed_streamline_endpoints implements Case 2 (Fig 8b):
# a streamline that PASSES a singularity (angle ~180 deg to a reference out-going
# streamline, within 45 deg) but does not terminate there is CUT onto it. The
# upstream __init__ defines it but never calls it, so wrap separatrices that
# graze the opposite blade tip never connect and the O-grid ring never closes.
# We reinstate the call between prepare and merge (no edit to domain_partition).
def _streamline_merging_init(self, mesh, verbose=True):
    self.verbose = verbose
    self.Streamlines, self.Singularities = self.prepare_streamlines(mesh)
    self.Streamlines, self.Singularities = \
        self.find_missed_streamline_endpoints(self.Streamlines, self.Singularities)
    self.new_streamlines = self.merge_streamlines(self.Streamlines)


StreamlineMerging.__init__ = _streamline_merging_init


def _robust_find_containing_face(self, point, mesh):
    """Robust point location. The upstream version only tests faces incident to
    the single nearest node, so a streamline stepping into a non-incident face
    drops out of the mesh mid-domain and the separatrix dies early (-> dangling
    stubs that never reach their target). We test faces of the k nearest nodes,
    then fall back to a brute-force scan. Returns a face index or None."""
    p = point.detach().numpy() if hasattr(point, "detach") else np.asarray(point)
    nodes = mesh.x[:, 0:2].numpy()
    faces = mesh.faces.numpy()
    d = np.sum((nodes - p) ** 2, axis=1)
    near = np.argsort(d)[:6]
    cand = []
    for nid in near:
        cand.extend(mesh.nodes_faces_ids.get(int(nid), []))
    for fi in dict.fromkeys(cand):
        if _point_in_tri(p, nodes[faces[:, fi]]):
            return fi
    for fi in range(faces.shape[1]):  # brute fallback
        if _point_in_tri(p, nodes[faces[:, fi]]):
            return fi
    return None


def _point_in_tri(p, tri, eps=1e-9):
    v0 = tri[2] - tri[0]
    v1 = tri[1] - tri[0]
    v2 = p - tri[0]
    d00 = v0 @ v0; d01 = v0 @ v1; d02 = v0 @ v2
    d11 = v1 @ v1; d12 = v1 @ v2
    den = d00 * d11 - d01 * d01
    if abs(den) < 1e-18:
        return False
    u = (d11 * d02 - d01 * d12) / den
    v = (d00 * d12 - d01 * d02) / den
    return (u >= -eps) and (v >= -eps) and (u + v <= 1 + eps)


StreamlineGenerator_v2.find_containing_face = _robust_find_containing_face


def _add_cross_at_boundaries_fixed(self):
    """Corrected FrameField.add_cross_at_boundaries.

    The upstream version stores the boundary Dirichlet cross at array position
    ``i`` (the loop counter over boundary nodes), but the solver in
    ``compute_initial_frame_field`` reads ``frame_field_coords[node_id]``. These
    only agree when boundary node ids happen to be 0..K-1 (true for the NACA
    gmsh meshes, false for a welded STL with scattered ids) -> scrambled BC ->
    spurious singularities. We store at the node id instead. Geometry/angle
    logic is verbatim from the upstream method.
    """
    pi = torch.pi
    mesh = self.mesh
    num_nodes = mesh.x.size(0)
    mask_boundaryEdges = mesh.edge_attr == 1
    idx_boundaryNodes = torch.unique(mesh.edge_index[0, mask_boundaryEdges])
    boundary_edges = mesh.edge_index[:, mask_boundaryEdges]

    frame_field_angle = torch.zeros((num_nodes), dtype=torch.float)
    frame_field_coords = torch.zeros((num_nodes, 2), dtype=torch.float)

    for idx_current_node in idx_boundaryNodes:
        nid = int(idx_current_node)
        boundary_edges_of_node = torch.where(boundary_edges[0, :] == idx_current_node)[0]
        neighbours_idx = boundary_edges[1, boundary_edges_of_node]
        source_node = mesh.x[idx_current_node, 0:2]
        destination_node0 = mesh.x[neighbours_idx[0], 0:2]
        destination_node1 = mesh.x[neighbours_idx[1], 0:2]

        edge0 = destination_node0 - source_node
        edge1 = destination_node1 - source_node
        edge0_normalized = edge0 / torch.norm(edge0, p=2)
        edge1_normalized = edge1 / torch.norm(edge1, p=2)

        angle0 = torch.atan2(edge0_normalized[1], edge0_normalized[0]) % (2 * pi)
        angle1 = torch.atan2(edge1_normalized[1], edge1_normalized[0]) % (2 * pi)

        # Average in the REPRESENTATIVE (4-theta) space, not in tangent space.
        # The upstream code averages the two wall tangents (edge0+edge1) and only
        # THEN maps to the cross representative. At a slanted (non-90deg) corner
        # that tangent-space mean disagrees with both adjacent walls and forces a
        # spurious singularity right next to each c0 corner. The cross field lives
        # in representative space, so the mean must be taken there: map each wall
        # tangent to its representative, then circular-mean them. For smooth
        # boundary nodes the two edges are ~collinear and both methods coincide.
        r0 = self.map_cross_vectors_to_reference_vector(angle0)
        r1 = self.map_cross_vectors_to_reference_vector(angle1)
        ref_angle = torch.atan2(torch.sin(r0) + torch.sin(r1),
                                torch.cos(r0) + torch.cos(r1))

        frame_field_angle[nid] = ref_angle
        frame_field_coords[nid, 0] = torch.cos(ref_angle)
        frame_field_coords[nid, 1] = torch.sin(ref_angle)

    self.mesh.frame_field_angle = frame_field_angle
    self.mesh.frame_field_coords = frame_field_coords


FrameField.add_cross_at_boundaries = _add_cross_at_boundaries_fixed


# --- optional soft boundary conditions (the "energy" knob) -----------------
# BC_WEIGHT = None reproduces the upstream HARD Dirichlet BC (row-replacement).
# A finite w adds a penalty term  w*|u - u_wall|^2  to the harmonic energy
# instead of hard-fixing the boundary DOFs: the harmonic smoothing then governs
# the interior more strongly and tends to merge/annihilate nearby opposite cross
# singularities. Lower w = softer boundary = smoother interior (but weaker wall
# alignment). Note generate_cross_field re-imposes the exact BC into
# mesh.frame_field; mesh.u (used downstream for singularities/separatrices) keeps
# the soft solve, which is exactly the relaxed interior we want.
BC_WEIGHT = None

# --- pitchwise periodicity (Phase 2) ---------------------------------------
# When True and mesh.periodic_pairs is non-empty, the two theta walls are NOT
# wall-aligned Dirichlet boundaries; instead the seam is glued: each slave node
# DOF is identified with its master partner (u_slave = u_master) and the slave's
# triangle stiffness is welded onto the master, giving a true periodic Laplacian.
# This removes the artificial wall holonomy on the pitch boundaries (Kowalski Eq
# 14 budget) that forces extra interior singularities. Kept behind a flag so the
# wall-BC result stays reproducible for comparison.
PERIODIC = True


def set_bc_weight(w):
    global BC_WEIGHT
    BC_WEIGHT = w


def set_periodic(flag):
    global PERIODIC
    PERIODIC = bool(flag)


def _seam_weld_map(mesh, num_nodes):
    """Return (weld, seam, slave_of) for the periodic seam.

    weld: (num_nodes,) int array remapping each slave node id to its master.
    seam: set of all node ids on the seam (masters + slaves).
    slave_of: {slave_id: master_id}."""
    weld = np.arange(num_nodes)
    seam, slave_of = set(), {}
    pairs = getattr(mesh, "periodic_pairs", None)
    if not PERIODIC or pairs is None or len(pairs) == 0:
        return weld, seam, slave_of
    pp = pairs.numpy() if hasattr(pairs, "numpy") else np.asarray(pairs)
    for m, s in pp:
        m, s = int(m), int(s)
        weld[s] = m
        slave_of[s] = m
        seam.add(m); seam.add(s)
    return weld, seam, slave_of


def _compute_initial_frame_field_soft(self):
    """FrameField.compute_initial_frame_field with optional soft (penalty) BC and
    optional pitchwise periodicity.

    BC_WEIGHT: None -> hard Dirichlet wall BC, finite -> soft penalty BC.
    PERIODIC: glue the theta-seam (weld slave stiffness onto master + constrain
    u_slave = u_master), and skip wall-Dirichlet on seam nodes."""
    num_nodes = self.mesh.x.shape[0]
    num_elements = self.mesh.faces.shape[1]
    nodes = self.mesh.x[:, 0:2]
    elements = self.mesh.faces.T
    mask_boundaryEdges = self.mesh.edge_attr == 1
    boundary_nodes_indices = torch.unique(self.mesh.edge_index[0, mask_boundaryEdges])
    num_dofs = num_nodes * 2
    b = np.zeros(num_dofs)

    weld, seam, slave_of = _seam_weld_map(self.mesh, num_nodes)
    self._seam_nodes = seam

    rows, cols, vals = [], [], []
    for e in range(num_elements):
        nodes_indices = elements[e]
        coords = nodes[nodes_indices]
        A_e = self.compute_local_stiffness_matrix(coords)
        dof_indices = np.empty(6, dtype=int)
        for i in range(3):
            wn = int(weld[int(nodes_indices[i])])  # periodic remap (identity if none)
            dof_indices[2 * i] = 2 * wn
            dof_indices[2 * i + 1] = 2 * wn + 1
        for i_local in range(6):
            for j_local in range(6):
                rows.append(dof_indices[i_local])
                cols.append(dof_indices[j_local])
                vals.append(A_e[i_local, j_local])

    A = coo_matrix((vals, (rows, cols)), shape=(num_dofs, num_dofs)).tocsr().tolil()
    for idx in boundary_nodes_indices:
        idx = int(idx)
        if idx in seam:           # periodic seam: no wall Dirichlet
            continue
        dof_x, dof_y = 2 * idx, 2 * idx + 1
        bx = float(self.mesh.frame_field_coords[idx, 0])
        by = float(self.mesh.frame_field_coords[idx, 1])
        if BC_WEIGHT is None:
            A.rows[dof_x] = [dof_x]; A.data[dof_x] = [1.0]
            A.rows[dof_y] = [dof_y]; A.data[dof_y] = [1.0]
            b[dof_x] = bx; b[dof_y] = by
        else:
            w = float(BC_WEIGHT)
            A[dof_x, dof_x] += w
            A[dof_y, dof_y] += w
            b[dof_x] += w * bx
            b[dof_y] += w * by

    # periodic constraint rows: u_slave - u_master = 0 (slave stiffness already
    # welded onto master above, so the master row carries the periodic stencil).
    for s, m in slave_of.items():
        for comp in (0, 1):
            ds, dm = 2 * s + comp, 2 * m + comp
            A.rows[ds] = [ds, dm]
            A.data[ds] = [1.0, -1.0]
            b[ds] = 0.0

    A_sparse = A.tocsr()
    u = spsolve(A_sparse, b)
    return A_sparse, b, u


FrameField.compute_initial_frame_field = _compute_initial_frame_field_soft


def _generate_cross_field_seam_aware(self):
    """FrameField.generate_cross_field, but the boundary re-imposition skips the
    periodic seam nodes. The upstream method overwrites every boundary node in
    mesh.frame_field with its wall-aligned Dirichlet cross; on the periodic seam
    that would reintroduce the wall holonomy into singularity detection. We keep
    the solved periodic value on the seam instead."""
    A, b, u = self.compute_initial_frame_field()
    u_new = self.Linearization_Norm_Constraint(A, b, u)

    u_init = torch.cat((torch.from_numpy(u[::2]).unsqueeze(1),
                        torch.from_numpy(u[1::2]).unsqueeze(1)), dim=1)
    u_final = torch.cat((torch.from_numpy(u_new[::2]).unsqueeze(1),
                         torch.from_numpy(u_new[1::2]).unsqueeze(1)), dim=1)

    seam = getattr(self, "_seam_nodes", set())
    vector_field = u_final.detach().clone().to(torch.float)
    bidx = torch.where(self.mesh.x[:, 2] != 2)[0].tolist()
    for i in bidx:
        if i in seam:             # periodic seam keeps the solved value
            continue
        vector_field[i, :] = self.mesh.frame_field_coords[i, :]
    self.mesh.u_init = u_init
    self.mesh.u = u_final
    self.mesh.frame_field = vector_field


FrameField.generate_cross_field = _generate_cross_field_seam_aware


def partition(stl_path, verbose=True, bc_weight=None):
    """Run cross-field block partition on the unwrapped surface.

    bc_weight: None = hard Dirichlet boundary (default); a finite float = soft
    penalty BC (the energy knob, see set_bc_weight) that relaxes the interior to
    merge spurious cross singularities.

    Returns (block_mesh, mesh, transform):
      block_mesh.x      (Nc,2) block corner nodes in normalized [0,1]^2
      block_mesh.faces  (4,B)  quad block connectivity
      mesh              the triangulated mesh w/ field, separatrices, streamlines
      transform         affine map normalized [0,1]^2 -> (s,t)
    """
    set_bc_weight(bc_weight)
    mesh, transform = build_dp_data(stl_path)

    ff = FrameField(mesh)
    m = detect_singularities(ff.mesh)
    n_sing = int((m.singularities != 0).sum())

    sl = StreamlineGenerator_v2(ff.mesh)
    _snap_separatrix_endpoints(sl.mesh, radius=0.045)
    pp = StreamlinePostProcessor(sl.mesh, verbose=False)
    block_mesh = pp.block_mesh
    n_blocks = block_mesh.faces.shape[1] if block_mesh.faces is not None else 0
    if verbose:
        print(f"[partition] singularities={n_sing}  "
              f"separatrices={len(sl.mesh.separatrices)}  "
              f"-> {n_blocks} quad blocks ({block_mesh.x.shape[0]} corners)")
    return block_mesh, sl.mesh, transform


def _validate_separatrix_counts(mesh):
    """Check the per-singularity separatrix-count invariant (Kowalski 2015): a
    field singularity of representation index +1 must keep EXACTLY 3 incident
    separatrices, index -1 EXACTLY 5. A violation after postprocessing is the
    direct cause of a non-four-sided region (irregular interior block node /
    Euler defect). Returns (all_ok, report_lines)."""
    sing_fids = list(getattr(mesh, "singularities_coords", {}).keys())
    if not sing_fids:
        return True, ["separatrix-count: no field singularities"]
    n_sing = len(sing_fids)
    got = {k: 0 for k in range(n_sing)}
    for sd in mesh.separatrices:
        for key in ("origin_sing", "target_sing"):
            kk = sd.get(key) if isinstance(sd, dict) else None
            if kk is not None and 0 <= kk < n_sing:
                got[kk] += 1
    all_ok = True
    lines = ["separatrix-count per singularity (Kowalski 3/5 invariant):"]
    for k, fid in enumerate(sing_fids):
        idx = int(mesh.singularities[fid].item())
        exp = int(mesh.expected_separatrices.get(fid, 4))
        g = got[k]
        ok = (g == exp)
        all_ok = all_ok and ok
        lines.append(f"  fid={fid} idx={idx:+d} expected={exp} got={g} "
                     f"{'OK' if ok else 'MISMATCH'}")
    lines.append(f"separatrix-count invariant: {'PASS' if all_ok else 'FAIL'}")
    return all_ok, lines


def validate(block_mesh, mesh, out_path=None):
    """Run domain_partition's QuadPartitionValidator on the finished block_mesh.

    Validation is only defined AFTER postprocessing has produced quad faces; if
    postproc yielded nothing, that is itself a validation failure (reported, not
    crashed). frame_field is the per-node cross representation mesh.u (cos4t/sin4t),
    needed for the Kowalski singularity-efficiency metric."""
    lines = []
    if (block_mesh is None or block_mesh.faces is None
            or block_mesh.faces.numel() == 0):
        lines.append("VALIDATION FAIL: postprocessing produced no quad blocks.")
        report = "\n".join(lines)
        if out_path:
            Path(out_path).write_text(report)
        print(report)
        return None

    v = QuadPartitionValidator(block_mesh, mesh, frame_field=mesh.u, strict=True)
    is_valid = v.is_valid()
    qs = v.quality_score()
    soft = v.passes_soft_thresholds()
    diag = v.diagnostics()

    lines.append(f"is_valid (hard checks): {is_valid}")
    lines.append(f"passes_soft_thresholds: {soft}")
    lines.append(f"blocks: {block_mesh.faces.shape[1]}  "
                 f"corners: {block_mesh.x.shape[0]}")
    exp = getattr(v, "_expected_singularities", None)
    act = getattr(v, "_actual_singularities", None)
    if exp is not None:
        eff = qs.get("singularity_efficiency")
        lines.append(f"singularities expected/actual: {exp}/{act}  "
                     f"efficiency: {eff}")
    lines.append("quality_score:")
    for k, val in qs.items():
        lines.append(f"  {k}: {val}")
    lines.append("diagnostics:")
    lines.extend(f"  - {d}" for d in (diag or ["(none)"]))

    sep_ok, sep_lines = _validate_separatrix_counts(mesh)
    lines.extend(sep_lines)

    report = "\n".join(lines)
    if out_path:
        Path(out_path).write_text(report)
        print(f"wrote {out_path}")
    print(report)
    return v


def _termination_nodes(mesh):
    """Targets a separatrix may end on: singularities + c0 corners."""
    pts = [np.array(c) for c in mesh.singularities_coords.values()]
    corners = mesh.x[mesh.x[:, 2] == 0, 0:2].numpy()
    pts.extend(list(corners))
    return np.array(pts) if pts else np.zeros((0, 2))


def _project_to_polyline(p, poly):
    """Closest point on a polyline. Returns (seg_index, t, dist, proj_point)."""
    best = (None, 0.0, np.inf, None)
    for i in range(len(poly) - 1):
        a, b = poly[i], poly[i + 1]
        ab = b - a
        L2 = ab @ ab
        t = 0.0 if L2 < 1e-18 else np.clip((p - a) @ ab / L2, 0.0, 1.0)
        proj = a + t * ab
        d = np.linalg.norm(p - proj)
        if d < best[2]:
            best = (i, t, d, proj)
    return best


def _min_boundary_dist(p, boundary):
    """Distance from point p to the nearest boundary polyline (incl. blade loop)."""
    p = np.asarray(p, float)
    best = np.inf
    for poly in boundary:
        d = _project_to_polyline(p, poly)[2]
        if d < best:
            best = d
    return best


def _snap_separatrix_endpoints(mesh, radius=0.045, bnd_radius=0.05, end_frac=0.3):
    """Kowalski post-processing: every separatrix must end on a singularity, a
    c0 corner, or the boundary dOmega (outer *and* the inner blade loop). The RK
    integrator drifts past targets and the block graph only treats streamline
    *endpoints* as nodes, so we:
      (1) truncate+snap a separatrix to the first point-target (singularity/
          corner) it enters, else
      (2) snap its dangling end onto the nearest boundary polyline and *split*
          that boundary streamline at the snap point, so the new T-junction is a
          shared graph node and blocks can close against it (notably the blade).
    """
    nodes = _termination_nodes(mesh)
    n_sing = len(mesh.singularities_coords)  # nodes[:n_sing] are singularities
    n_boundary = len(mesh.streamlines) - len(mesh.separatrices)
    boundary = [np.asarray(s, float) for s in mesh.streamlines[:n_boundary]]
    seps = [np.asarray(s, float) for s in mesh.streamlines[n_boundary:]]
    sep_dicts = list(mesh.separatrices)

    def _origin_sing(si):
        """node index of a separatrix's origin if it is a field singularity."""
        if si >= len(sep_dicts) or int(sep_dicts[si].get("face_id", -1)) < 0:
            return None  # blade-tip / corner origin, not a field singularity
        oc = np.asarray(sep_dicts[si]["singularity_coords"], float)
        if n_sing == 0:
            return None
        k = int(np.argmin(np.linalg.norm(nodes[:n_sing] - oc, axis=1)))
        return k

    # --- pass 1: snap each separatrix, record metadata (no boundary split yet) -
    recs = []  # per-sep: dict(poly, origin_sing, target_sing, bnd=(bi,segpos,proj))
    n_pt, n_bnd = 0, 0
    for si, s in enumerate(seps):
        rec = {"poly": s, "origin_sing": _origin_sing(si),
               "target_sing": None, "bnd": None}
        if s.ndim != 2 or len(s) < 2:
            recs.append(rec)
            continue
        origin = s[0]
        origin_node = int(np.argmin(np.linalg.norm(nodes - origin, axis=1))) \
            if len(nodes) else -1
        # (1) point-target snapping along the path.
        # Blade-vs-neighbour guard: the two tip singularities of a cluster sit
        # ~0.03 apart (below `radius`). A short separatrix emanating from one tip
        # toward the blade ends ON the blade loop right between the two tips, but
        # the neighbouring tip singularity is within `radius` of that end, so the
        # naive nearest-node snap truncates it at the neighbour instead of letting
        # it dock on the blade. Skip such a singularity snap when the approach
        # point is closer to a boundary than to that singularity, lies within
        # bnd_radius of the boundary, AND is near the END of the separatrix
        # (little arclength remaining). The last clause is essential: a runaway
        # multi-wrap spiral around the blade also hugs the blade, but it grazes
        # the neighbour at its START (lots of path remaining) -- there we DO keep
        # the sing snap so the spiral is discarded, not kept.
        seg_all = np.linalg.norm(np.diff(s, axis=0), axis=1)
        total_len = float(seg_all.sum())
        cut = False
        for j in range(1, len(s)):
            if len(nodes) == 0:
                break
            d = np.linalg.norm(nodes - s[j], axis=1)
            k = int(np.argmin(d))
            if d[k] < radius and k != origin_node:
                if k < n_sing and total_len > 1e-9:
                    dbl = _min_boundary_dist(s[j], boundary)
                    rem = float(seg_all[j:].sum()) / total_len
                    if dbl < d[k] and dbl < bnd_radius and rem < end_frac:
                        continue               # terminates on the blade, not S_k
                s = np.vstack([s[:j], nodes[k]])
                cut = True
                n_pt += 1
                if k < n_sing:                 # snapped onto a field singularity
                    rec["target_sing"] = k
                break
        # (2) boundary snapping of the (still dangling) end
        if not cut:
            end = s[-1]
            best = (None, None, np.inf, None)
            for bi, poly in enumerate(boundary):
                seg, t, dist, proj = _project_to_polyline(end, poly)
                if dist < best[2]:
                    best = (bi, (seg, t), dist, proj)
            if best[0] is not None and best[2] < bnd_radius:
                bi, (seg, t), _, proj = best
                s = np.vstack([s[:-1], proj])
                rec["bnd"] = (bi, seg + t, proj)
                n_bnd += 1
        rec["poly"] = s
        recs.append(rec)

    # --- blend separatrices joining two singularities (Kowalski 2015, sec 4.3) -
    # A separatrix S0->S1 between two field singularities is integrated TWICE:
    # once from S0 (rec i: origin S0, target S1) and once from S1 (rec j: origin
    # S1, target S0). Round-off makes the two paths differ. Kowalski Eq 28 does
    # NOT delete one of them; it linearly BLENDS them into a single, more accurate
    # curve, trusting each path near the singularity it spawned from:
    #     gamma_b(s) = (1 - s) * gamma1(s) + s * gamma2(s),   s in [0, 1]
    # with gamma1: S0->S1 and gamma2: S1->S0 (read reversed, i.e. S0->S1). We keep
    # the blended curve as rec i and drop rec j. This yields exactly one clean
    # separatrix per singularity-singularity connection (correct 3/5 valence) with
    # no duplicate and no residual integration error.
    def _resample(poly, n):
        poly = np.asarray(poly, float)
        if len(poly) < 2:
            return None
        seg = np.linalg.norm(np.diff(poly, axis=0), axis=1)
        d = np.concatenate([[0.0], np.cumsum(seg)])
        if d[-1] < 1e-12:
            return None
        ss = np.linspace(0.0, 1.0, n)
        return np.column_stack([np.interp(ss, d / d[-1], poly[:, 0]),
                                np.interp(ss, d / d[-1], poly[:, 1])])

    N_BLEND = 60
    drop = set()
    for i, ri in enumerate(recs):
        if i in drop:
            continue
        A, B = ri["origin_sing"], ri["target_sing"]
        if A is None or B is None or A == B:
            continue
        # matching reverse separatrix: origin B, target A (Case 1, Fig 8a)
        j = next((k for k, rk in enumerate(recs)
                  if k not in drop and k != i
                  and rk["origin_sing"] == B and rk["target_sing"] == A), None)
        if j is None:
            continue
        g1 = _resample(ri["poly"], N_BLEND)          # S0 -> S1
        g2 = _resample(recs[j]["poly"], N_BLEND)     # S1 -> S0
        if g1 is None or g2 is None:
            continue
        g2 = g2[::-1]                                # read as S0 -> S1
        s = np.linspace(0.0, 1.0, N_BLEND)[:, None]
        ri["poly"] = (1.0 - s) * g1 + s * g2         # Kowalski Eq 28 / Xiao Eq 10
        drop.add(j)

    # --- Xiao 2020 Case 2 (Fig 8b): geometric reverse-duplicate wraps -----------
    # A wrap along one side of the blade connects a TE tip to an LE tip. Because
    # each tip is a *cluster* of two singularities ~0.03 apart, the same physical
    # wrap gets emitted from BOTH ends but the two integrations snap to DIFFERENT
    # cluster-neighbours (e.g. S0->S2 one way, S3->S0 the other). The Case-1 blend
    # above only matches the mutual A->B / B->A pair, so it misses this: the two
    # copies survive and over-count their shared tip (Kowalski 3/5 invariant then
    # fails, e.g. a -1 tip gets 6 separatrices instead of 5). Detect the duplicate
    # geometrically (near-identical curve, orientation-independent Hausdorff) and
    # drop ONE copy -- specifically the copy whose removal keeps BOTH its endpoint
    # singularities at or above their required count (Xiao/Kowalski count-safe
    # dedup). Blend the survivor with the reversed dropped copy for accuracy.
    sfids = list(mesh.singularities_coords.keys())

    def _expected(k):
        if k is None or k >= len(sfids):
            return 4
        fid = sfids[k]
        idx = int(mesh.singularities[fid].item())
        return int(getattr(mesh, "expected_separatrices", {}).get(
            fid, 3 if idx == 1 else 5 if idx == -1 else 4))

    def _got(exclude=()):
        g = {}
        for ii, rr in enumerate(recs):
            if ii in drop or ii in exclude:
                continue
            for key in ("origin_sing", "target_sing"):
                kk = rr[key]
                if kk is not None:
                    g[kk] = g.get(kk, 0) + 1
        return g

    def _set_hausdorff(a, b):
        # orientation-independent max of nearest-neighbour distances
        d = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
        return max(d.min(axis=1).max(), d.min(axis=0).max())

    DUP_TOL = 0.08
    sing_recs = [i for i, r in enumerate(recs)
                 if i not in drop and r["origin_sing"] is not None
                 and r["target_sing"] is not None
                 and r["origin_sing"] != r["target_sing"]]
    for a_i in range(len(sing_recs)):
        i = sing_recs[a_i]
        if i in drop:
            continue
        gi = _resample(recs[i]["poly"], N_BLEND)
        if gi is None:
            continue
        for b_i in range(a_i + 1, len(sing_recs)):
            j = sing_recs[b_i]
            if j in drop:
                continue
            gj = _resample(recs[j]["poly"], N_BLEND)
            if gj is None or _set_hausdorff(gi, gj) > DUP_TOL:
                continue
            got = _got()
            exp = {k: _expected(k) for k in got}

            def _safe_to_drop(m):
                for key in ("origin_sing", "target_sing"):
                    k = recs[m][key]
                    if k is not None and got.get(k, 0) - 1 < exp.get(k, 4):
                        return False
                return True

            victim = i if _safe_to_drop(i) else (j if _safe_to_drop(j) else None)
            if victim is None:
                continue
            keep = j if victim == i else i
            gk = _resample(recs[keep]["poly"], N_BLEND)
            gv = _resample(recs[victim]["poly"], N_BLEND)
            if gk is not None and gv is not None:
                if _set_hausdorff(gk[:1], gv[:1]) > _set_hausdorff(gk[:1], gv[-1:]):
                    gv = gv[::-1]              # align victim to survivor direction
                ss = np.linspace(0.0, 1.0, N_BLEND)[:, None]
                recs[keep]["poly"] = (1.0 - ss) * gk + ss * gv
            drop.add(victim)
            if victim == i:
                break                          # i gone, move to next survivor

    kept = [i for i in range(len(recs)) if i not in drop]

    # stamp origin/target singularity back onto the surviving separatrix dicts
    # so the count validator (and any downstream graph) can read the connectivity
    for i in kept:
        sep_dicts[i]["origin_sing"] = recs[i]["origin_sing"]
        sep_dicts[i]["target_sing"] = recs[i]["target_sing"]

    # --- pass 2: boundary T-junction splits, only from KEPT separatrices -------
    splits = {i: [] for i in range(n_boundary)}
    for i in kept:
        if recs[i]["bnd"] is not None:
            bi, segpos, proj = recs[i]["bnd"]
            splits[bi].append((segpos, proj))

    new_boundary = []
    for bi, poly in enumerate(boundary):
        pts = sorted(splits[bi])
        if not pts:
            new_boundary.append(poly)
            continue
        cuts = []
        cur = [poly[0]]
        ptr = 0
        for i in range(len(poly) - 1):
            cur.append(poly[i + 1])
            while ptr < len(pts) and seg_floor(pts[ptr][0]) == i:
                proj = pts[ptr][1]
                cur[-1] = proj  # end this sub-streamline at the T-junction
                cuts.append(np.array(cur))
                cur = [proj, poly[i + 1]]
                ptr += 1
        cuts.append(np.array(cur))
        new_boundary.extend(c for c in cuts if len(c) >= 2)

    mesh.streamlines = new_boundary + [recs[i]["poly"] for i in kept]
    mesh.separatrices = [sep_dicts[i] for i in kept]
    print(f"[snap] point-snapped {n_pt}, boundary-snapped {n_bnd}; "
          f"deduped {len(drop)} matching prongs; "
          f"boundary {n_boundary}->{len(new_boundary)} segments, "
          f"separatrices {len(seps)}->{len(kept)}")


def seg_floor(x):
    return int(np.floor(x))


def _plot_separatrices(mesh, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    xy = mesh.x[:, 0:2].numpy()
    tris = mesh.faces.T.numpy()
    n_b = len(mesh.streamlines) - len(mesh.separatrices)
    fig, ax = plt.subplots(figsize=(9, 8))
    ax.triplot(xy[:, 0], xy[:, 1], tris, lw=0.1, color="0.92")
    for i, s in enumerate(mesh.streamlines):
        s = np.asarray(s, dtype=float)
        if s.ndim != 2 or len(s) < 2:
            continue
        ax.plot(s[:, 0], s[:, 1], "0.4" if i < n_b else "C3",
                lw=1.0 if i < n_b else 1.4)
    sc = xy[tris[mesh.singularities.numpy() != 0]].mean(axis=1)
    ax.scatter(sc[:, 0], sc[:, 1], c="blue", s=70, zorder=6, label="singularity")
    cm = mesh.x[:, 2].numpy() == 0
    ax.scatter(xy[cm, 0], xy[cm, 1], c="green", s=60, marker="^", zorder=6,
               label="c0 corner")
    ax.set_aspect("equal")
    ax.legend()
    ax.set_title(f"{len(mesh.separatrices)} separatrices, "
                 f"{int((mesh.singularities.numpy()!=0).sum())} singularities")
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    print(f"wrote {out_png}")


def _block_annotations(bx, bf, ratio_max=10.0):
    """Per-block-mesh validation annotations for plotting:
    inverted quads, interior nodes with valence != 4, high-aspect quads."""
    # edge topology
    from collections import defaultdict
    cnt = defaultdict(int)
    deg = np.zeros(len(bx), dtype=int)
    for quad in bf:
        for k in range(4):
            a, b = quad[k], quad[(k + 1) % 4]
            cnt[(min(a, b), max(a, b))] += 1
    for (a, b), c in cnt.items():
        deg[a] += 1
        deg[b] += 1
    boundary_nodes = set()
    for (a, b), c in cnt.items():
        if c == 1:
            boundary_nodes.add(a)
            boundary_nodes.add(b)
    irregular = [i for i in range(len(bx))
                 if i not in boundary_nodes and deg[i] != 4]
    inverted, high_aspect = [], []
    for qi, quad in enumerate(bf):
        p = bx[quad]
        # signed area (shoelace); negative => inverted ordering
        area = 0.5 * np.sum(p[:, 0] * np.roll(p[:, 1], -1)
                            - np.roll(p[:, 0], -1) * p[:, 1])
        if area <= 0:
            inverted.append(qi)
        elens = [np.linalg.norm(p[(k + 1) % 4] - p[k]) for k in range(4)]
        if min(elens) > 1e-12 and max(elens) / min(elens) > ratio_max:
            high_aspect.append(qi)
    return irregular, inverted, high_aspect


def _plot_blocks(block_mesh, mesh, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    xy = mesh.x[:, 0:2].numpy()
    tris = mesh.faces.T.numpy()
    bx = block_mesh.x.numpy()
    bf = block_mesh.faces.numpy().T  # (B,4)
    irregular, inverted, high_aspect = _block_annotations(bx, bf)
    fig, ax = plt.subplots(figsize=(9, 7))
    ax.triplot(xy[:, 0], xy[:, 1], tris, lw=0.1, color="0.92")
    for qi, quad in enumerate(bf):
        ring = bx[list(quad) + [quad[0]]]
        if qi in inverted:
            ax.fill(ring[:, 0], ring[:, 1], color="red", alpha=0.4)
            ax.plot(ring[:, 0], ring[:, 1], "red", lw=2.0)
        elif qi in high_aspect:
            ax.fill(ring[:, 0], ring[:, 1], color="orange", alpha=0.3)
            ax.plot(ring[:, 0], ring[:, 1], "darkorange", lw=1.8)
        else:
            ax.fill(ring[:, 0], ring[:, 1], alpha=0.25)
            ax.plot(ring[:, 0], ring[:, 1], "C0", lw=1.3)
    ax.scatter(bx[:, 0], bx[:, 1], c="k", s=8, zorder=5)
    if irregular:
        ax.scatter(bx[irregular, 0], bx[irregular, 1], facecolors="none",
                   edgecolors="red", s=160, linewidths=2.0, zorder=6,
                   label=f"irregular interior node (valence!=4): {len(irregular)}")
    if high_aspect:
        ax.plot([], [], "darkorange", lw=1.8,
                label=f"aspect>10: {len(high_aspect)}")
    if inverted:
        ax.plot([], [], "red", lw=2.0, label=f"inverted: {len(inverted)}")
    if irregular or high_aspect or inverted:
        ax.legend(loc="upper left", fontsize=8)
    ax.set_aspect("equal")
    ax.set_title(f"Hub quad block partition: {bf.shape[0]} blocks  "
                 f"(irregular={len(irregular)}, inverted={len(inverted)}, "
                 f"aspect>10={len(high_aspect)})")
    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    print(f"wrote {out_png}")


if __name__ == "__main__":
    stl = sys.argv[1] if len(sys.argv) > 1 else \
        "/root/repos/block_structured_meshing/T1_9_hub_raw.stl"
    out = Path("/root/repos/block_structured_meshing/output/T1_9/hub_stage1")
    out.mkdir(parents=True, exist_ok=True)
    block_mesh, mesh, tf = partition(stl)
    _plot_blocks(block_mesh, mesh, out / "blocks_2d.png")
    _plot_separatrices(mesh, out / "separatrices.png")
    validate(block_mesh, mesh, out / "validation_report.txt")
