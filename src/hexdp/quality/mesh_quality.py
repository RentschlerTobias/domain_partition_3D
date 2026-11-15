"""CFD mesh quality metrics for hexahedral meshes.

The metrics a CFD solver actually cares about, computed on our own meshes
rather than quoted from a table. Definitions and default limits follow
OpenFOAM's `checkMesh` / `meshQualityDict`, because that is the solver this
project feeds (`etc/caseDicts/meshQualityDict`: `maxNonOrtho 65`,
`maxBoundarySkewness 20`, `maxInternalSkewness 4`, `minVolRatio 0.01`,
`minDeterminant 0.001`, `minFaceWeight 0.05`, `minTwist 0.02`).

Why this exists separately from `HexBlockValidator`: that one checks the
*block structure* (is every block a topological cuboid). This one checks the
*mesh* a solver will integrate over. A perfect block structure can still
carry cells no solver will tolerate, and the TFI stage is where that is
decided — so the gate has to exist before TFI, not after.

Every metric is written as a cell or face field into a VTK so it can be
thresholded in ParaView, exactly like the showcase steps.
"""

import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "experimentell" / "hex3d_algohex"))

import base_complex as bc                                             # noqa: E402
import clean_blocks as cb                                             # noqa: E402

HF = bc.HEX_FACES
HE = bc.HEX_EDGES

# OpenFOAM meshQualityDict defaults, for reporting against
LIMITS = {
    "non_orthogonality_deg": 65.0,      # maxNonOrtho
    "skewness_internal": 4.0,           # maxInternalSkewness
    "skewness_boundary": 20.0,          # maxBoundarySkewness
    "vol_ratio_min": 0.01,              # minVolRatio
    "determinant_min": 0.001,           # minDeterminant
    "face_weight_min": 0.05,            # minFaceWeight
    "aspect_ratio_max": 1000.0,         # checkMesh warning threshold
}


def _face_loops(hexes):
    """(face key -> [cell, ...]) and the ordered loop of each face."""
    f2h = defaultdict(list)
    loop = {}
    for hi, c in enumerate(hexes):
        for fc in HF:
            lp = [int(c[i]) for i in fc]
            k = frozenset(lp)
            f2h[k].append(hi)
            loop.setdefault(k, lp)
    return f2h, loop


def face_geometry(P, loop):
    """Centre, area vector and unit normal of a quad, via its two triangles.

    A quad in 3D is generally NOT planar; splitting it consistently is what
    OpenFOAM does too, and it is why `face_flatness` below is a metric in its
    own right rather than an afterthought."""
    Q = P[np.asarray(loop)]
    c = Q.mean(1) if Q.ndim == 3 else Q.mean(0)
    return c


def metrics(P, hexes, verbose=True):
    """All per-cell metrics. Returns a dict of arrays of len(hexes)."""
    P = np.asarray(P, float)
    H = np.asarray(hexes)
    n = len(H)
    f2h, loop = _face_loops(H)

    # cell centres and volumes
    ctr = P[H].mean(1)
    import ovm_io
    vol = np.array([ovm_io._hex_volume(P[c]) for c in H])

    # per-face geometry
    keys = list(f2h)
    L = np.array([loop[k] for k in keys])
    Q = P[L]                                   # (F, 4, 3)
    fc = Q.mean(1)
    # area vector as the sum of the two triangle normals (OpenFOAM-style)
    n1 = np.cross(Q[:, 1] - Q[:, 0], Q[:, 2] - Q[:, 0]) / 2.0
    n2 = np.cross(Q[:, 2] - Q[:, 0], Q[:, 3] - Q[:, 0]) / 2.0
    Sf = n1 + n2
    area = np.linalg.norm(Sf, axis=1)
    nf = Sf / np.maximum(area, 1e-30)[:, None]

    non_ortho = np.zeros(n)
    skew = np.zeros(n)
    vol_ratio = np.ones(n)
    face_weight = np.ones(n)
    flatness = np.ones(n)

    for i, k in enumerate(keys):
        hs = f2h[k]
        # face flatness: the two triangles' normals against the mean
        cosang = float(np.dot(n1[i], n2[i]) /
                       max(np.linalg.norm(n1[i]) * np.linalg.norm(n2[i]), 1e-30))
        fl = max(0.0, cosang)
        if len(hs) == 2:
            a, b = hs
            d = ctr[b] - ctr[a]
            nd = np.linalg.norm(d)
            if nd > 0:
                # non-orthogonality: angle between d and the face normal
                ca = abs(float(np.dot(d / nd, nf[i])))
                ang = np.degrees(np.arccos(np.clip(ca, 0, 1)))
                # skewness: distance from the face centre to where the
                # centre-to-centre line pierces the face, over |d|.
                # The face normal's sign comes from the loop order and is
                # arbitrary, so clamping the denominator with max(x, 1e-30)
                # divided by ~0 whenever it came out negative -- that pushed
                # skewness to 1e28 on half the faces.
                den = float(np.dot(d, nf[i]))
                if abs(den) < 1e-12:
                    continue
                t = float(np.dot(fc[i] - ctr[a], nf[i]) / den)
                pierce = ctr[a] + t * d
                sk = float(np.linalg.norm(pierce - fc[i]) / nd)
                fw = min(abs(t), abs(1 - t))
                for h in (a, b):
                    non_ortho[h] = max(non_ortho[h], ang)
                    skew[h] = max(skew[h], sk)
                    face_weight[h] = min(face_weight[h], fw)
                r = min(abs(vol[a]), abs(vol[b])) / max(abs(vol[a]), abs(vol[b]), 1e-30)
                vol_ratio[a] = min(vol_ratio[a], r)
                vol_ratio[b] = min(vol_ratio[b], r)
        for h in hs:
            flatness[h] = min(flatness[h], fl)

    # aspect ratio: longest over shortest edge of the cell
    E = P[H[:, np.array([e[0] for e in HE])]] - P[H[:, np.array([e[1] for e in HE])]]
    elen = np.linalg.norm(E, axis=2)
    aspect = elen.max(1) / np.maximum(elen.min(1), 1e-30)

    sj = cb.scaled_jacobians(P, H)

    out = {"scaled_jacobian": sj, "non_orthogonality_deg": non_ortho,
           "skewness": skew, "aspect_ratio": aspect, "vol_ratio": vol_ratio,
           "face_weight": face_weight, "face_flatness": flatness,
           "volume": vol, "edge_min": elen.min(1), "edge_max": elen.max(1)}
    if verbose:
        report(out)
    return out


def report(m):
    print(f"{'metric':24s} {'min':>10s} {'median':>10s} {'max':>10s}  limit / violations")
    checks = [
        ("scaled_jacobian", "min", 0.0, "> 0"),
        ("non_orthogonality_deg", "max", LIMITS["non_orthogonality_deg"], "<= 65 deg"),
        ("skewness", "max", LIMITS["skewness_internal"], "<= 4"),
        ("aspect_ratio", "max", LIMITS["aspect_ratio_max"], "<= 1000"),
        ("vol_ratio", "min", LIMITS["vol_ratio_min"], ">= 0.01"),
        ("face_weight", "min", LIMITS["face_weight_min"], ">= 0.05"),
        ("face_flatness", "min", 0.8, ">= 0.8"),
    ]
    for name, side, lim, txt in checks:
        v = m[name]
        bad = int((v <= lim).sum()) if side == "min" else int((v > lim).sum())
        flag = "OK" if bad == 0 else f"{bad} cells"
        print(f"{name:24s} {v.min():10.4f} {np.median(v):10.4f} {v.max():10.4f}"
              f"  {txt:>12s}  {flag}")


def render(path, P, H, m, limits=None, title="mesh quality"):
    """One PNG + interactive HTML per metric, same camera as the showcase.

    Cells that violate the OpenFOAM limit are drawn opaque on top of a
    translucent hull, so a figure shows WHERE a metric fails rather than only
    that it does."""
    import showcase as sc
    path = Path(path)
    hull = sc.hex_grid(P, H).extract_surface()
    for name, v in m.items():
        lim = (limits or {}).get(name)
        g = sc.hex_grid(P, H, **{name: np.asarray(v, float)})
        acts = [dict(mesh=hull, color="#DDDDDD", opacity=0.15, label="mesh"),
                dict(mesh=g.extract_surface(), scalars=name, cmap="viridis",
                     colorscale="Viridis", label=name)]
        if lim is not None:
            side, thr = lim
            bad = np.where(v > thr)[0] if side == "max" else np.where(v <= thr)[0]
            if len(bad):
                acts.append(dict(mesh=sc.hex_grid(P, H[bad]), color="#EE6677",
                                 label=f"{len(bad)} outside {side} {thr}"))
        stem = str(path).replace(".vtk", f"_{name}")
        sc._render(acts, stem + ".png", f"{title}: {name}")
        sc._plotly(acts, stem + ".html", f"{title}: {name}")
    print(f"[mesh_quality] rendered {len(m)} metrics as png + html")


VIOLATION = {"scaled_jacobian": ("min", 0.0),
             "non_orthogonality_deg": ("max", LIMITS["non_orthogonality_deg"]),
             "skewness": ("max", LIMITS["skewness_internal"]),
             "aspect_ratio": ("max", LIMITS["aspect_ratio_max"]),
             "vol_ratio": ("min", LIMITS["vol_ratio_min"]),
             "face_weight": ("min", LIMITS["face_weight_min"]),
             "face_flatness": ("min", 0.8)}


def write_vtk(path, P, H, m, title="mesh quality"):
    import export_vtk as ev
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # one file per metric keeps every field thresholdable on its own scale
    # integers, because export_vtk writes an int CELL_DATA scalar. Each
    # metric gets a factor that keeps its useful range visible.
    FACTOR = {"scaled_jacobian": 1000, "vol_ratio": 1000, "face_weight": 1000,
              "face_flatness": 1000, "skewness": 1000,
              "non_orthogonality_deg": 100, "aspect_ratio": 100,
              "volume": 10 ** 9, "edge_min": 10 ** 6, "edge_max": 10 ** 6}
    for name, v in m.items():
        v = np.nan_to_num(np.asarray(v, float), nan=0.0,
                          posinf=np.iinfo(np.int32).max, neginf=0.0)
        scaled = np.clip(np.round(v * FACTOR.get(name, 1000)),
                         np.iinfo(np.int32).min, np.iinfo(np.int32).max
                         ).astype(np.int64)
        f = str(path).replace(".vtk", f"_{name}.vtk")
        ev.write_vtk(f, P, H, [12] * len(H), scaled, name, f"{title}: {name}")
    return path


if __name__ == "__main__":
    import argparse
    import meshio
    ap = argparse.ArgumentParser()
    ap.add_argument("mesh", nargs="?",
                    default=str(REPO / "output" / "hex3d_algohex" / "deliverable"
                               / "T1_9_blocks_v11_full.vtk"))
    ap.add_argument("--out", default=str(REPO / "output" / "hex3d_algohex"
                                        / "quality" / "quality.vtk"))
    ap.add_argument("--no-render", action="store_true",
                    help="write the VTK fields only, skip png/html")
    a = ap.parse_args()
    mm = meshio.read(a.mesh)
    H = np.vstack([b.data for b in mm.cells if b.type == "hexahedron"])
    print(f"[mesh_quality] {a.mesh}: {len(H)} cells\n")
    m = metrics(mm.points, H)
    write_vtk(a.out, mm.points, H, m)
    if not a.no_render:
        render(a.out, mm.points, H, m, VIOLATION)
    print(f"\n[mesh_quality] fields written next to {a.out}")


# --------------------------------------------------------------------------
# the same metrics on a mixed-cell mesh, for comparing against the source
# --------------------------------------------------------------------------

# faces of each cell type, by vertex count, in VTK ordering and wound
# consistently; the sign is fixed per cell against its own centroid below
CELL_FACES = {
    4: ((0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2)),                # tet
    5: ((0, 3, 2, 1), (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)),  # pyramid
    6: ((0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4),
        (2, 0, 3, 5)),                                              # prism
    8: HF,                                                          # hex
}


def mixed_metrics(P, cells, verbose=True):
    """`metrics` for a mesh of mixed cell types.

    The whole point of this branch is the claim that a structured hex mesh is
    better for the solver than the unstructured source. Testing that claim
    with two different implementations would prove nothing, so this is the one
    used for BOTH sides: the same non-orthogonality, skewness, face-weight and
    volume-ratio formulas as `metrics`, over faces enumerated per cell type.

    `cells` is a list of index lists. Volume comes from the divergence
    theorem over the cell's own faces, which holds for any closed polyhedron;
    face normals are oriented outward against the cell centroid first."""
    P = np.asarray(P, float)
    n = len(cells)
    ctr = np.array([P[np.asarray(c)].mean(0) for c in cells])

    f2c = defaultdict(list)
    floop = {}
    vol = np.zeros(n)
    for ci, c in enumerate(cells):
        fcs = CELL_FACES.get(len(c))
        if fcs is None:
            raise ValueError(f"unsupported cell with {len(c)} vertices")
        v = 0.0
        for fc in fcs:
            lp = [int(c[i]) for i in fc]
            Q = P[np.asarray(lp)]
            fcen = Q.mean(0)
            Sf = np.zeros(3)
            for t in range(1, len(lp) - 1):
                Sf = Sf + np.cross(Q[t] - Q[0], Q[t + 1] - Q[0]) / 2.0
            if np.dot(Sf, fcen - ctr[ci]) < 0:      # orient outward
                Sf = -Sf
            v += float(np.dot(fcen, Sf)) / 3.0
            k = frozenset(lp)
            f2c[k].append(ci)
            floop.setdefault(k, lp)
        vol[ci] = v

    non_ortho = np.zeros(n)
    skew = np.zeros(n)
    vol_ratio = np.ones(n)
    face_weight = np.ones(n)
    for k, hs in f2c.items():
        if len(hs) != 2:
            continue
        a, b = hs
        lp = floop[k]
        Q = P[np.asarray(lp)]
        fcen = Q.mean(0)
        Sf = np.zeros(3)
        for t in range(1, len(lp) - 1):
            Sf = Sf + np.cross(Q[t] - Q[0], Q[t + 1] - Q[0]) / 2.0
        na = np.linalg.norm(Sf)
        if na < 1e-30:
            continue
        nf = Sf / na
        d = ctr[b] - ctr[a]
        nd = float(np.linalg.norm(d))
        if nd <= 0:
            continue
        ang = np.degrees(np.arccos(np.clip(abs(float(np.dot(d / nd, nf))), 0, 1)))
        den = float(np.dot(d, nf))
        if abs(den) < 1e-12:
            continue
        t = float(np.dot(fcen - ctr[a], nf) / den)
        sk = float(np.linalg.norm((ctr[a] + t * d) - fcen) / nd)
        fw = min(abs(t), abs(1 - t))
        for h in (a, b):
            non_ortho[h] = max(non_ortho[h], ang)
            skew[h] = max(skew[h], sk)
            face_weight[h] = min(face_weight[h], fw)
        r = min(abs(vol[a]), abs(vol[b])) / max(abs(vol[a]), abs(vol[b]), 1e-30)
        vol_ratio[a] = min(vol_ratio[a], r)
        vol_ratio[b] = min(vol_ratio[b], r)

    aspect = np.zeros(n)
    for ci, c in enumerate(cells):
        Q = P[np.asarray(c)]
        dd = np.linalg.norm(Q[:, None, :] - Q[None, :, :], axis=-1)
        pos = dd[dd > 1e-30]
        aspect[ci] = pos.max() / pos.min() if len(pos) else 1.0

    out = {"non_orthogonality_deg": non_ortho, "skewness": skew,
           "vol_ratio": vol_ratio, "face_weight": face_weight,
           "aspect_ratio": aspect, "volume": vol}
    if verbose:
        print(f"[mesh_quality] mixed mesh: {n} cells, "
              f"{sum(1 for v in f2c.values() if len(v) == 2)} internal faces")
    return out
