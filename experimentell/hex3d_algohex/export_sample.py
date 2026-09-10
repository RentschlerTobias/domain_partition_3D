"""Stage 8 (T9): one neutral geometry file per sample.

The seam of decision C: this repository writes neutral geometry, meshtron owns
the ML format. Everything that will be tuned during training -- coordinate
frame, ordering, point-cloud density, normalisation, quantization -- belongs on
the cheap side of the seam, because tuning it must never re-run AlgoHex. So
this file is deliberately over-complete: anything missing from it costs a
13-55 minute re-run to add.

Format is `.npz`, not `.pt`. Decision C keeps this repository free of torch,
and npz needs numpy alone on both banks. Dicts (`params`, `quality`,
`provenance`) ride along as JSON strings in the same archive, so a sample is
one file.

Run:
    PY=/root/repos/duty/quadmesh/.venv/bin/python
    $PY experimentell/hex3d_algohex/export_sample.py \
        output/hex3d_algohex/deliverable/T1_9_blocks_v11.vtk \
        --tet data/T1_9/T1_9_tet_v5.vtk --out /tmp/v11.npz --check
"""

import json
import subprocess
from collections import defaultdict
from pathlib import Path

import numpy as np

import block_edges as be
import tfi

REPO = Path(__file__).resolve().parent.parent.parent

# VTK hex corner order, as (i, j, k) in {0, 1}^3 -- the same convention
# `tfi.block_cells` uses, so a block written here refills without a relabel.
CORNERS = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
           (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
# the six faces of a hex in that order, as corner slots
FACES = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
         (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]


# --------------------------------------------------------------------------
# the block complex, reduced to its corners
# --------------------------------------------------------------------------

def block_corners(lat):
    """(corner vertex ids per block, in VTK order) for every lattice block."""
    out = {}
    for r, (_dims, vert) in lat.items():
        idx = []
        for c in CORNERS:
            sl = tuple(0 if s == 0 else -1 for s in c)
            idx.append(int(vert[sl]))
        out[int(r)] = idx
    return out


def quad_shell(blocks):
    """The block faces used by exactly one block.

    Block-level quads, not fine-mesh quads: four corners each. A face used
    twice is interior. Note this is a topological test on corner ids, so a
    face whose two sides did not weld shows up as TWO one-sided faces -- the
    same signature T12 is chasing, and the reason `quality` carries the count.
    """
    seen = defaultdict(list)
    for r, c in blocks.items():
        for f in FACES:
            quad = tuple(c[i] for i in f)
            seen[frozenset(quad)].append((r, quad))
    shell, interior = [], 0
    for _k, v in seen.items():
        if len(v) == 1:
            shell.append(v[0][1])
        else:
            interior += 1
    return shell, interior


def axis_class_of_edge(lat, classes, owners):
    """Direction-class id per undirected edge chain.

    An edge runs along one axis of each block that owns it, and every such
    (block, axis) is in the same class by construction -- that is what
    `tfi.direction_classes` computes. So the class is well defined per edge,
    and disagreement means the classes are wrong, which is worth knowing.
    """
    cls_of = {}
    for ci, members in enumerate(classes):
        for key in members:
            cls_of[(int(key[0]), int(key[1]))] = ci
    out, bad = {}, 0
    for chain, own in owners.items():
        ids = {cls_of.get((r, ax), -1) for r, ax in own}
        if len(ids) > 1:
            bad += 1
        out[chain] = min(ids)
    return out, bad


# --------------------------------------------------------------------------
# the surface, stored as a triangulation rather than as a point cloud
# --------------------------------------------------------------------------

def surface_arrays(tet_vtk):
    """(points, triangles, surface id per triangle) from the AlgoHex input.

    **Deviation from decision C's field list, and it is deliberate.** C asks
    for `surface_points [N,3]` and `surface_label [N]`, one label per point.
    That cannot be written honestly: a vertex on a ring between two surfaces
    belongs to both, and picking a winner would bake a tie-break into the
    EXPENSIVE side of the seam -- the one thing decision C exists to prevent.
    The triangulation carries the labels unambiguously, per triangle, and any
    point cloud at any density with any tie-break follows from it in seconds
    at load time, which is where decision E was deferred to anyway.
    """
    import clean_blocks as cb
    P, tri, tid = cb.read_input_surface(str(tet_vtk))
    return np.asarray(P, float), np.asarray(tri, np.int64), np.asarray(tid, int)


# --------------------------------------------------------------------------

def _git_sha():
    try:
        return subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"],
                              capture_output=True, text=True,
                              check=True).stdout.strip()
    except Exception:
        return "unknown"


def build(blocks_vtk, tet_vtk=None, params_json=None, status_json=None,
          target_h=None, n=None, collapse_rounds=None, verbose=True):
    """Everything the sample file holds, as plain arrays and dicts."""
    blocks_vtk = Path(blocks_vtk)
    P, H, B, f2h = tfi.load_blocks(str(blocks_vtk))
    lat, missing = tfi.lattices(P, H, B, f2h, verbose=verbose)
    classes = tfi.direction_classes(lat, f2h, H, B, verbose=verbose)

    corners = block_corners(lat)
    owners = be.edge_chains(lat)
    chains = sorted(owners)
    recs = be.fit_edges(P, chains)
    fit = be.summary(recs)
    cls_of_edge, cls_disagree = axis_class_of_edge(lat, classes, owners)

    # global vertex id -> compact corner index. Corners only: the fine mesh is
    # not part of the sample, the block complex is (decision A).
    used = sorted({v for c in corners.values() for v in c}
                  | {chain[0] for chain in chains}
                  | {chain[-1] for chain in chains})
    remap = {v: i for i, v in enumerate(used)}
    vertices = P[used]

    blocks = np.array([corners[r] for r in sorted(corners)], np.int64)
    blocks = np.vectorize(remap.get)(blocks)
    shell, n_interior = quad_shell({r: [remap[v] for v in c]
                                    for r, c in corners.items()})
    quad_faces = (np.array(shell, np.int64) if shell
                  else np.zeros((0, 4), np.int64))

    # DIRECTED edges: both traversals, per decision D2. The reversed direction
    # reverses the control points, so a consumer never has to know which way a
    # stored edge was fitted.
    edges, edge_ctrl, poly, poly_off, dir_class, edge_len = [], [], [], [0], [], []
    for r in recs:
        chain = tuple(r["chain"])
        pts = P[list(chain)]
        B1, B2 = (np.asarray(c, float) for c in r["ctrl"])
        ci = cls_of_edge[chain]
        for fwd in (True, False):
            q = pts if fwd else pts[::-1]
            edges.append([remap[chain[0]], remap[chain[-1]]] if fwd
                         else [remap[chain[-1]], remap[chain[0]]])
            edge_ctrl.append([B1, B2] if fwd else [B2, B1])
            poly.append(q)
            poly_off.append(poly_off[-1] + len(q))
            dir_class.append(ci)
            edge_len.append(len(q))

    counts = None
    if target_h is not None:
        # {class index: count}, so it lines up with dir_class by construction
        solved = tfi.solve_block_divisions(lat, classes, P, target_h,
                                           verbose=verbose)
        counts = np.array([solved[c] for c in range(len(classes))], np.int64)

    quality = {
        "blocks": len(lat),
        "blocks_without_lattice": len(missing),
        "quad_shell_faces": len(shell),
        "interior_faces": n_interior,
        "direction_classes": len(classes),
        "class_disagreements": cls_disagree,
        "edges_undirected": len(recs),
        "fit_residual_median": fit["err_rel_median"],
        "fit_residual_p95": fit["err_rel_p95"],
        "fit_residual_max": fit["err_rel_max"],
        "fit_mode_cubic": fit["fit_mode_cubic"],
        "fit_mode_quadratic": fit["fit_mode_quadratic"],
        "fit_mode_chord": fit["fit_mode_chord"],
        "chord_min": fit["chord_min"],
        "chord_median": fit["chord_median"],
        "degenerate_edges": fit["degenerate_edges"],
        "inflection_edges": fit["edges_with_inflection"],
        "planarity_p95": fit["planarity_p95"],
        # the two fields a load-time filter actually wants (decision G stores
        # quality rather than gating on it): arc/chord predicts which edges a
        # single cubic cannot represent, and edges_over_5pct counts them
        "arc_over_chord_p95": fit["arc_over_chord_p95"],
        "arc_over_chord_max": fit["arc_over_chord_max"],
        "edges_over_5pct": fit["edges_over_5pct"],
    }
    if status_json and Path(status_json).exists():
        st = json.loads(Path(status_json).read_text())
        quality["pipeline_status"] = st
        # `status.json` can disagree with the file it describes: on cand_001 it
        # says 22 blocks while cand_001_blocks.vtk carries 21 contiguous block
        # ids, all of them lattices. The VTK is the artifact, so quality
        # ["blocks"] is authoritative -- but both numbers are in here now, and
        # a filter that reads the wrong one fails silently. This flag makes the
        # disagreement visible instead.
        if "blocks" in st and int(st["blocks"]) != len(lat):
            quality["status_blocks_disagree"] = [int(st["blocks"]), len(lat)]
            if verbose:
                print(f"[export_sample] status.json says {st['blocks']} blocks,"
                      f" the VTK has {len(lat)} -- flagged, VTK wins")

    params = {}
    if params_json and Path(params_json).exists():
        params = json.loads(Path(params_json).read_text())

    provenance = {
        "blocks_vtk": str(blocks_vtk),
        "tet_vtk": str(tet_vtk) if tet_vtk else None,
        "name": blocks_vtk.stem,
        "n": n,
        "collapse_rounds": collapse_rounds,
        "target_h": target_h,
        "git_sha": _git_sha(),
        "exporter": "export_sample.py",
    }

    out = {
        "vertices": vertices,
        "blocks": blocks,
        "quad_faces": quad_faces,
        "edges": np.asarray(edges, np.int64),
        "edge_ctrl": np.asarray(edge_ctrl, float),
        "edge_polyline": np.concatenate(poly) if poly
        else np.zeros((0, 3), float),
        "edge_polyline_offset": np.asarray(poly_off, np.int64),
        "dir_class": np.asarray(dir_class, np.int64),
        "params": json.dumps(params),
        "quality": json.dumps(quality),
        "provenance": json.dumps(provenance),
    }
    if counts is not None:
        out["dir_class_count"] = counts
    if tet_vtk:
        sp, st, sl = surface_arrays(tet_vtk)
        out["surface_points"] = sp
        out["surface_tris"] = st
        out["surface_tri_label"] = sl
    return out, quality


def write(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **data)
    return path


def load(path):
    """The sample back as a dict, with the JSON fields decoded."""
    z = np.load(path, allow_pickle=False)
    out = {k: z[k] for k in z.files}
    for k in ("params", "quality", "provenance"):
        if k in out:
            out[k] = json.loads(str(out[k]))
    return out


def polyline(sample, i):
    """Edge `i`'s raw polyline out of the flat array."""
    o = sample["edge_polyline_offset"]
    return sample["edge_polyline"][o[i]:o[i + 1]]


# --------------------------------------------------------------------------
# the gate
# --------------------------------------------------------------------------

def roundtrip(sample_path, blocks_vtk, verbose=True):
    """Rebuild the block complex from the sample ALONE, compare to the source.

    What is compared, and what deliberately is not. The sample stores the
    block complex, not the fine hex mesh -- that is decision A, and refilling
    it back to cells is T7's job, not a round-trip. So the test is: same
    blocks, same corner coordinates, same 8-corner sets, same undirected edge
    set, and every stored control point reproducing its own polyline to the
    residual the fit reported. Topology exact, geometry to floating point.
    """
    s = load(sample_path)
    P, H, B, f2h = tfi.load_blocks(str(blocks_vtk))
    lat, _missing = tfi.lattices(P, H, B, f2h, verbose=False)
    corners = block_corners(lat)
    owners = be.edge_chains(lat)

    fails = []

    # blocks, as sets of corner coordinates rounded to the VTK file's own
    # precision -- comparing ids would only test the remap against itself
    def kset(coords):
        return frozenset(tuple(np.round(c, 9)) for c in coords)

    src_blocks = {kset(P[c]) for c in corners.values()}
    got_blocks = {kset(s["vertices"][b]) for b in s["blocks"]}
    if src_blocks != got_blocks:
        fails.append(f"blocks differ: {len(src_blocks)} source, "
                     f"{len(got_blocks)} sample, "
                     f"{len(src_blocks ^ got_blocks)} not shared")

    # edges: the directed set must be exactly each undirected edge twice
    src_edges = {kset([P[ch[0]], P[ch[-1]]]) for ch in owners}
    got_edges = {kset(s["vertices"][e]) for e in s["edges"]}
    if src_edges != got_edges:
        fails.append(f"edge endpoints differ: {len(src_edges)} source, "
                     f"{len(got_edges)} sample")
    if len(s["edges"]) != 2 * len(owners):
        fails.append(f"expected {2 * len(owners)} directed edges, "
                     f"got {len(s['edges'])}")

    # geometry: the stored control points must still fit their stored polyline
    worst = 0.0
    for i in range(len(s["edges"])):
        pts = polyline(s, i)
        B1, B2 = s["edge_ctrl"][i]
        curve = be.bezier_points(pts[0], B1, B2, pts[-1], 200)
        chord = float(np.linalg.norm(pts[-1] - pts[0]))
        if chord > 0:
            worst = max(worst, be.max_dist_to_curve(pts, curve) / chord)
    q = s["quality"]
    if worst > q["fit_residual_max"] + 1e-9:
        fails.append(f"control points do not reproduce their polyline: "
                     f"{worst:.3%} against a reported max of "
                     f"{q['fit_residual_max']:.3%}")

    # the reversed direction must be the same curve, not a different one
    rev = 0
    for i in range(0, len(s["edges"]), 2):
        a, b = s["edge_ctrl"][i], s["edge_ctrl"][i + 1]
        if not (np.allclose(a[0], b[1]) and np.allclose(a[1], b[0])):
            rev += 1
    if rev:
        fails.append(f"{rev} directed pairs are not reverses of each other")

    if verbose:
        print(f"\n[export_sample] round-trip against {Path(blocks_vtk).name}")
        print(f"  blocks           {len(got_blocks)} / {len(src_blocks)}"
              f"   {'exact' if src_blocks == got_blocks else 'MISMATCH'}")
        print(f"  edge endpoints   {len(got_edges)} / {len(src_edges)}"
              f"   {'exact' if src_edges == got_edges else 'MISMATCH'}")
        print(f"  directed edges   {len(s['edges'])} "
              f"(= 2 x {len(owners)})")
        print(f"  curve residual   worst {worst * 100:.3f} % "
              f"<= reported {q['fit_residual_max'] * 100:.3f} %")
        print(f"  verdict          "
              f"{'PASS' if not fails else 'FAIL: ' + '; '.join(fails)}")
    return not fails, fails


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("blocks")
    ap.add_argument("--tet", default=None,
                    help="AlgoHex input VTK, for the labelled surface")
    ap.add_argument("--params", default=None, help="dtOO params.json")
    ap.add_argument("--status", default=None, help="pipeline status.json")
    ap.add_argument("--target-h", type=float, default=None,
                    help="solve the conforming division counts at this cell "
                         "size and store them per direction class")
    ap.add_argument("--n", type=int, default=None,
                    help="the -n AlgoHex ran with, for provenance")
    ap.add_argument("--collapse-rounds", type=int, default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--check", action="store_true",
                    help="round-trip the written file against the source")
    a = ap.parse_args()

    out = a.out or str(Path(a.blocks).with_suffix("")) + "_sample.npz"
    data, quality = build(a.blocks, tet_vtk=a.tet, params_json=a.params,
                          status_json=a.status, target_h=a.target_h,
                          n=a.n, collapse_rounds=a.collapse_rounds)
    p = write(out, data)
    size = p.stat().st_size
    print(f"\n[export_sample] wrote {p}  ({size / 1e6:.2f} MB)")
    print(f"  vertices {len(data['vertices'])}, blocks {len(data['blocks'])}, "
          f"shell {len(data['quad_faces'])}, "
          f"directed edges {len(data['edges'])}")
    print(f"  fit residual median {quality['fit_residual_median'] * 100:.3f} %,"
          f" modes {quality['fit_mode_cubic']}/"
          f"{quality['fit_mode_quadratic']}/{quality['fit_mode_chord']}")
    if "surface_tris" in data:
        print(f"  surface {len(data['surface_points'])} points, "
              f"{len(data['surface_tris'])} tris, "
              f"{len(np.unique(data['surface_tri_label']))} labels")
    if a.check:
        ok, _f = roundtrip(p, a.blocks)
        raise SystemExit(0 if ok else 3)
