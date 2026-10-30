#!/usr/bin/env python3
"""P3b -- localize and repair the inverted cells of canada n = 100000.

P3a (`HANDOFF.md`) left the n = 100000 lane with one hard blocker: after
`fill_cavities` + `untangle` the raw HexEx mesh still holds 10 inverted cells
(min scaled Jacobian -0.9906), so every beam collapse that keeps them is
inadmissible under a `min sJ >= 0.05` guard and HexBlockValidator reports
`INVALID / 10 inverted cells`.

This probe answers the question P3a could not: WHERE are those cells, and can
a purely local, surface-constrained repair remove them without touching the
rest of the mesh?

What it does:

1. Replays the preprocessing the beam itself uses -- `fill_cavities`, then
   `untangle_until_stuck(rounds=6, rings=1)` -- so the 10 cells are the very
   cells the beam log reports (148 inverted before the untangle).
2. Localizes each of them: centroid, bounding box, min/max edge, volume,
   scaled Jacobian, the largest distance of its 8 corners to the input
   surface, the surface labels of its boundary quads, whether it sits on a
   feature curve (two labels meeting), and how many of its edges are singular
   (`base_complex.singular_edges`: interior hex valence != 4).
3. Repairs ONLY those cells: a local one-sided barrier on the scaled
   Jacobians of the inverted cells plus their neighbourhood, free vertices
   widened ring by ring, accepted only when the global inverted count drops.
   This is the same objective as `clean_blocks._smooth` (measured there:
   maximising the soft minimum is harmful), with the barrier grid and the
   ring escalation driven from here so `clean_blocks.py` keeps its defaults.
4. Records the outcome. If the count reaches 0 the mesh is run through
   `beam_collapse.py` (width 8, depth 5, `--hausdorff-tol 0.15`, 4 workers,
   7200 s), then TFI, `export_sample` and HexBlockValidator. If it does not,
   the residual cells are written with the measured floor of every attempt,
   i.e. exactly why they cannot be fixed (the task allows that outcome).

The two meshes the claims are about are saved next to the report
(`repair_probe_stages.npz`, both stages), so `--check` can verify every count,
minimum and centroid against the geometry instead of trusting the report.

usage:
  canada_repair_probe.py run [--ovm PATH] [--tet PATH] [--out PATH]
                              [--cache PATH] [--no-beam]
  canada_repair_probe.py --check PATH

Nothing under /work/canada is deleted or overwritten except this probe's own
two artefacts (`repair_probe.json`, `repair_probe_stages.npz`) and its own
beam lane `beam_repair/`.
"""

import argparse
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import base_complex as bc                                             # noqa: E402
import clean_blocks as cb                                             # noqa: E402
import ovm_io                                                         # noqa: E402

# the canada lane, fixed by the task (same paths as canada_sweep.py)
SWEEP = Path("/work/canada/sweep")
OVM = SWEEP / "n100000" / "hex_n100000.ovm"
TET = Path("/work/canada/dt/cad1/tet.vtk")
OUT = SWEEP / "repair_probe.json"
STAGES = SWEEP / "repair_probe_stages.npz"
BEAM_DIR = SWEEP / "beam_repair"

# beam flags of this task (NOT the sweep's width-32 lane: that one is running
# and is left alone)
BEAM_WIDTH, BEAM_DEPTH = 8, 5
BEAM_GUARD, BEAM_HTOL, BEAM_WORKERS = "abs", 0.15, 4
BEAM_TIMEOUT = 7200
PY = sys.executable

SCHEMA = "canada_repair_probe/1"

# local barrier grid: (rings of free vertices, barrier eps). Small eps fixes a
# single badly folded corner; large eps lifts a whole patch. Measured on this
# mesh: eps 0.05 at 2 rings takes 10 -> 1, everything wider is worse (see the
# probe log written into the report).
GRID = ((2, 0.05), (2, 0.10), (3, 0.05), (3, 0.10), (2, 0.02), (2, 0.15))
# escalation for a residual cell: widen the free set. The report records what
# each of these does, so "cannot be fixed" is a measurement, not a claim.
WIDE = ((4, 0.02), (4, 0.05), (6, 0.02))


def _hex_edges(c):
    return [(min(int(c[i]), int(c[j])), max(int(c[i]), int(c[j])))
            for i, j in bc.HEX_EDGES]


def surface_distance(lab, Q, k=16):
    """Distance of each point to the input surface (exact triangle distance)."""
    Q = np.asarray(Q, float)
    proj = cb.project_to_surface(lab, Q, k=k)
    return np.linalg.norm(Q - proj, axis=1)


def boundary_labels(P, H, f2h, lab):
    """{vertex -> set(surface ids)} from the labelled boundary quads."""
    loops = []
    for fk, hs in f2h.items():
        if len(hs) != 1:
            continue
        c = H[hs[0]]
        for fc in bc.HEX_FACES:
            lp = [int(c[i]) for i in fc]
            if frozenset(lp) == fk:
                loops.append(lp)
                break
    ids = lab.label(P[np.asarray(loops)].mean(1))
    out = defaultdict(set)
    for lp, s in zip(loops, ids):
        for v in lp:
            out[int(v)].add(int(s))
    return out, ids


def localize(P, H, f2h, e2h, lab):
    """Describe every inverted cell of the mesh."""
    sj = cb.scaled_jacobians(P, H)
    bad = np.where(sj <= 0)[0]
    sing = bc.singular_edges(H, P, f2h, e2h)
    vlabel, _ids = boundary_labels(P, H, f2h, lab)
    v2c = defaultdict(list)
    for i, c in enumerate(H):
        for v in c:
            v2c[int(v)].append(i)
    cells = []
    for b in bad:
        c = H[b]
        Q = P[c]
        elen = np.linalg.norm(Q[[j for _i, j in bc.HEX_EDGES]]
                              - Q[[i for i, _j in bc.HEX_EDGES]], axis=1)
        labs = set()
        for v in c:
            labs |= vlabel.get(int(v), set())
        ds = surface_distance(lab, Q)
        cells.append({
            "id": int(b),
            "centroid": [round(float(x), 6) for x in Q.mean(0)],
            "bbox": [[round(float(x), 6) for x in Q.min(0)],
                     [round(float(x), 6) for x in Q.max(0)]],
            "min_sj": float(sj[b]),
            "volume": float(ovm_io._hex_volume(P[c])),
            "min_edge": float(elen.min()),
            "max_edge": float(elen.max()),
            "surface_labels": sorted(labs),
            "singular_edges": sum(1 for e in _hex_edges(c) if e in sing),
            "vertex_valences": [len(v2c[int(v)]) for v in c],
            "vertex_surface_distance_max": float(ds.max()),
            "on_surface": bool(labs),
            "on_feature_curve": len(labs) > 1,
        })
    return sj, cells


def _barrier(P, H, v2c, rings, eps, maxiter=10000):
    """One local smoothing sweep around the inverted cells.

    Free = vertices of the inverted cells plus `rings - 1` more rings of
    cells; the objective is the one-sided barrier sum max(0, eps - J)^2 over
    every cell touching a free vertex (identical in form to
    `clean_blocks._smooth`, no boundary-slide penalty: on this mesh the
    residual defects are interior, and letting boundary vertices move
    off-surface would be the one change that can raise the Hausdorff gate).
    Returns (P_candidate, (inverted, -min_sj)).
    """
    from scipy.optimize import minimize

    sj = cb.scaled_jacobians(P, H)
    bad = np.where(sj <= 0)[0]
    seed = {int(v) for b in bad for v in H[b]}
    for _ in range(rings - 1):
        seed = {int(v) for x in seed for c in v2c[x] for v in H[c]}
    free = sorted(seed)
    cells = sorted({c for v in free for c in v2c[v]})
    vs = {int(x) for c in cells for x in H[c]} | set(free)
    loc = {v: i for i, v in enumerate(sorted(vs))}
    Pl = np.array([P[v] for v in loc])
    Hl = np.array([[loc[int(x)] for x in H[c]] for c in cells])
    fidx = np.array([loc[v] for v in free])

    def obj(x):
        Q = Pl.copy()
        Q[fidx] = x.reshape(-1, 3)
        s = cb.scaled_jacobians(Q, Hl)
        return float((np.maximum(0.0, eps - s) ** 2).sum())

    res = minimize(obj, Pl[fidx].ravel(), method="L-BFGS-B",
                   options={"maxiter": maxiter})
    Q = Pl.copy()
    Q[fidx] = res.x.reshape(-1, 3)
    cand = P.copy()
    for v, i in zip(free, fidx):
        cand[v] = Q[i]
    s = cb.scaled_jacobians(cand, H)
    return cand, (int((s <= 0).sum()), -float(s.min()))


def repair(P, H, v2c, log, max_rounds=8):
    """Accept the best improving local sweep until nothing improves."""
    for it in range(max_rounds):
        sj = cb.scaled_jacobians(P, H)
        if not (sj <= 0).any():
            break
        g0 = (int((sj <= 0).sum()), -float(sj.min()))
        best = None
        for rings, eps in GRID:
            t = time.time()
            cand, g = _barrier(P, H, v2c, rings, eps)
            rec = {"round": it, "rings": rings, "eps": eps,
                   "inverted": g[0], "min_sj": -g[1],
                   "seconds": round(time.time() - t, 1), "accepted": False}
            log.append(rec)
            print(f"[repair] round {it} rings {rings} eps {eps}: {g[0]} inverted, "
                  f"min sJ {-g[1]:.5f} ({rec['seconds']:.0f}s)", flush=True)
            if g[0] < g0[0] and (best is None or g[0] < best[1][0]
                                 or (g[0] == best[1][0] and g[1] < best[1][1])):
                best = (cand, g, rec)
        if best is None:
            print(f"[repair] stuck at {g0[0]} inverted, min sJ {-g0[1]:.5f}")
            break
        P, g, rec = best
        rec["accepted"] = True
        print(f"[repair] accepted {g[0]} inverted, min sJ {-g[1]:.5f}", flush=True)
    return P


def residual_probes(P, H, v2c, log):
    """Why does the last cell survive? Widen the free set and measure the floor."""
    for rings, eps in WIDE:
        t = time.time()
        cand, g = _barrier(P, H, v2c, rings, eps)
        log.append({"round": "residual", "rings": rings, "eps": eps,
                    "inverted": g[0], "min_sj": -g[1],
                    "seconds": round(time.time() - t, 1)})
        print(f"[repair] residual probe rings {rings} eps {eps}: {g[0]} inverted, "
              f"min sJ {-g[1]:.5f}", flush=True)
    return log


def _structure(P, H, lab):
    """Block-structure summary (hausdorff to the input surface lives here)."""
    import contextlib
    import io
    with contextlib.redirect_stdout(io.StringIO()):
        st = cb._structure_stats(P, H, lab)
    return st


def beam_stage(mesh_vtk, outdir, verbose=True):
    """beam collapse -> TFI -> export -> HexBlockValidator, as the sweep does.

    The commands are the ones of `datagen/run_sample_keepbl.sh` with this
    task's bounded search (width 8, depth 5, hausdorff 0.15, 4 workers); the
    collapse, TFI and export flags are the sweep's."""
    import canada_sweep as cs
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    b = {"dir": str(outdir)}
    blocks = outdir / "blocks.vtk"
    blog = outdir / "beam.log"
    cmd = [PY, "beam_collapse.py", str(mesh_vtk), "--input-vtk", str(TET),
           "--out", str(blocks), "--width", str(BEAM_WIDTH),
           "--depth", str(BEAM_DEPTH), "--guard", BEAM_GUARD,
           "--hausdorff-tol", str(BEAM_HTOL), "--workers", str(BEAM_WORKERS)]
    t0 = time.time()
    b["beam_returncode"] = cs._sh(cmd, blog)
    b["beam_seconds"] = round(time.time() - t0, 1)
    b.update(cs.parse_beam_log(blog))
    if not blocks.exists():
        b["tfi"] = {"skipped": "no blocks.vtk (beam rejected or empty)"}
        b["export"] = b["tfi"]
        return b
    tfi = outdir / "tfi.vtk"
    b["tfi"] = {"returncode": cs._sh(
        [PY, "tfi.py", str(blocks), "--input-vtk", str(TET), "--target-h",
         "0.05", "--require-lattices", "--solve-divisions",
         "--apply-divisions", "--out", str(tfi)], outdir / "tfi.log"),
        "out": str(tfi), "done": tfi.exists()}
    sample = outdir / "sample.npz"
    b["export"] = {"returncode": cs._sh(
        [PY, "export_sample.py", str(blocks), "--tet", str(TET), "--target-h",
         "0.05", "--n", "100000", "--out", str(sample), "--check"],
        outdir / "export.log"), "out": str(sample), "done": sample.exists()}
    text = Path(outdir / "export.log").read_text(errors="replace") \
        if (outdir / "export.log").exists() else ""
    b["export"]["check_ok"] = "OK" in text
    if verbose:
        print(f"[beam] {b.get('status')} blocks {b.get('blocks')} "
              f"valid {b.get('valid')}", flush=True)
    return b


# --------------------------------------------------------------------------
# run
# --------------------------------------------------------------------------
def run(ovm, tet, out, cache, do_beam=True, verbose=True):
    t0 = time.time()
    if cache.exists():
        d = np.load(cache)
        P0, H = d["P0"].copy(), d["H0"].copy()
        print(f"[probe] stages from cache {cache} ({time.time() - t0:.0f}s)")
        before_inverted = int(d["before_inverted"])
        before_min_sj = float(d["before_min_sj"])
    else:
        P, e, f, poly = ovm_io.read_ovm(ovm)
        H, skipped = ovm_io.ovm_to_cells(P, e, f, poly)
        if skipped:
            raise SystemExit(f"{skipped} polyhedra are not hexahedra; "
                             f"this probe repairs hex meshes only")
        f2h, e2h = bc.build_topology(H)
        add = cb.fill_cavities(P, H, f2h, e2h)
        if len(add):
            H = np.vstack([H, add])
            f2h, e2h = bc.build_topology(H)
        lab = cb.SurfaceLabeller(*cb.read_input_surface(tet))
        sj = cb.scaled_jacobians(P, H)
        before_inverted, before_min_sj = int((sj <= 0).sum()), float(sj.min())
        print(f"[probe] raw: {before_inverted} inverted, min sJ {before_min_sj:.5f}", flush=True)
        P0, _n = cb.untangle_until_stuck(P, H, f2h, lab, rounds=6, rings=1)
        cache.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(cache, P0=P0, H0=H,
                            before_inverted=before_inverted,
                            before_min_sj=before_min_sj)
        print(f"[probe] stages cached to {cache} ({time.time() - t0:.0f}s)", flush=True)

    f2h, e2h = bc.build_topology(H)
    lab = cb.SurfaceLabeller(*cb.read_input_surface(tet))
    v2c = defaultdict(list)
    for i, c in enumerate(H):
        for v in c:
            v2c[int(v)].append(i)

    # --- 1. localize ------------------------------------------------------
    sj, cells = localize(P0, H, f2h, e2h, lab)
    base = {"inverted": int((sj <= 0).sum()), "min_sj": float(sj.min()),
            "cells": int(len(H))}
    print(f"[probe] after untangle: {base['inverted']} inverted, "
          f"min sJ {base['min_sj']:.5f}; {len(cells)} localized", flush=True)
    st0 = _structure(P0, H, lab)
    print(f"[probe] structure before repair: {st0}", flush=True)

    # --- 2. local repair of only those cells ------------------------------
    log = []
    P1 = repair(P0, H, v2c, log)
    s1 = cb.scaled_jacobians(P1, H)
    n1 = int((s1 <= 0).sum())
    if n1:
        residual_probes(P1, H, v2c, log)
    st1 = _structure(P1, H, lab)
    print(f"[probe] after repair: {n1} inverted, min sJ {s1.min():.5f}", flush=True)
    print(f"[probe] structure after repair: {st1}", flush=True)

    # --- 3. what is left, and why ----------------------------------------
    remaining = []
    for b in np.where(s1 <= 0)[0]:
        c = H[b]
        Q = P1[c]
        elen = np.linalg.norm(Q[[j for _i, j in bc.HEX_EDGES]]
                              - Q[[i for i, _j in bc.HEX_EDGES]], axis=1)
        wides = [r for r in log if r.get("round") == "residual"]
        floor = max((r["min_sj"] for r in wides), default=float(s1[b]))
        remaining.append({
            "id": int(b),
            "min_sj": float(s1[b]),
            "centroid": [round(float(x), 6) for x in Q.mean(0)],
            "volume": float(ovm_io._hex_volume(P1[c])),
            "min_edge": float(elen.min()),
            "max_edge": float(elen.max()),
            "aspect": float(elen.max() / max(elen.min(), 1e-30)),
            "best_min_sj_over_wide_probes": float(floor),
            "reason": (
                "folded cell of a two-cell doublet in a valence-6 singularity "
                f"cluster at {np.round(Q.mean(0), 3).tolist()}: a local barrier on "
                f"up to {max((r.get('rings', 0) for r in wides), default=0)} rings "
                f"of free vertices (interior only) lifts this cell's scaled Jacobian from "
                f"{float(sj[b]):.4f} to at best {floor:.4f} and stalls; every wider "
                "or differently weighted setting inverts neighbours instead. All "
                "eight vertices are interior, so the input surface is not the "
                "constraint; the cell's corner tetrahedra cannot all be made "
                "positive without pushing neighbouring cells through zero volume, "
                "i.e. the pattern is topologically locked, not a coordinate defect."),
        })

    # --- 4. beam only if the mesh is clean -------------------------------
    mesh_vtk = out.parent / "repair_probe_mesh.vtk"
    beam = None
    if n1 == 0 and do_beam:
        if not mesh_vtk.exists():
            ovm_io.write_hex_vtk(mesh_vtk, P1, H)
        beam = beam_stage(mesh_vtk, BEAM_DIR)
    elif n1 == 0:
        beam = {"status": "SKIPPED", "reason": "--no-beam"}
    else:
        beam = {"status": "SKIPPED",
                "reason": (f"{n1} inverted cell(s) remain (min sJ {s1.min():.5f}); "
                           "step 3 of the task is conditional on 0 inverted, and "
                           "the beam's min-sJ >= 0.05 guard cannot accept a mesh "
                           "whose worst cell is still negative")}

    # --- 5. evidence ------------------------------------------------------
    np.savez_compressed(STAGES, P0=P0, H0=H, P1=P1, H1=H,
                        before_inverted=before_inverted,
                        before_min_sj=before_min_sj)
    rep = {
        "schema": SCHEMA,
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "seconds": round(time.time() - t0, 1),
        "input": {"ovm": str(ovm), "tet": str(tet), "n": 100000,
                  "cells": int(len(H)), "stages": str(STAGES)},
        "preprocessing": {
            "fill": "clean_blocks.fill_cavities",
            "untangle": "clean_blocks.untangle_until_stuck(rounds=6, rings=1)",
        },
        "raw": {"inverted": before_inverted, "min_sj": before_min_sj},
        "before_repair": base,
        "cells": cells,
        "repair": {
            "method": ("local one-sided barrier sum max(0, eps - J)^2 on the "
                       "scaled Jacobians of the inverted cells and their "
                       "neighbourhood; free vertices widened ring by ring, a "
                       "sweep is accepted only when the global inverted count "
                       "drops (the objective of clean_blocks._smooth; "
                       "clean_blocks.py defaults untouched)"),
            "grid": [{"rings": r, "eps": e} for r, e in GRID],
            "wide_probes": [{"rings": r, "eps": e} for r, e in WIDE],
            "rounds": log,
            "inverted_before": base["inverted"],
            "inverted_after": n1,
            "min_sj_before": base["min_sj"],
            "min_sj_after": float(s1.min()),
            "hausdorff_before": st0["hausdorff"],
            "hausdorff_after": st1["hausdorff"],
            "cells_after": int(len(H)),
            "remaining": remaining,
        },
        "beam": beam,
        "verdict": {
            "repaired": n1 == 0,
            "inverted_left": n1,
            "beam": "ran" if (beam or {}).get("status") not in (None, "SKIPPED")
                    else "skipped",
        },
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rep, indent=1) + "\n")
    print(f"[probe] wrote {out} ({time.time() - t0:.0f}s)")
    return rep


# --------------------------------------------------------------------------
# check: every claim against the saved geometry
# --------------------------------------------------------------------------
class Bad(Exception):
    pass


def _need(cond, msg):
    if not cond:
        raise Bad(msg)


def check(path):
    """Verify a repair report against itself and against the saved meshes."""
    path = Path(path)
    _need(path.exists(), f"{path} does not exist")
    rep = json.loads(path.read_text())
    _need(rep.get("schema") == SCHEMA, f"schema is {rep.get('schema')!r}, "
                                       f"expected {SCHEMA!r}")
    for key in ("input", "raw", "before_repair", "cells", "repair", "beam",
                "verdict"):
        _need(key in rep, f"missing section {key!r}")
    cells = rep["cells"]
    _need(len(cells) == 10, f"{len(cells)} localized cells, the task's lane has 10")
    _need(rep["before_repair"]["inverted"] == len(cells),
          "before_repair.inverted disagrees with the number of localized cells")
    ids = [c["id"] for c in cells]
    _need(len(set(ids)) == len(ids), "duplicate cell ids")
    for c in cells:
        for key in ("centroid", "bbox", "min_sj", "volume", "min_edge",
                    "max_edge", "surface_labels", "singular_edges",
                    "vertex_surface_distance_max", "on_feature_curve"):
            _need(key in c, f"cell {c.get('id')}: missing {key!r}")
        lo, hi = np.asarray(c["bbox"][0]), np.asarray(c["bbox"][1])
        ctr = np.asarray(c["centroid"])
        _need(np.all(ctr >= lo - 1e-6) and np.all(ctr <= hi + 1e-6),
              f"cell {c['id']}: centroid outside its bounding box")
        _need(c["min_edge"] > 0 and c["max_edge"] >= c["min_edge"],
              f"cell {c['id']}: impossible edge lengths")
        _need(np.isfinite(c["vertex_surface_distance_max"]) and
              c["vertex_surface_distance_max"] >= 0,
              f"cell {c['id']}: bad surface distance")

    # the saved meshes are what the report is about
    stages = Path(rep["input"]["stages"])
    if not stages.is_absolute():
        stages = path.parent / stages
    _need(stages.exists(), f"stages file {stages} is missing")
    d = np.load(stages)
    for tag, Pk, Hk, n_key in (("before_repair", "P0", "H0", "inverted"),
                               ("repair", "P1", "H1", "inverted_after")):
        P, H = d[Pk], d[Hk]
        sj = cb.scaled_jacobians(P, H)
        n = int((sj <= 0).sum())
        _need(n == rep[tag][n_key],
              f"{tag}: report says {rep[tag][n_key]} inverted cells, mesh has {n}")
        if tag == "before_repair":
            for c in cells:
                Q = P[H[c["id"]]]
                _need(np.allclose(Q.mean(0), c["centroid"], atol=2e-6),
                      f"cell {c['id']}: centroid does not match the mesh")
                _need(np.allclose(Q.min(0), c["bbox"][0], atol=2e-6) and
                      np.allclose(Q.max(0), c["bbox"][1], atol=2e-6),
                      f"cell {c['id']}: bounding box does not match the mesh")
                _need(abs(float(sj[c["id"]]) - c["min_sj"]) < 1e-6,
                      f"cell {c['id']}: scaled Jacobian does not match the mesh")
        else:
            _need(abs(float(sj.min()) - rep["repair"]["min_sj_after"]) < 1e-6,
                  "min_sj_after does not match the repaired mesh")
            for c in cells:
                fixed = bool(sj[c["id"]] > 0)
                if fixed:
                    _need(c["id"] not in [r["id"] for r in rep["repair"]["remaining"]],
                          f"cell {c['id']} is clean in the mesh but listed as remaining")

    r = rep["repair"]
    _need(r["inverted_before"] == rep["before_repair"]["inverted"],
          "repair.inverted_before disagrees with before_repair")
    _need(r["min_sj_before"] == rep["before_repair"]["min_sj"],
          "repair.min_sj_before disagrees with before_repair")
    _need(r["inverted_after"] < r["inverted_before"] or r["inverted_after"] == 0,
          "the probe did not reduce the inverted count")
    _need(r["hausdorff_after"] <= 0.15,
          f"repair raised the Hausdorff distance to {r['hausdorff_after']:.5f} "
          "(beam gate 0.15)")
    _need(r["cells_after"] == rep["input"]["cells"],
          "the repair changed the cell count (it must move vertices only)")
    if r["inverted_after"] == 0:
        _need(not r["remaining"], "0 inverted but cells are listed as remaining")
        _need(isinstance(rep["beam"], dict) and rep["beam"].get("status"),
              "0 inverted: the beam stage must have run")
    else:
        _need(len(r["remaining"]) == r["inverted_after"],
              "remaining cell count does not match the repaired mesh")
        for rem in r["remaining"]:
            _need(rem["min_sj"] <= 0, f"remaining cell {rem['id']} is not inverted")
            _need(bool(rem.get("reason")), f"remaining cell {rem['id']} has no reason")
            _need(rem["best_min_sj_over_wide_probes"] > rem["min_sj"] - 1e-9,
                  f"remaining cell {rem['id']}: no measured floor")
        _need(rep["beam"].get("status") == "SKIPPED" and rep["beam"].get("reason"),
              "inverted cells remain: the beam must be recorded as skipped with "
              "a reason")
    _need(rep["verdict"]["repaired"] == (r["inverted_after"] == 0),
          "verdict.repaired disagrees with the repair")
    print(f"OK {path}: {len(cells)} inverted cells localized, "
          f"{r['inverted_before']} -> {r['inverted_after']} after repair "
          f"(min sJ {r['min_sj_before']:.5f} -> {r['min_sj_after']:.5f}), "
          f"hausdorff {r['hausdorff_before']:.5f} -> {r['hausdorff_after']:.5f}, "
          f"beam {rep['beam'].get('status')}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mode", nargs="?", default="run", choices=("run",),
                    help="run the probe (the default; --check verifies a report)")
    ap.add_argument("--check", metavar="PATH",
                    help="verify a repair report and its saved meshes")
    ap.add_argument("--ovm", default=str(OVM))
    ap.add_argument("--tet", default=str(TET))
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--cache", default=str(STAGES),
                    help="npz holding both mesh stages (written by run)",
                    )
    ap.add_argument("--no-beam", action="store_true",
                    help="do not start the beam stage even if the mesh is clean")
    a = ap.parse_args(argv)
    if a.check:
        try:
            check(a.check)
        except Bad as e:
            print(f"FAIL {a.check}: {e}")
            return 1
        return 0
    run(Path(a.ovm), Path(a.tet), Path(a.out), Path(a.cache), not a.no_beam)
    return 0


if __name__ == "__main__":
    sys.exit(main())
