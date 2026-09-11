#!/usr/bin/env python
"""Generate machine candidates for the dataset: sampler -> params.json -> dtOO mesh.

T14 of the dataset pipeline, implemented per
docs/decisions/2026-09-11-dataset-sampler-strategy.md (decisions Q1-Q4):

* Q1: Sobol over the 30-dimensional design box, seeded, incremental — the
  first k points of a scrambled Sobol sequence stay valid, so a killed run
  resumes by drawing k+1..N. ``--strategy`` is the plugin seam for the later
  optimizer: its DE population becomes another strategy, not a rewrite.
* Q2: state lives in ``sampler.json`` (seed, strategy, bounds hash, counts);
  a machine directory's presence IS the idempotency marker. Failed exports
  are retried at the back of the queue, same pattern as sample_one.sh.
* Q3: layout ``data/dataset/<strategy>/machine_00NN/`` — ``machine_00NN`` is
  the draw number and the load-bearing address. The 21 reference candidates
  stay untouched under data/random_tistos/investigated/.
* Q4: the driver lives here; the dtOO world is reached through ONE import
  seam, the sibling ``eigenfrequencies`` repo. The whole pipeline is
  destined to be integrated into that repo's optimization loop, so its
  container/import is used as-is, never copied.

The output dirs are CANDS-compatible with scripts/batch_samples.slurm:
``mesh.msh`` + ``params.json`` per machine, exactly what that manifest
branch consumes.

Usage:
    # draw + write params.json only, no dtOO, no cluster needed:
    python generate_machines.py --count 64 --preview

    # on the cluster, inside the dtOO container: also export mesh.msh
    python generate_machines.py --count 64 --export

Resume is automatic: rerun with the same --seed/--strategy/--count and the
existing machine dirs are skipped.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import math
import os
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]  # quadmesh/domain_partition_3D
E = REPO / "experimentell" / "hex3d_algohex"

# Q4: the single import seam. eigenfrequencies is a sibling repo; on the
# cluster it is mounted into the dtOO container per
# eigenfrequencies/cluster/enroot_dtoo_import.md.
EF_ROOT = Path(os.environ.get("EIGENFREQUENCIES_ROOT", REPO.parent.parent / "eigenfrequencies"))
EF_SRC = EF_ROOT / "src"
DEFAULT_YAML = EF_ROOT / "adapters" / "machines" / "tistos.yaml"


def load_adapter(yaml_path: Path):
    """Import DtooAdapter through Q4's path insert. Importable without dtOO:
    the heavy SWIG import only happens inside export_mesh."""
    sys.path.insert(0, str(EF_SRC))
    from eigenfrequencies.adapters.dtoo.adapter import DtooAdapter  # noqa: E402

    return DtooAdapter(yaml_path)


def bounds_digest(bounds: dict[str, tuple[float, float]]) -> str:
    """A hash over the ordered label:min:max list — sampler.json uses it to
    fail loudly when the machine YAML's bounds changed under an existing
    dataset root instead of silently mixing two spaces."""
    payload = "\n".join(f"{k}:{v[0]!r}:{v[1]!r}" for k, v in sorted(bounds.items()))
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def draw_points(
    strategy: str, count: int, seed: int, bounds: dict[str, tuple[float, float]]
) -> list[dict[str, float]]:
    """N design vectors as {label: value}. count MUST be a power of two for
    Sobol: the sequence's balance guarantees only hold at those sizes. The
    single seam a later optimizer strategy plugs into."""
    if strategy != "sobol":  # future: "de" / "lhs" plugins
        raise SystemExit(f"unknown strategy: {strategy}")
    if count & (count - 1) or count < 2:
        raise SystemExit("--count must be a power of two >= 2 (Sobol balance)")
    from scipy.stats import qmc

    labels = sorted(bounds)
    d = len(labels)
    sop = qmc.Sobol(d=d, scramble=True, seed=seed)
    unit = sop.random(count)  # (count, d) in [0, 1)
    lo = [bounds[k][0] for k in labels]
    hi = [bounds[k][1] for k in labels]
    scaled = qmc.scale(unit, lo, hi)
    return [
        dict(zip(labels, [round(float(v), 6) for v in row])) for row in scaled
    ]


def machine_dirs(root: Path) -> list[Path]:
    return sorted(d for d in root.glob("machine_00*") if d.is_dir())


def pending(dirs: list[Path], want_mesh: bool) -> list[Path]:
    """Q2: presence is the marker. --export only looks for mesh.msh; the
    params draw itself is already durable."""
    if want_mesh:
        return [d for d in dirs if not (d / "mesh.msh").exists()]
    return []


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--count", type=int, default=64, help="power-of-two draw count")
    p.add_argument("--seed", type=int, default=20260911, help="sobol scramble seed")
    p.add_argument("--strategy", default="sobol", help="sampler strategy (seam, Q1)")
    p.add_argument("--yaml", type=Path, default=DEFAULT_YAML, help="machine YAML")
    p.add_argument(
        "--root",
        type=Path,
        default=None,
        help="dataset root (default: data/dataset/<strategy>)",
    )
    p.add_argument("--preview", action="store_true", help="draw + params.json only")
    p.add_argument("--export", action="store_true", help="also run dtOO (container only)")
    p.add_argument(
        "--only",
        default=None,
        help="export only this machine (e.g. machine_0007); the batch loop "
             "spawns one fresh python process per machine so dtOO SWIG state "
             "cannot accumulate (segfault risk, see dtoo_cfd_build.py)",
    )
    args = p.parse_args()

    if not EF_SRC.is_dir():
        raise SystemExit(f"eigenfrequencies not found at {EF_ROOT} (set EIGENFREQUENCIES_ROOT)")

    adapter = load_adapter(args.yaml)
    bounds = adapter.design_bounds()
    if not bounds:
        raise SystemExit(f"no design bounds in {args.yaml}")

    root = args.root or (E.parent.parent / "data" / "dataset" / args.strategy)
    root.mkdir(parents=True, exist_ok=True)

    digest = bounds_digest(bounds)
    state_path = root / "sampler.json"
    if state_path.exists():
        state = json.loads(state_path.read_text())
        for key, want in (("seed", args.seed), ("strategy", args.strategy), ("bounds_sha256_16", digest)):
            if state.get(key) != want:
                raise SystemExit(
                    f"{state_path} was generated with {key}={state.get(key)!r}, "
                    f"now asked for {want!r}. Point --root elsewhere or match."
                )
    else:
        state = {
            "strategy": args.strategy,
            "seed": args.seed,
            "bounds_sha256_16": digest,
            "yaml": str(args.yaml),
            "labels": sorted(bounds),
            "created_by": "generate_machines.py",
        }

    # On the cluster the --export pass runs INSIDE the dtOO container, whose
    # python3.13 provably has dtOOPythonSWIG (the smoke test) but NOT
    # necessarily scipy.stats (only the draw needs that). When every machine
    # directory already carries its params.json the draw is redundant -- the
    # values are durable in the filesystem itself (Q2) -- so skip it and keep
    # the container dependency down to the adapter.
    drawn = machine_dirs(root)
    if not (args.export
            and (root / "sampler.json").exists()
            and drawn
            and all((d / "params.json").exists() for d in drawn)
            and len(drawn) >= args.count):
        points = draw_points(args.strategy, args.count, args.seed, bounds)
        new = 0
        for i, values in enumerate(points):
            d = root / f"machine_{i + 1:04d}"  # draw NN, 1-based (decision Q3)
            pj = d / "params.json"
            if pj.exists():  # Q2: presence is the idempotency marker
                continue
            d.mkdir(parents=True, exist_ok=True)
            pj.write_text(json.dumps(values, indent=2, sort_keys=True) + "\n")
            new += 1

        state["count_target"] = args.count
        state["machines"] = len(machine_dirs(root))
        state["updated"] = datetime.datetime.now(datetime.timezone.utc).isoformat(
            timespec="seconds")
        state_path.write_text(json.dumps(state, indent=2, sort_keys=True) + "\n")

        print(f"[sampler] {args.strategy} seed={args.seed} target={args.count} "
              f"bounds={digest}: {new} new, {args.count - new} existing")

        if not args.preview and not args.export:
            print("[sampler] drawn only -- pass --preview as reminder, or --export on the cluster")
    else:
        print(f"[sampler] draw skipped: all {len(drawn)} machines carry params.json "
              f"(export pass, scipy not needed here)")

    if args.export:
        todo = pending(machine_dirs(root), True)
        if args.only:
            todo = [d for d in todo if d.name == args.only]
            if not todo:
                raise SystemExit(
                    f"--only {args.only}: not a pending machine under {root} "
                    "(unknown name or mesh.msh already present)")
        # A SIGTERM from SLURM mid-export leaves a torn mesh.msh.suffix next
        # to the finished ones. Presence of mesh.msh is Q2's idempotency
        # marker, so the suffix is not just cosmetic: without it the next
        # resubmit would feed a torn mesh into tet_prep. Write to the suffix,
        # rename only on success, sweep stale suffixes before anything runs.
        for d in machine_dirs(root):
            (d / "mesh.msh.part").unlink(missing_ok=True)
        print(f"[sampler] exporting {len(todo)} meshes (dtOO, container only)")
        failed: list[str] = []
        for d in todo:
            values = json.loads((d / "params.json").read_text())
            try:
                # export_mesh(design_values) -> str has no output kwarg
                # (adapter.py:35); location comes from DTOO_OUTPUT_MSH, which
                # export.py reads. Older checkouts ignore the env and return
                # their default path, which is moved (cross-device safe).
                part = d / "mesh.msh.part"
                os.environ["DTOO_OUTPUT_MSH"] = str(part)
                mesh = adapter.export_mesh(values)
                if part.is_file():
                    os.replace(part, d / "mesh.msh")
                elif mesh:
                    shutil.move(str(mesh), str(d / "mesh.msh"))
                else:
                    raise RuntimeError("export_mesh returned no mesh path")
                (d / "export_error.txt").unlink(missing_ok=True)
                print(f"[sampler] {d.name} mesh OK")
            except Exception as e:  # noqa: BLE001 - yield measurement needs every failure
                failed.append(d.name)
                (d / "export_error.txt").write_text(f"{type(e).__name__}: {e}\n")
                print(f"[sampler] {d.name} export FAILED: {type(e).__name__}")
        if failed:
            print(f"[sampler] yield note: {len(failed)} failures recorded in export_error.txt; "
                  "a resubmit retries them, .part tearing excluded by rename-on-success")

    if args.preview:
        # Coverage proof: realized per-dimension span next to the bounds. The
        # point of the preview is to eyeball exactly this.
        labels = sorted(bounds)
        vals = [json.loads((d / "params.json").read_text()) for d in machine_dirs(root)]
        print(f"[sampler] coverage over {len(vals)} machines:")
        for k in labels:
            got = [v[k] for v in vals]
            lo, hi = bounds[k]
            print(f"  {k:34s} bounds [{lo:+.4g},{hi:+.4g}] "
                  f"realized [{min(got):+.4g},{max(got):+.4g}] span {(max(got)-min(got))/(hi-lo):.2f}")


if __name__ == "__main__":
    main()
