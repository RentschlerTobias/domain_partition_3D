#!/bin/bash
# One sample of the dataset per candidate geometry, strictly sequential.
#
#   generate_dataset.sh data/random_tistos/investigated/cand_0*  [-- extra AlgoHex flags]
#
# ~43 min per sample on an idle box: tet_prep 3, AlgoHex 10, clean_blocks 29,
# TFI 1. Two AlgoHex runs must never overlap (~6 GB peak each), and this box
# throttles after hours of sustained full load -- see HANDOFF.md. One at a
# time is therefore not a suggestion.
#
# A sample is REJECTED rather than repaired when a block has no lattice: the
# passthrough then pins whole direction classes through the complex (measured
# on cand_000: one 725-cell block forced the largest class, 20 axes, 21 %
# above what the free solve wants, and the refill came out with 3 inverted
# cells). For training data, a skipped sample is cheaper than a distorted one.
#
# Per candidate, writes into output/hex3d_algohex/gen/<name>/:
#   <name>_tet.vtk        AlgoHex input, exact surface labels
#   <name>_blocks.vtk     block structure  <- the training target
#   <name>_blocks.msh     same, with block edges as 1D elements
#   <name>_blocks_h05.vtk refilled at h=0.05, only if the structure passes
#   <name>_full.vtk       core + blade O-grid + both wall layers, CFD-ready
#   status.json           every stage's outcome and timing
set -u
cd /root/repos/duty/quadmesh/domain_partition_3D
PY=/root/repos/duty/quadmesh/.venv/bin/python
# unbuffered: otherwise a 29-minute stage shows an empty log the whole time
export PYTHONUNBUFFERED=1
E=experimentell/hex3d_algohex
H=${TARGET_H:-0.05}

for cand in "$@"; do
  name=$(basename "$cand")
  msh="$cand/mesh.msh"
  [ -f "$msh" ] || { echo "[$name] no mesh.msh, skipped"; continue; }
  d=output/hex3d_algohex/gen/$name
  mkdir -p "$d"
  t0=$(date +%s)
  echo "=== $(date -Is) $name ==="

  $PY $E/tet_prep_v5.py --msh "$msh" --out-dir "$d" --out ${name}_tet.vtk \
      > "$d/tet_prep.log" 2>&1 \
    || { echo "[$name] tet_prep FAILED"; echo '{"stage":"tet_prep","ok":false}' > "$d/status.json"; continue; }

  $PY $E/run_algohex.py --tag "$name" --prefix gen --in-vtk "$d/${name}_tet.vtk" \
      ${ALGOHEX_CPUS:+--cpus $ALGOHEX_CPUS} \
      -- -n 60000 > "$d/algohex.log" 2>&1 \
    || { echo "[$name] AlgoHex FAILED"; echo '{"stage":"algohex","ok":false}' > "$d/status.json"; continue; }

  $PY $E/clean_blocks.py output/hex3d_algohex/gen_hex_${name}.ovm \
      --input-vtk "$d/${name}_tet.vtk" \
      --collapse-rounds 5 --untangle --untangle-rounds 6 \
      --out "$d/${name}_blocks.vtk" > "$d/clean_blocks.log" 2>&1 \
    || { echo "[$name] clean_blocks FAILED"; echo '{"stage":"clean_blocks","ok":false}' > "$d/status.json"; continue; }

  # the CFD-ready variant: the AlgoHex core is the REDUCED domain, with the
  # blade O-grid and both wall layers cut out. reattach puts them back -- the
  # O-grid verbatim from the source MSH (5 blocks, dtOO's own split) and the
  # boundary layer extruded from the core's wall faces, which yields one block
  # per core block per wall because it has to stay conforming to the partition
  # above it. Both artifacts are kept: the core is what AlgoHex actually
  # produces and what a model would have to predict, the full domain is what a
  # solver can run.
  $PY $E/reattach.py "$d/${name}_blocks.vtk" --layers "${LAYERS:-17}" \
      --msh "$msh" --input-vtk "$d/${name}_tet.vtk" \
      --out "$d/${name}_full.vtk" > "$d/reattach.log" 2>&1 \
    || echo "[$name] reattach FAILED (core artifacts are still valid)"

  $PY $E/tfi.py "$d/${name}_blocks.vtk" --input-vtk "$d/${name}_tet.vtk" \
      --target-h "$H" --require-lattices --apply-divisions \
      --out "$d/${name}_blocks_h${H}.vtk" > "$d/tfi.log" 2>&1
  rc=$?
  dt=$(( $(date +%s) - t0 ))
  $PY - "$d" "$name" "$rc" "$dt" <<'PYEOF'
import json, re, sys
from pathlib import Path
d, name, rc, dt = Path(sys.argv[1]), sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
def grab(f, pat, cast=float):
    # re.M: the block counts are anchored per line, and without it the two
    # most important fields came out null
    m = re.search(pat, (d / f).read_text(), re.M) if (d / f).exists() else None
    return cast(m.group(1)) if m else None
out = {"name": name, "runtime_s": dt,
       "refill": {0: "ok", 3: "rejected: block without lattice"}.get(rc, f"failed rc={rc}"),
       "feature_edges": grab("tet_prep.log", r"feature edges: (\d+)", int),
       "blocks": grab("clean_blocks.log", r"^  (\d+) blocks,", int),
       "cuboids": grab("clean_blocks.log", r"^  \d+ blocks, (\d+) cuboids", int),
       "tiny": grab("clean_blocks.log", r"(\d+) blocks < 10 cells", int),
       "inverted": grab("clean_blocks.log", r"scaled Jacobian: (\d+) cells <= 0", int),
       "min_sj": grab("clean_blocks.log", r"min (-?[\d.]+), mean"),
       "hausdorff": grab("clean_blocks.log", r"boundary Hausdorff to input surface: ([\d.]+)"),
       "valid": "VALID" in (d / "clean_blocks.log").read_text(),
       "refill_cells": grab("tfi.log", r"-> (\d+) cells"),
       "refill_inverted": grab("tfi.log", r"mean [\d.]+, (\d+) inverted", int),
       "full_cells": grab("reattach.log", r"ASSEMBLED: (\d+) hexes", int),
       "full_blocks": grab("reattach.log", r"hexes, (\d+) blocks", int),
       "full_inverted": grab("reattach.log", r"mean [\d.]+, (\d+) inverted", int),
       "ogrid_gap_median": grab("reattach.log", r"SURFACE: median ([\d.]+)")}
(d / "status.json").write_text(json.dumps(out, indent=1))
print(f"[{name}] {out['blocks']} blocks, {out['tiny']} tiny, "
      f"{out['inverted']} inverted, min sJ {out['min_sj']}, "
      f"{'VALID' if out['valid'] else 'INVALID'}, refill {out['refill']}, {dt}s")
PYEOF
done
echo "=== $(date -Is) done ==="
