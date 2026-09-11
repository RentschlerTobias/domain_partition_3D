#!/usr/bin/env bash
# One sample, end to end. The unit of work that batch_samples.slurm runs many
# of in parallel.
#
#   sample_one.sh <tet.vtk> <n> <outdir> [params.json]
#
# AlgoHex -> clean_blocks -> tfi -> export_sample, all logs inside <outdir>.
# Writes <outdir>/sample.npz and nothing outside <outdir>.
#
# IDEMPOTENT: if <outdir>/sample.npz already exists, it exits 0 immediately.
# That is what makes a killed batch job cheap to resubmit -- the survivors are
# skipped in milliseconds.
#
# Exit codes are the batch summary's data: 0 ok, 2 algohex, 3 clean_blocks,
# 4 tfi, 5 export, 6 bad arguments.

set -u

TET="${1:?usage: sample_one.sh <tet.vtk> <n> <outdir> [params.json]}"
N="${2:?missing -n}"
D="${3:?missing outdir}"
PARAMS="${4:-}"

REPO="${REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
E="$REPO/experimentell/hex3d_algohex"
PY="${PY:-$REPO/.venv/bin/python}"
H="${TARGET_H:-0.05}"
IMG="${ALGOHEX_IMAGE:-algohex}"
BACKEND="${BACKEND:-enroot}"
# Keep the AlgoHex .ovm? Decision J says no at dataset scale: 13 MB per sample
# at n=60000 is 130 GB across 10 000 samples. Worth keeping for a showcase set.
KEEP_OVM="${KEEP_OVM:-0}"

[ -f "$TET" ] || { echo "[$D] no such tet mesh: $TET"; exit 6; }
mkdir -p "$D" || exit 6

if [ -f "$D/sample.npz" ]; then
    echo "[$(basename "$D")] already done, skipping"
    exit 0
fi

t0=$(date +%s)
name=$(basename "$D")
log() { echo "[$name] $*"; }

log "start  tet=$(basename "$TET") n=$N h=$H"

"$PY" "$E/run_algohex.py" --backend "$BACKEND" --image "$IMG" \
    --prefix hex --tag "$name" --in-vtk "$TET" --out-dir "$D" \
    -- -n "$N" > "$D/algohex.log" 2>&1 \
  || { log "FAIL algohex, see $D/algohex.log"; exit 2; }

"$PY" "$E/clean_blocks.py" "$D/hex_hex_${name}.ovm" --input-vtk "$TET" \
    --collapse-rounds 5 --untangle --untangle-rounds 6 \
    --out "$D/blocks.vtk" > "$D/clean_blocks.log" 2>&1 \
  || { log "FAIL clean_blocks, see $D/clean_blocks.log"; exit 3; }

# --require-lattices: a block without a lattice pins whole direction classes
# through the complex, and the refill then comes out distorted. For training
# data a skipped sample is cheaper than a distorted one (generate_dataset.sh
# has the measurement on cand_000).
"$PY" "$E/tfi.py" "$D/blocks.vtk" --input-vtk "$TET" --target-h "$H" \
    --require-lattices --solve-divisions > "$D/tfi.log" 2>&1
rc=$?
[ $rc -eq 0 ] || { log "FAIL tfi rc=$rc, see $D/tfi.log"; exit 4; }

"$PY" "$E/export_sample.py" "$D/blocks.vtk" --tet "$TET" --target-h "$H" \
    --n "$N" ${PARAMS:+--params "$PARAMS"} \
    --out "$D/sample.npz" --check > "$D/export_sample.log" 2>&1 \
  || { log "FAIL export, see $D/export_sample.log"; exit 5; }

# Retention, decision J: the sample plus blocks.vtk are the recovery boundary,
# everything else is derivable from them in seconds.
[ "$KEEP_OVM" = 1 ] || rm -f "$D"/*.ovm
rm -f "$D"/blocks_nfaces.vtk "$D"/blocks_quality.vtk "$D"/blocks.msh

sz=$(stat -c %s "$D/sample.npz" 2>/dev/null || echo 0)
log "OK     $(( $(date +%s) - t0 ))s, sample $(( sz / 1024 )) KB, kept $(du -sk "$D" | cut -f1) KB total"
