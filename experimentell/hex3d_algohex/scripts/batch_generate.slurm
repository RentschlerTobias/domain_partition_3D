#!/bin/bash
# One production job for the whole dataset: dtOO meshes -> hex samples.
#
#   sbatch experimentell/hex3d_algohex/scripts/batch_generate.slurm
#   sbatch -t 04:00:00 .../batch_generate.slurm          # a shorter request backfills sooner
#   NS="2000 8000 60000" PARALLEL=32 \
#       sbatch .../batch_generate.slurm                  # larger/future runs
#
# WHY cpu_il. The dev queues cap at 30 minutes, which fits exactly one wave of
# samples; every batch so far was a resubmit chain. This job is the escape
# hatch the cluster notes describe: cpu_il allows 3 days, so ONE queue wait
# covers the whole generation instead of a dozen 30-minute slots. The price is
# queue time (about a week, measured) -- so this is the job for a full or big
# generation, not for finishing a handful of missing samples. Those are still
# fastest as dev_cpu_il resubmits.
#
# WHY ONE NODE, for now. The measured optimum on 64 cores / 128 GB is 32
# workers with OMP_NUM_THREADS=2 (see the batch_samples.slurm header): roughly
# 80-100 samples per hour, so one 8-hour slot covers ~600-800 samples. More
# nodes would parallelize the samples linearly, but the algohex container runs
# with --root --rw from a shared unpacked rootfs, so concurrent nodes need one
# container copy each (plus a per-node share of the manifest). Deliberately not
# built yet; revisit when the target dataset outgrows a single node.
#
# The two phase scripts are chained unchanged on purpose. Both are idempotent
# (a finished mesh.msh / sample.npz is skipped), each keeps its per-job log
# under logs/, and their summary blocks -- the dtOO yield and the sample
# throughput -- are the measurements DATASET_PIPELINE.md hangs off.
#SBATCH --job-name=hex3d_generate
#SBATCH --output=hex3d_generate_%j.out
#SBATCH --partition=cpu_il
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G

set -u

WS="$(ws_find hex3d 2>/dev/null)"
REPO="${REPO:-${SLURM_SUBMIT_DIR:-$WS/domain_partition_3D}}"
[ -n "$WS" ] && [ -d "$REPO/experimentell/hex3d_algohex" ] || {
    echo "[gen] not the hex3d checkout: WS=$WS REPO=$REPO"; exit 1; }

LOGDIR="$REPO/experimentell/hex3d_algohex/logs"
mkdir -p "$LOGDIR"
LOG="$LOGDIR/batch_generate_${SLURM_JOB_ID:-manual}.log"
exec > >(tee "$LOG") 2>&1

cd "$REPO" || exit 1
E=experimentell/hex3d_algohex

# Thread budget and concurrency; both phase scripts read these from the
# environment. 32 x 2 is the measured optimum (48 x 1 produced fewer samples
# per slot), and sample_one.sh already forwards the pin into the container.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
export PARALLEL="${PARALLEL:-32}"   # sample phase: xargs workers
export PAR="${PAR:-32}"             # dtOO phase: parallel exports
export NS="${NS:-2000 8000}"
export CANDS="${CANDS:-data/dataset/sobol/machine_00*}"

echo "[gen] host $(hostname), job ${SLURM_JOB_ID:-none}"
echo "[gen] NS='$NS' PAR=$PAR PARALLEL=$PARALLEL OMP_NUM_THREADS=$OMP_NUM_THREADS"

echo "[gen] ================ phase 1/2: dtOO meshes ================"
bash "$E/scripts/batch_dtoo_export.slurm"
RC_DTOO=$?
echo "[gen] phase 1 done, rc=$RC_DTOO"

echo "[gen] ================ phase 2/2: hex samples ================"
bash "$E/scripts/batch_samples.slurm"
RC_SAMPLES=$?
echo "[gen] phase 2 done, rc=$RC_SAMPLES"

echo "[gen] ================ run complete ================"
echo "[gen] dtOO rc=$RC_DTOO, samples rc=$RC_SAMPLES (0 = nothing missing)"
echo "[gen] log: $LOG"
# let the tee drain before SLURM sees the job as done
wait
[ "$RC_SAMPLES" -eq 0 ]
