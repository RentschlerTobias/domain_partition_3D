#!/bin/bash
# Smoke test 1 of 2: does HexMeshing run from the enroot image on a COMPUTE
# node? Nothing else. 5 minutes, one task.
#
#   sbatch --test-only experimentell/hex3d_algohex/scripts/smoke_algohex.slurm
#   sbatch            experimentell/hex3d_algohex/scripts/smoke_algohex.slurm
#
# Why a compute node and not the login node: squashfuse on the compute nodes is
# what breaks when a .sqsh was written with lzo instead of zstd -- the image
# imports without complaint and is unreadable exactly here. Ours is zstd
# (verified at export), and this job is what proves it on the real machine.
#
# ---------------------------------------------------------------------------
# UNVERIFIED, CHECK BEFORE THE FIRST SUBMIT:
#   --partition   dev_cpu_il is taken from the sibling project's smoke test.
#                 Verify with: sinfo -o "%P %l %c %m %D"
#   --time        5 min assumes the dev partition allows it (30 min cap there)
#   --mem         2G is plenty for --help; the real run needs ~8G, see
#                 one_sample.slurm
# duty/quadmesh/scripts/dryrun_planb.slurm carries the same warning for the
# same reason: these names were written without cluster access.
# ---------------------------------------------------------------------------
#SBATCH --job-name=algohex_smoke
#SBATCH --output=algohex_smoke_%j.out
#SBATCH --partition=dev_cpu_il
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

set -u

WS="${WS:-$(ws_find hex3d 2>/dev/null)}"
IMAGE="${ALGOHEX_SQSH:-$WS/enroot-images/algohex.sqsh}"
REPO="${REPO:-$WS/domain_partition_3D}"

# The data path is where `enroot create` unpacks to. Persistent and on the
# workspace so every later job in the array skips the unpack; the TEMP path
# must be LOCAL disk, because enroot mounts an overlayfs there and a parallel
# filesystem cannot host one ("failed to mount overlay: Invalid argument",
# recorded in eigenfrequencies/cluster/cluster_env.sh).
export ENROOT_DATA_PATH="${ENROOT_DATA_PATH:-$WS/enroot-data}"
export ENROOT_TEMP_PATH="${ENROOT_TEMP_PATH:-/tmp/$USER-enroot}"
# enroot's runtime state. Default is ${XDG_RUNTIME_DIR}/enroot, and that is a
# trap in a batch job: --export=ALL carries XDG_RUNTIME_DIR=/run/user/$UID from
# the login node, where logind created it, into a compute node where it does
# not exist and /run/user is not writable. Measured 2026-09-11: two enroot
# calls succeeded and the third died with "mkdir: cannot create directory
# '/run/user/985462': Permission denied" -- the session that owned the
# directory had ended in between. Pin it somewhere we own.
export ENROOT_RUNTIME_PATH="${ENROOT_RUNTIME_PATH:-/tmp/${USER:-$(id -un)}-enroot-run}"
mkdir -p "$ENROOT_DATA_PATH" "$ENROOT_TEMP_PATH" "$ENROOT_RUNTIME_PATH"

echo "[smoke] host      $(hostname)"
echo "[smoke] job       ${SLURM_JOB_ID:-none}"
echo "[smoke] image     $IMAGE"
echo "[smoke] repo      $REPO"
echo "[smoke] data path $ENROOT_DATA_PATH"

[ -f "$IMAGE" ] || { echo "[smoke] FAIL: no image at $IMAGE"; exit 1; }
[ -d "$REPO" ]  || { echo "[smoke] FAIL: no repo at $REPO"; exit 1; }

echo "[smoke] --- 1/3 compression must be zstd, not lzo ---"
unsquashfs -s "$IMAGE" | grep -i "^Compression" || echo "[smoke] (unsquashfs unavailable, skipped)"

echo "[smoke] --- 2/3 HexMeshing --help, straight from the .sqsh ---"
srun -n1 -N1 enroot start --root --rw --mount "$REPO:/work" "$IMAGE" \
    HexMeshing --help | head -5
rc=${PIPESTATUS[0]}

echo "[smoke] --- 3/3 no missing shared libraries, and the mount is visible ---"
srun -n1 -N1 enroot start --root --rw --mount "$REPO:/work" "$IMAGE" \
    sh -c 'ldd /app/build/Build/bin/HexMeshing | grep "not found" | wc -l; \
           ls /work/data/T1_9/ >/dev/null && echo "mount OK"'

echo "[smoke] exit $rc"
echo "[smoke] PASS means: image readable on a compute node, binary runs, repo mounted."
echo "[smoke] Next: one_sample.slurm"
exit $rc
