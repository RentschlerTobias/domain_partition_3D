#!/bin/bash
# Phase 2 of 2: compile AlgoHex and export the image. COMPUTE node, no network
# needed -- build_algohex_enroot.sh already downloaded everything.
#
#   sbatch --test-only experimentell/hex3d_algohex/scripts/build_algohex.slurm
#   sbatch            experimentell/hex3d_algohex/scripts/build_algohex.slurm
#
# Output: $WS/enroot-images/algohex.sqsh, which is exactly what
# run_algohex.py --backend enroot consumes.
#
# ---------------------------------------------------------------------------
# UNVERIFIED, CHECK BEFORE THE FIRST SUBMIT:
#   --partition   cpu_il: 64 cores, 256 GiB, documented and recommended in
#                 eigenfrequencies/docs/cluster-resource-sizing.md. NOT a dev
#                 partition -- dev_cpu_il caps at 30 min and this build has
#                 never run to completion anywhere.
#   --time        4 h is deliberate slack. On 64 cores the coinbrew chain plus
#                 ninja is plausibly 20-40 min, but MUMPS/Ipopt/Bonmin from
#                 source is the fragile part of this whole pipeline. A job that
#                 dies at the walltime wastes the whole build.
#   --cpus        ninja -j$(nproc) and coinbrew --parallel-jobs both scale.
# ---------------------------------------------------------------------------
#SBATCH --partition=cpu_il
#SBATCH --job-name=algohex_build
#SBATCH --output=algohex_build_%j.out
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=32G

set -u

CONTAINER="${CONTAINER:-algohex-build}"
WS="${WS:-$(ws_find hex3d 2>/dev/null || echo "$HOME")}"
OUT="${OUT:-$WS/enroot-images/algohex.sqsh}"

export ENROOT_DATA_PATH="${ENROOT_DATA_PATH:-$WS/enroot-data}"
export ENROOT_TEMP_PATH="${ENROOT_TEMP_PATH:-/tmp/$USER-enroot}"
export ENROOT_SQUASH_OPTIONS="${ENROOT_SQUASH_OPTIONS:--comp zstd -noD}"
mkdir -p "$ENROOT_TEMP_PATH" "$(dirname "$OUT")"

J="${SLURM_CPUS_PER_TASK:-$(nproc)}"
echo "[build-2] host $(hostname), job ${SLURM_JOB_ID:-none}, -j$J"

enroot list 2>/dev/null | grep -qxF "$CONTAINER" || {
    echo "[build-2] FAIL: no container '$CONTAINER'."
    echo "[build-2] Run build_algohex_enroot.sh on the login node first."
    exit 1; }

t0=$(date +%s)
echo "[build-2] --- 1/3 coin-or: Ipopt then Bonmin ---"
enroot start --root --rw --env NJ="$J" "$CONTAINER" bash -eux <<'INNER'
cd /usr/src/coin-or
./coinbrew build Ipopt@3.14.19  --verbosity 2 --skip-update \
    --prefix=/opt/coin-or --parallel-jobs "$NJ" --tests none
./coinbrew build Bonmin@master  --verbosity 2 --skip-update \
    --prefix=/opt/coin-or --parallel-jobs "$NJ" --tests none
# AlgoHex's cmake looks for coin/, the install puts coin-or/
ln -sfn /opt/coin-or/include/coin-or /opt/coin-or/include/coin
INNER
rc=$?
echo "[build-2] coin-or exit $rc after $(( $(date +%s) - t0 ))s"
[ $rc -eq 0 ] || { echo "[build-2] FAIL in the coinbrew chain -- the fragile part."; exit $rc; }

echo "[build-2] --- 2/3 AlgoHex: cmake Release, ninja -j$J ---"
enroot start --root --rw --env NJ="$J" "$CONTAINER" bash -eux <<'INNER'
export CoinUtils_DIR=/opt/coin-or IPOPT_HOME=/opt/coin-or \
       CBC_DIR=/opt/coin-or CLP_DIR=/opt/coin-or
mkdir -p /app/build && cd /app/build
# Release, not upstream's RelWithDebInfo: the one change from the local
# Dockerfile diff worth keeping (external_patches/algohex-Dockerfile.patch).
cmake -G Ninja -D CMAKE_BUILD_TYPE=Release -D BONMIN_ROOT_DIR=/opt/coin-or ..
# -j$NJ, not upstream's commented-out -j1. That -j1 is why the build was moved
# into a volume on a 2-vCPU box in the first place; here it is the point.
ninja -j"$NJ"

# A wrapper, not a bare symlink. enroot does not carry the Dockerfile's ENV,
# so a symlinked binary starts and cannot find its own .so files. This is what
# makes plain `HexMeshing` work from run_algohex.py with no env setup.
cat > /usr/local/bin/HexMeshing <<'WRAP'
#!/bin/sh
export LD_LIBRARY_PATH=/app/build/Build/lib:/opt/coin-or/lib:${LD_LIBRARY_PATH:-}
exec /app/build/Build/bin/HexMeshing "$@"
WRAP
chmod +x /usr/local/bin/HexMeshing
sed 's|bin/HexMeshing|bin/LocalMeshabilityCheck|' /usr/local/bin/HexMeshing \
    > /usr/local/bin/LocalMeshabilityCheck
chmod +x /usr/local/bin/LocalMeshabilityCheck

echo "--- self-test, the same gate T2 uses ---"
HexMeshing --help > /dev/null
ldd /app/build/Build/bin/HexMeshing | grep "not found" && exit 1
echo "self-test OK: HexMeshing runs, no missing libraries"
INNER
rc=$?
[ $rc -eq 0 ] || { echo "[build-2] FAIL in the AlgoHex build."; exit $rc; }

echo "[build-2] --- 3/3 export ---"
enroot export -f -o "$OUT" "$CONTAINER" || exit 1
ls -la "$OUT"
unsquashfs -s "$OUT" 2>/dev/null | grep -i "^Compression" || true

echo "[build-2] total $(( $(date +%s) - t0 ))s"
cat <<EOF
[build-2] done: $OUT

Provenance of this image: AlgoHex at the pinned commit, built from public
source by these two scripts. Nothing came from a Docker volume or another
machine.

Next: sbatch experimentell/hex3d_algohex/scripts/smoke_algohex.slurm
EOF
