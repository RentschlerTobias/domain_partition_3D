#!/bin/bash
# Phase 2 of 2: compile AlgoHex and export the image, on a compute node.
#
# AlgoHex's cmake DOWNLOADS its own dependencies -- ALGOHEX_DOWNLOAD_MISSING_DEPS
# is on by default when building standalone, per upstream's README -- and that
# works here, because compute nodes on this cluster have internet (measured
# 2026-09-11, job 6856932, COMPLETED). build_algohex_enroot.sh has already
# fetched the coin-or sources, so this job mostly compiles.
#
#   STAGE=ipopt   sbatch .../build_algohex.slurm
#   STAGE=bonmin  sbatch .../build_algohex.slurm
#   STAGE=algohex sbatch .../build_algohex.slurm
#   STAGE=all     sbatch --time=04:00:00 -p cpu_il .../build_algohex.slurm
#
# cluster_setup.sh submits the three stages chained with --dependency=afterok.
#
# ---------------------------------------------------------------------------
# WHY IN STAGES, AND WHY RESUMABLE
#
# dev_cpu_il caps at 30 minutes. cpu_il has room but measures about a WEEK of
# queue time, which makes a 4-hour single job the slowest possible route to a
# built image. So the build is cut at its natural seams -- Ipopt, Bonmin,
# AlgoHex -- and each piece is meant to fit a 30-minute dev slot.
#
# Every stage is idempotent in two ways, which is what makes the 30-minute cap
# survivable even if an estimate is wrong:
#
#   1. It checks for its own artifact first and exits immediately if present.
#      Re-running a finished stage costs seconds.
#   2. If it is KILLED at the walltime, resubmit the same stage. coinbrew drives
#      make and AlgoHex drives ninja; both keep their object files in the
#      container's rootfs, so a re-run picks up where the kill happened rather
#      than starting over.
#
# This is the same reason the original VPS build was made resumable: a
# non-resumable single-shot build against a hard limit loses everything at the
# limit.
#
# UNVERIFIED: whether each stage actually fits 30 min on 64 cores. Nobody has
# run this chain to completion anywhere. If a stage keeps hitting the wall,
# either resubmit it (it resumes) or run STAGE=all once on cpu_il and wait out
# the queue.
# ---------------------------------------------------------------------------
#SBATCH --job-name=algohex_build
#SBATCH --output=algohex_build_%j.out
#SBATCH --partition=dev_cpu_il
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=32G

set -u

STAGE="${STAGE:-all}"
CONTAINER="${CONTAINER:-algohex-build}"
WS="${WS:-$(ws_find hex3d 2>/dev/null || echo "$HOME")}"
OUT="${OUT:-$WS/enroot-images/algohex.sqsh}"

export ENROOT_DATA_PATH="${ENROOT_DATA_PATH:-$WS/enroot-data}"
export ENROOT_TEMP_PATH="${ENROOT_TEMP_PATH:-/tmp/${USER:-$(id -un)}-enroot}"
export ENROOT_SQUASH_OPTIONS="${ENROOT_SQUASH_OPTIONS:--comp zstd -noD}"
# enroot's runtime state. Default is ${XDG_RUNTIME_DIR}/enroot, and that is a
# trap in a batch job: --export=ALL carries XDG_RUNTIME_DIR=/run/user/$UID from
# the login node, where logind created it, into a compute node where it does
# not exist and /run/user is not writable. Measured 2026-09-11: two enroot
# calls succeeded and the third died with "mkdir: cannot create directory
# '/run/user/985462': Permission denied" -- the session that owned the
# directory had ended in between. Pin it somewhere we own.
export ENROOT_RUNTIME_PATH="${ENROOT_RUNTIME_PATH:-/tmp/${USER:-$(id -un)}-enroot-run}"
mkdir -p "$ENROOT_TEMP_PATH" "$ENROOT_RUNTIME_PATH" "$(dirname "$OUT")"

J="${SLURM_CPUS_PER_TASK:-$(nproc)}"
T0=$(date +%s)
echo "[build] stage=$STAGE host=$(hostname) job=${SLURM_JOB_ID:-none} -j$J"
echo "[build] walltime for this slot: ${SLURM_JOB_END_TIME:+ends $(date -d @"$SLURM_JOB_END_TIME" 2>/dev/null)}"

enroot list 2>/dev/null | grep -qxF "$CONTAINER" || {
    echo "[build] FAIL: no container '$CONTAINER'."
    echo "[build] Run build_algohex_enroot.sh on the login node first."
    exit 1; }

inside() { enroot start --root --rw --env NJ="$J" "$CONTAINER" bash -eux; }
elapsed() { echo "[build] ${1} after $(( $(date +%s) - T0 ))s"; }

# --------------------------------------------------------------------- ipopt
if [ "$STAGE" = ipopt ] || [ "$STAGE" = all ]; then
    if enroot start --root "$CONTAINER" test -f /opt/coin-or/lib/libipopt.so \
         2>/dev/null; then
        echo "[build] ipopt: already installed, skipping"
    else
        echo "[build] --- ipopt (with MUMPS) ---"
        inside <<'INNER'
cd /usr/src/coin-or
./coinbrew build Ipopt@3.14.19 --verbosity 2 --skip-update \
    --prefix=/opt/coin-or --parallel-jobs "$NJ" --tests none
test -f /opt/coin-or/lib/libipopt.so
INNER
        rc=$?; elapsed "ipopt exit $rc"
        [ $rc -eq 0 ] || { echo "[build] FAIL in Ipopt. Resubmit STAGE=ipopt; make resumes."; exit $rc; }
    fi
fi

# -------------------------------------------------------------------- bonmin
if [ "$STAGE" = bonmin ] || [ "$STAGE" = all ]; then
    if enroot start --root "$CONTAINER" test -f /opt/coin-or/lib/libbonmin.so \
         2>/dev/null; then
        echo "[build] bonmin: already installed, skipping"
    else
        echo "[build] --- bonmin (the moving part: upstream takes master) ---"
        inside <<'INNER'
cd /usr/src/coin-or
./coinbrew build Bonmin@master --verbosity 2 --skip-update \
    --prefix=/opt/coin-or --parallel-jobs "$NJ" --tests none
# AlgoHex's cmake looks for coin/, the install puts coin-or/
ln -sfn /opt/coin-or/include/coin-or /opt/coin-or/include/coin
test -f /opt/coin-or/lib/libbonmin.so
INNER
        rc=$?; elapsed "bonmin exit $rc"
        [ $rc -eq 0 ] || { echo "[build] FAIL in Bonmin -- the fragile part. Resubmit STAGE=bonmin."; exit $rc; }
    fi
fi

# ------------------------------------------------------- algohex, then export
if [ "$STAGE" = algohex ] || [ "$STAGE" = all ]; then
    echo "[build] --- AlgoHex: cmake Release, ninja -j$J ---"
    inside <<'INNER'
export CoinUtils_DIR=/opt/coin-or IPOPT_HOME=/opt/coin-or \
       CBC_DIR=/opt/coin-or CLP_DIR=/opt/coin-or
mkdir -p /app/build && cd /app/build
# Release, not upstream's RelWithDebInfo: the one change from the local
# Dockerfile diff worth keeping (external_patches/algohex-Dockerfile.patch).
# cmake is cheap and idempotent, so it re-runs on every resubmission.
cmake -G Ninja -D CMAKE_BUILD_TYPE=Release -D BONMIN_ROOT_DIR=/opt/coin-or ..
# -j$NJ, not upstream's commented-out -j1. That -j1 is why the build was moved
# into a volume on a 2-vCPU box in the first place; here it is the point.
# Killed at the walltime? ninja keeps its objects -- resubmit and it continues.
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
    rc=$?; elapsed "algohex exit $rc"
    [ $rc -eq 0 ] || { echo "[build] FAIL in the AlgoHex build. Resubmit STAGE=algohex; ninja resumes."; exit $rc; }

    echo "[build] --- export ---"
    enroot export -f -o "$OUT" "$CONTAINER" || exit 1
    ls -la "$OUT"
    unsquashfs -s "$OUT" 2>/dev/null | grep -i "^Compression" || true
    elapsed "export done"

    cat <<EOF
[build] done: $OUT

Provenance: AlgoHex at the pinned commit, built from public source by these two
scripts. Nothing came from a Docker volume or another machine.

Next: sbatch experimentell/hex3d_algohex/scripts/smoke_algohex.slurm
EOF
fi

elapsed "stage=$STAGE complete"
