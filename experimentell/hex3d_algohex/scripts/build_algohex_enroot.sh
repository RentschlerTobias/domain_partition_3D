#!/usr/bin/env bash
# Phase 1 of 2: prepare the build container. Run it on the login node.
#
# NOT because compute nodes lack internet -- they have it. Measured 2026-09-11:
# AlgoHex's cmake downloaded OpenVolumeMesh, libHexEx, CoMISo, GMM, TinyAD,
# QGP3D, Eigen and CLI11 from inside a dev_cpu_il job and the build completed.
# An earlier version of this file claimed the opposite, taken from a note in
# the sibling eigenfrequencies repo rather than from this cluster, and it would
# have split the build across an unnecessary manual step.
#
# The reason to download here is politeness and speed: the login node is not
# billed against a 30-minute slot, and doing it once means the batch job spends
# its whole slot compiling.
#
#   bash experimentell/hex3d_algohex/scripts/build_algohex_enroot.sh
#
# Builds AlgoHex from source with NO Docker anywhere: enroot runs upstream's
# Dockerfile steps as what they are, a list of shell commands. Verified on the
# VPS 2026-09-11 in miniature -- apt-get install inside `enroot start --root
# --rw` works (enroot even ships 10-aptfix.sh for Debian rootfs), the changes
# persist across invocations, and `enroot export` captured them (87 -> 174 MB).
#
# Split where upstream's Dockerfile splits: everything that DOWNLOADS runs
# here, everything that COMPILES runs in build_algohex.slurm on a node with
# cores. A convenience, not a constraint -- see above.
#
# Provenance after this route: a pinned public commit plus these two scripts.
# No Docker volume, no image copied from anyone's disk, no registry.
#
# ---------------------------------------------------------------------------
# Cost: 10-20 min, mostly download. Light on CPU -- fine for a login node.
# Space: ~6 GB in $ENROOT_DATA_PATH (coinbrew fetches 481+232+84 MB of source,
#        plus apt cache and the AlgoHex clone). Put it on the WORKSPACE.
# ---------------------------------------------------------------------------

set -euo pipefail

ALGOHEX_COMMIT="${ALGOHEX_COMMIT:-3519289}"   # pinned; the clone this repo used
BASE_IMAGE="${BASE_IMAGE:-docker://debian:trixie-20251117}"  # upstream's FROM
CONTAINER="${CONTAINER:-algohex-build}"

WS="${WS:-$(ws_find hex3d 2>/dev/null || echo "$HOME")}"
export ENROOT_DATA_PATH="${ENROOT_DATA_PATH:-$WS/enroot-data}"
export ENROOT_CACHE_PATH="${ENROOT_CACHE_PATH:-$WS/enroot-cache}"
# LOCAL disk, never the workspace: enroot flattens image layers through an
# overlayfs and a parallel filesystem cannot host one -- the import dies with
# "failed to mount overlay: Invalid argument" AFTER the download finished.
# Recorded in eigenfrequencies/cluster/cluster_env.sh.
export ENROOT_TEMP_PATH="${ENROOT_TEMP_PATH:-/tmp/$USER-enroot}"
# zstd, not lzo: an lzo image imports without complaint and is then unreadable
# on every compute node.
export ENROOT_SQUASH_OPTIONS="${ENROOT_SQUASH_OPTIONS:--comp zstd -noD}"
# enroot's runtime state. Default is ${XDG_RUNTIME_DIR}/enroot, and that is a
# trap in a batch job: --export=ALL carries XDG_RUNTIME_DIR=/run/user/$UID from
# the login node, where logind created it, into a compute node where it does
# not exist and /run/user is not writable. Measured 2026-09-11: two enroot
# calls succeeded and the third died with "mkdir: cannot create directory
# '/run/user/985462': Permission denied" -- the session that owned the
# directory had ended in between. Pin it somewhere we own.
export ENROOT_RUNTIME_PATH="${ENROOT_RUNTIME_PATH:-/tmp/${USER:-$(id -un)}-enroot-run}"
mkdir -p "$ENROOT_DATA_PATH" "$ENROOT_CACHE_PATH" "$ENROOT_TEMP_PATH" \
         "$ENROOT_RUNTIME_PATH" \
         "$WS/enroot-images"

echo "[build-1] workspace   $WS"
echo "[build-1] container   $CONTAINER"
echo "[build-1] AlgoHex     $ALGOHEX_COMMIT"
df -h "$WS" | tail -1

if enroot list 2>/dev/null | grep -qxF "$CONTAINER"; then
    echo "[build-1] container exists already."
    echo "[build-1] To start over: enroot remove -f $CONTAINER"
    exit 0
fi

echo "[build-1] --- 1/3 base image ---"
enroot import -o "$WS/enroot-images/debian-base.sqsh" "$BASE_IMAGE"
enroot create --name "$CONTAINER" "$WS/enroot-images/debian-base.sqsh"

echo "[build-1] --- 2/3 packages, exactly upstream's list ---"
enroot start --root --rw "$CONTAINER" bash -eux <<'INNER'
export DEBIAN_FRONTEND=noninteractive
apt-get update
apt-get install -y --no-install-recommends \
    binutils build-essential ca-certificates g++ cmake curl git libc-dev \
    liblapack-dev libopenblas64-serial-dev libopenmpi-dev libtool locales \
    ninja-build time tzdata wget pkg-config libgmp-dev libsuitesparse-dev \
    gfortran
INNER

echo "[build-1] --- 3/3 fetch sources: coinbrew, MUMPS, Ipopt, Bonmin, AlgoHex ---"
enroot start --root --rw --env ALGOHEX_COMMIT="$ALGOHEX_COMMIT" \
    "$CONTAINER" bash -eux <<'INNER'
mkdir -p /usr/src/coin-or /opt/coin-or
cd /usr/src/coin-or
wget -q https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
chmod +x coinbrew
# Upstream pins MUMPS and Ipopt but takes Bonmin from master -- the one moving
# part in this chain, and the first suspect if a future build stops matching.
./coinbrew fetch https://github.com/coin-or-tools/ThirdParty-Mumps@3.0.11
./coinbrew fetch Ipopt@3.14.19 --skip-update
./coinbrew fetch Bonmin@master --skip-update
# AlgoHex itself, at the pinned commit, with submodules -- it has .gitmodules
# and a shallow clone without them fails at cmake, not at ninja.
git clone https://github.com/cgg-bern/AlgoHex.git /app
cd /app
git checkout "$ALGOHEX_COMMIT"
git submodule update --init --recursive
git log --oneline -1
INNER

cat <<EOF
[build-1] done. Sources are in the container, nothing is compiled yet.

Next, on a COMPUTE node (this is the part that needs cores):
  sbatch --test-only experimentell/hex3d_algohex/scripts/build_algohex.slurm
  sbatch            experimentell/hex3d_algohex/scripts/build_algohex.slurm
EOF
