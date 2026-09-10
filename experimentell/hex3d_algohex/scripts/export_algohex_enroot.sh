#!/usr/bin/env bash
# Export the AlgoHex Docker image as an enroot squashfs for bwUniCluster 3.0.
#
#   bash experimentell/hex3d_algohex/scripts/export_algohex_enroot.sh [OUT_DIR]
#
# Default OUT_DIR: output/hex3d_algohex, so the .sqsh sits where
# run_algohex.py --backend enroot looks for it by default.
#
# Modelled on eigenfrequencies/cluster/export_dtoo_enroot.sh, with one
# deliberate difference: that script produces a docker-save tarball for the
# cluster to import, this one produces the .sqsh directly. Two reasons.
# `enroot import dockerd://` needs no intermediate 2.7 GB tar (this box has
# 10 GiB free), and shipping the .sqsh means the cluster does no import work at
# all -- on Lustre at the 5-10 MB/s measured in
# eigenfrequencies/docs/cluster-dtoo-enroot-befund-v2.md, import is the
# expensive step, not the transfer.
#
# THE TRAP THAT MATTERS, inherited from cluster/enroot_dtoo_import.md:
# ENROOT_SQUASH_OPTIONS must ask for zstd. Without it mksquashfs defaults to
# lzo, and an lzo image imports without complaint and is then UNREADABLE ON
# EVERY COMPUTE NODE. It is set below rather than left to the environment.

set -euo pipefail

IMAGE="${ALGOHEX_IMAGE:-algohex:portable}"
REPO=/root/repos/duty/quadmesh/domain_partition_3D
OUT_DIR="${1:-$REPO/output/hex3d_algohex}"
SQSH="$OUT_DIR/algohex.sqsh"

export ENROOT_SQUASH_OPTIONS="-comp zstd -noD"
export ENROOT_DATA_PATH="${ENROOT_DATA_PATH:-/root/enroot/data}"
export ENROOT_CACHE_PATH="${ENROOT_CACHE_PATH:-/root/enroot/cache}"
export ENROOT_TEMP_PATH="${ENROOT_TEMP_PATH:-/root/enroot/tmp}"
mkdir -p "$OUT_DIR" "$ENROOT_DATA_PATH" "$ENROOT_CACHE_PATH" "$ENROOT_TEMP_PATH"

echo "[export] image        $IMAGE"
echo "[export] out          $SQSH"
echo "[export] squash opts  $ENROOT_SQUASH_OPTIONS"

docker image inspect "$IMAGE" >/dev/null 2>&1 || {
    echo "[export] ERROR: $IMAGE not in the local docker daemon."
    echo "[export] Build it from external_patches/Dockerfile.portable first."
    exit 1
}

# Self-test before spending 20 minutes on a squashfs: an image whose binary
# cannot run is exactly the failure T2 exists to prevent.
docker run --rm "$IMAGE" HexMeshing --help >/dev/null 2>&1 || {
    echo "[export] ERROR: '$IMAGE HexMeshing --help' fails with no volume."
    echo "[export] That is T2's gate. Do not ship this image."
    exit 1
}
echo "[export] self-test OK: HexMeshing runs from the image, no volume"

df -h "$OUT_DIR" | tail -1

rm -f "$SQSH"
# dockerd:// and not docker://: this box has IPv4 forwarding disabled, so the
# registry route cannot resolve DNS (see /root/enroot/import.log). The image is
# already local anyway.
enroot import -o "$SQSH" "dockerd://$IMAGE"

ls -la "$SQSH"
echo "[export] compression check (must say zstd, NOT lzo):"
unsquashfs -s "$SQSH" 2>/dev/null | grep -i compress || true

cat <<EOF
[export] done.

Cluster side, with cluster_env.sh sourced so \$ENROOT_IMAGES is set:
  1. scp $SQSH bwunicluster:"\$ENROOT_IMAGES/"
  2. no import needed -- run_algohex.py --backend enroot starts straight from
     the .sqsh, so there is no 'enroot create' and no rootfs unpack. That is
     the step that costs 14-27 min per job on Lustre.
  3. smoke test one job:
     enroot start --root --rw --mount \$PWD:/work \\
         "\$ENROOT_IMAGES/algohex.sqsh" HexMeshing --help
EOF
