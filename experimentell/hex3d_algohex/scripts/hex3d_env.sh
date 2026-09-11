#!/usr/bin/env bash
# Per-login setup for the hex3d pipeline on bwUniCluster. SOURCE it, do not
# execute it:
#
#     source experimentell/hex3d_algohex/scripts/hex3d_env.sh
#
# Add that line to ~/.bashrc and every login is ready. Idempotent.
#
# WHY THIS EXISTS. The submit scripts set every ENROOT_* variable themselves,
# so they work. A shell you type in has none of them, and then enroot fails in
# two ways that look like broken containers and are not:
#
#   enroot start --root "$ENROOT_IMAGE" ...
#     [ERROR] Invalid argument              <- the variable is empty
#
#   enroot start --root dtOO ...
#     [ERROR] No such file or directory: ~/.local/share/enroot/dtOO
#                                           <- ENROOT_DATA_PATH unset, so
#                                              enroot looked in $HOME
#
# Modelled on eigenfrequencies/cluster/cluster_env.sh, which solves the same
# problem for that project.

# ── workspace ─────────────────────────────────────────────────────────────
if [ -z "${WS:-}" ]; then
    WS="$(ws_find "${WS_NAME:-hex3d}" 2>/dev/null)"
fi
if [ -z "$WS" ]; then
    echo "hex3d_env: no workspace '${WS_NAME:-hex3d}' — run: ws_allocate hex3d 60" >&2
else
    export WS
    export REPO="${REPO:-$WS/domain_partition_3D}"
fi

# ── enroot ────────────────────────────────────────────────────────────────
# Shared across nodes and jobs: this is what makes "unpack once, ever" work.
export ENROOT_DATA_PATH="${ENROOT_DATA_PATH:-$WS/enroot-data}"
export ENROOT_CACHE_PATH="${ENROOT_CACHE_PATH:-$WS/enroot-cache}"
# Node-local, both of them. TEMP because enroot flattens image layers through
# an overlayfs and a parallel filesystem cannot host one. RUNTIME because the
# default is ${XDG_RUNTIME_DIR}/enroot, i.e. /run/user/$UID, which exists on a
# login node and not on a compute node -- see docs/cluster-enroot-findings.md.
export ENROOT_TEMP_PATH="${ENROOT_TEMP_PATH:-/tmp/${USER:-$(id -un)}-enroot}"
export ENROOT_RUNTIME_PATH="${ENROOT_RUNTIME_PATH:-/tmp/${USER:-$(id -un)}-enroot-run}"
# zstd, never lzo: an lzo image imports without complaint and is then
# unreadable on every compute node.
export ENROOT_SQUASH_OPTIONS="${ENROOT_SQUASH_OPTIONS:--comp zstd -noD}"
mkdir -p "$ENROOT_TEMP_PATH" "$ENROOT_RUNTIME_PATH" 2>/dev/null

# Handy for typing by hand; the submit scripts set their own.
export ALGOHEX_SQSH="${ALGOHEX_SQSH:-$WS/enroot-images/algohex.sqsh}"

# ── python ────────────────────────────────────────────────────────────────
# The system python3 is 3.9 and the pipeline needs >= 3.11. The venv holds
# numpy/scipy/meshio/gmsh and deliberately NOT torch.
PYTHON_MODULE="${PYTHON_MODULE:-devel/python/3.13.3-gnu-14.2}"
if command -v module >/dev/null 2>&1; then
    module load "$PYTHON_MODULE" 2>/dev/null \
        || echo "hex3d_env: could not load $PYTHON_MODULE" >&2
fi
if [ -n "${REPO:-}" ] && [ -x "$REPO/.venv/bin/python" ]; then
    export PATH="$REPO/.venv/bin:$PATH"
fi

# ── status, one line ──────────────────────────────────────────────────────
printf 'hex3d_env: WS=%s  image=%s  containers=%s\n' \
    "${WS:-<none>}" \
    "$([ -f "$ALGOHEX_SQSH" ] && echo present || echo MISSING)" \
    "$(enroot list 2>/dev/null | tr '\n' ' ')"
