#!/bin/bash
# The whole chain on ONE input, at -n 2000, from the fixtures archive.
#
#   sbatch --test-only experimentell/hex3d_algohex/scripts/smoke_pipeline.slurm
#   sbatch            experimentell/hex3d_algohex/scripts/smoke_pipeline.slurm
#
# AlgoHex -> clean_blocks -> tfi -> export_sample. Starts from
# data/T1_9/T1_9_tet_v5.vtk, which comes out of
# experimentell/hex3d_algohex/fixtures/blocks_core.tar.zst -- 5.5 MiB, IN GIT.
#
# So this test needs NO geometry transfer from anywhere: git clone plus the
# fixtures archive plus the image built by build_algohex.slurm. It is the
# cheapest end-to-end proof that the cluster can produce a sample.
#
# one_sample.slurm is the other variant: it starts at tet_prep from a
# candidate's mesh.msh, and those are 22 MB each and not in git.
#
# ---------------------------------------------------------------------------
# dev_cpu_il and its 30-minute CAP, deliberately: cpu_il has room but measures
# about a week of queue time, so a test that waits a week is not a test.
# Budget adds to ~20 min -- AlgoHex 16, clean_blocks 1-3, tfi and export under
# 1 each -- all upper bounds measured on a 2-vCPU VPS at -n 2000, so a 64-core
# node should be well inside. If it hits the wall anyway, drop to -n 500:
#   ALGOHEX_N=500 sbatch .../smoke_pipeline.slurm
# --mem 8G: AlgoHex peaks near 6 GB. Two concurrent runs OOMed a 7.7 GiB box.
# ---------------------------------------------------------------------------
#SBATCH --job-name=hex3d_smoke_pipe
#SBATCH --output=hex3d_smoke_pipe_%j.out
#SBATCH --partition=dev_cpu_il
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G

set -u

WS="${WS:-$(ws_find hex3d 2>/dev/null)}"
REPO="${REPO:-$WS/domain_partition_3D}"
N="${ALGOHEX_N:-2000}"
H="${TARGET_H:-0.05}"
TAG="smoke${N}"
PYTHON_MODULE="${PYTHON_MODULE:-devel/python/3.13.3-gnu-14.2}"

module load "$PYTHON_MODULE" 2>/dev/null || true
export PATH="$REPO/.venv/bin:$PATH"
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
mkdir -p "$ENROOT_TEMP_PATH" "$ENROOT_RUNTIME_PATH"

cd "$REPO" || exit 1
E=experimentell/hex3d_algohex
D=output/hex3d_algohex/$TAG
mkdir -p "$D"

echo "[smoke-pipe] host $(hostname), job ${SLURM_JOB_ID:-none}, -n $N"

# The input comes out of the committed fixtures archive.
TET=data/T1_9/T1_9_tet_v5.vtk
[ -f "$TET" ] || {
    echo "[smoke-pipe] restoring fixtures"
    tar -I zstd -xf $E/fixtures/blocks_core.tar.zst || exit 1; }
[ -f "$TET" ] || { echo "[smoke-pipe] FAIL: no $TET"; exit 1; }

python -c "import numpy, scipy, meshio" 2>/dev/null || {
    echo "[smoke-pipe] FAIL: need numpy, scipy, meshio in \$PATH python. NOT torch."
    exit 1; }

# The container: unpacked once into the persistent store, then used by NAME.
IMG=algohex
if ! enroot list 2>/dev/null | grep -qxF "$IMG"; then
    echo "[smoke-pipe] unpacking algohex.sqsh once"
    enroot create --name "$IMG" "$WS/enroot-images/algohex.sqsh" || exit 1
fi

# What smoke_algohex.slurm does as its own job, folded in here as five seconds
# of checks. dev_cpu_il allows only 4 queued jobs per user
# (QOSMaxSubmitJobPerUserLimit), and three of those are build stages, so a
# separate smoke job costs the slot this job needs.
echo "[smoke-pipe] --- 0/4 image usable on this compute node ---"
enroot start --root --rw --mount "$PWD:/work" "$IMG" \
    sh -c 'HexMeshing --help >/dev/null \
        && n=$(ldd /app/build/Build/bin/HexMeshing 2>/dev/null | grep -c "not found") \
        && echo "HexMeshing runs, $n missing libraries" && [ "$n" = 0 ]' \
  || { echo "[smoke-pipe] FAIL: the image does not work here."
       echo "[smoke-pipe] Run smoke_algohex.slurm alone for the detailed checks."
       exit 1; }

t0=$(date +%s)
echo "[smoke-pipe] --- 1/4 AlgoHex, enroot backend ---"
python $E/run_algohex.py --backend enroot --image "$IMG" \
    --prefix smoke --tag "$TAG" --in-vtk "$TET" --out-dir "$D" \
    -- -n "$N" || { echo "[smoke-pipe] FAIL: AlgoHex"; exit 1; }

echo "[smoke-pipe] --- 2/4 clean_blocks ---"
python $E/clean_blocks.py "$D/smoke_hex_${TAG}.ovm" --input-vtk "$TET" \
    --collapse-rounds 5 --untangle --untangle-rounds 6 \
    --out "$D/blocks.vtk" > "$D/clean_blocks.log" 2>&1 \
  || { echo "[smoke-pipe] FAIL: clean_blocks, see $D/clean_blocks.log"; exit 1; }
grep -E "^  [0-9]+ blocks,|Hausdorff|VALID|INVALID" "$D/clean_blocks.log" | tail -4

echo "[smoke-pipe] --- 3/4 conforming divisions ---"
python $E/tfi.py "$D/blocks.vtk" --input-vtk "$TET" --target-h "$H" \
    --require-lattices --solve-divisions > "$D/tfi.log" 2>&1
grep -E "direction classes|target h" "$D/tfi.log" | tail -2

echo "[smoke-pipe] --- 4/4 the sample, with the round-trip gate ---"
python $E/export_sample.py "$D/blocks.vtk" --tet "$TET" --target-h "$H" \
    --n "$N" --out "$D/sample.npz" --check
rc=$?

echo "[smoke-pipe] total $(( $(date +%s) - t0 ))s, exit $rc"
if [ $rc -eq 0 ]; then
    echo "[smoke-pipe] PASS -- the cluster produced a verified sample:"
    ls -la "$D/sample.npz"
    echo "[smoke-pipe] Next: transfer the candidate meshes, then one_sample.slurm"
else
    echo "[smoke-pipe] FAIL in export/round-trip. Newest log in $D."
fi
exit $rc
