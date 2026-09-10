#!/bin/bash
# Smoke test 2 of 2: ONE sample end to end, at -n 2000. The first job that
# exercises the Python side, so it is where a missing module shows up.
#
#   sbatch --test-only experimentell/hex3d_algohex/scripts/one_sample.slurm
#   sbatch            experimentell/hex3d_algohex/scripts/one_sample.slurm
#   CAND=cand_003 sbatch .../one_sample.slurm      # a different candidate
#
# Runs tet_prep -> AlgoHex -> clean_blocks -> tfi -> export_sample, with
# reattach skipped: it is the only stage that needs torch, and it produces the
# full.* CFD variant decision J deletes anyway.
#
# Time budget, measured on the 2-vCPU VPS at -n 2000, so treat as an upper
# bound on a real node: AlgoHex 16 min, clean_blocks 1-3 min, tfi and
# export_sample under 1 min each, tet_prep ~3 min. 30 min of walltime.
#
# ---------------------------------------------------------------------------
# UNVERIFIED, CHECK BEFORE THE FIRST SUBMIT:
#   --partition   dev_cpu_il, unverified. sinfo -o "%P %l %c %m %D"
#   --time        30 min; the dev partition on bwUniCluster caps at 30
#   --mem         8G. AlgoHex peaks at ~6 GB and two concurrent runs OOMed a
#                 7.7 GiB box, which is what killed v7 and v8. Do not lower.
#   PYTHON_MODULE the system python3 is 3.9 and this needs >= 3.11
# ---------------------------------------------------------------------------
#SBATCH --job-name=algohex_one
#SBATCH --output=algohex_one_%j.out
#SBATCH --partition=dev_cpu_il
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G

set -u

WS="${WS:-$(ws_find hex3d 2>/dev/null)}"
REPO="${REPO:-$WS/domain_partition_3D}"
CAND="${CAND:-cand_001}"
PYTHON_MODULE="${PYTHON_MODULE:-devel/python/3.13.3-gnu-14.2}"

module load "$PYTHON_MODULE" 2>/dev/null || echo "[one] module load failed, using \$PATH python"
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
mkdir -p "$ENROOT_DATA_PATH" "$ENROOT_TEMP_PATH" "$ENROOT_RUNTIME_PATH"

# Unpack once into the persistent store, then pass the container NAME. For a
# single job the .sqsh would do, but the array must not re-read the squashfs
# through Lustre per task: the sibling project measured 20 of 30 budgeted
# minutes as pure IO wait unpacking 9.6 GB that way (job 6822816, 20 s CPU).
# Ours is 2.1 GB, so scale that down -- but pay it once, here.
if ! enroot list 2>/dev/null | grep -qxF algohex; then
    echo "[one] unpacking algohex.sqsh into $ENROOT_DATA_PATH (once, ever)"
    enroot create "$WS/enroot-images/algohex.sqsh"
else
    echo "[one] container 'algohex' already unpacked"
fi

cd "$REPO"
echo "[one] host $(hostname), job ${SLURM_JOB_ID:-none}, candidate $CAND"
python -c "import numpy, scipy, meshio; print('[one] numpy', numpy.__version__,
'scipy', scipy.__version__, 'meshio', meshio.__version__)" || {
    echo "[one] FAIL: the Python environment is incomplete."
    echo "[one] Need numpy, scipy, meshio (+ gmsh for tet_prep). NOT torch."
    exit 1; }

t0=$(date +%s)
BACKEND=enroot ALGOHEX_IMAGE=algohex ALGOHEX_N=2000 \
SKIP_REATTACH=1 TARGET_H=0.05 \
    bash experimentell/hex3d_algohex/scripts/generate_dataset.sh \
        "data/random_tistos/investigated/$CAND"
rc=$?
echo "[one] pipeline exit $rc in $(( $(date +%s) - t0 ))s"

d="output/hex3d_algohex/gen/$CAND"
echo "[one] --- what landed ---"
ls -la "$d" 2>/dev/null | tail -12
echo "[one] --- the artifact that matters ---"
if [ -f "$d/${CAND}_sample.npz" ]; then
    python - "$d/${CAND}_sample.npz" <<'PYEOF'
import sys, json
sys.path.insert(0, "experimentell/hex3d_algohex")
import export_sample as ex
s = ex.load(sys.argv[1]); q = s["quality"]
print(f"[one] PASS  {len(s['blocks'])} blocks, {len(s['edges'])} directed edges, "
      f"{len(s['quad_faces'])} shell faces")
print(f"[one]       fit residual median {q['fit_residual_median']*100:.3f} %, "
      f"max {q['fit_residual_max']*100:.2f} %, "
      f"arc/chord max {q['arc_over_chord_max']:.3f}")
print(f"[one]       provenance: {json.dumps(s['provenance'])[:160]}")
PYEOF
else
    echo "[one] FAIL: no sample file. Check $d/*.log, newest first."
fi
echo "[one] Next: T10's array, once this passes."
exit $rc
