#!/bin/bash
#SBATCH --job-name=dtoo_bisect3
#SBATCH --partition=dev_cpu_il
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=dtoo_bisect3_%j.out
#
# Bisect3 -- job 6869858 (uc2n601) sat 18+ min with ZERO file activity and no
# container output. It was the first run with REAL two-field enroot mounts
# (all earlier runs carried the broken ":rw" fstype and bound nothing), and
# the first run that could actually reach generate_machines.py. The steps
# below isolate the two suspects:
#
#   P1   real mount, remapped form $REPO:/repo -- does a real rbind hang?
#   P1b  real mount, production form (same path on both sides, mirrors the
#        proven eigenfrequencies invocation)
#   P2   like P1 but WITHOUT "< /dev/null" (the batch script omits it)
#   P3   source OpenFOAM + dtOO env (proven healthy in bisect2)
#   P4a  container python packages (yaml, pydantic)
#   P4b  import eigenfrequencies.config
#   P4c  import eigenfrequencies.bc.builders
#   P4d  import the DtooAdapter (line ~131 of generate_machines.py)
#        -- the whole chain runs BEFORE the script's first print, so a
#           block here reproduces the 6869858 silence exactly
#   P5   single-machine export through DtooAdapter.export_mesh, the same
#        call generate_machines.py makes, but for one machine only
#   P6   the exact batch payload: generate_machines.py --export in the form
#        job 6869858 ran and never left a trace of
#
# Every step: "timeout T srun ... < /dev/null", rc + timestamp echoed.
#
# v2 fix: the mount array must carry the --mount flags themselves -- v1
# sent bare specs to enroot ("--rc ... /pfs/...:/repo ..."), enroot read
# them as stray arguments and died with "[ERROR] Invalid argument: ...".

set -u

module load devel/python/3.13.3-gnu-14.2 || true

WS=$(ws_find hex3d)
REPO="${REPO:-${SLURM_SUBMIT_DIR:-$WS/domain_partition_3D}}"
EIG_WS=$(ws_find eigenfreq)

export ENROOT_TEMP_PATH="${ENROOT_TEMP_PATH:-/tmp/$USER-enroot}"
export ENROOT_RUNTIME_PATH="${ENROOT_RUNTIME_PATH:-/tmp/$USER-enroot-run}"
export ENROOT_DATA_PATH="${ENROOT_DATA_PATH:-$EIG_WS/enroot-data}"

SQSH="$EIG_WS/enroot-images/dtOO.sqsh"
CONTAINER=dtOO

EF=""
for c in "$EIG_WS/eigenfrequencies" "$HOME/eigen/eigenfrequencies" "$HOME/eigenfrequencies"; do
    if [ -f "$c/src/eigenfrequencies/adapters/dtoo/adapter.py" ]; then EF="$c"; break; fi
done
[ -n "$EF" ] || { echo "[bisect3] no eigenfrequencies checkout"; exit 1; }
RC="$EF/cluster/enroot_rc.sh"
[ -f "$RC" ] || { echo "[bisect3] no --rc script: $RC"; exit 1; }
[ -f "$SQSH" ] || { echo "[bisect3] no image: $SQSH"; exit 1; }

LOGDIR="$REPO/experimentell/hex3d_algohex/logs"
mkdir -p "$LOGDIR"
LOG="$LOGDIR/dtoo_bisect3_${SLURM_JOB_ID:-manual}.log"
exec > >(tee "$LOG") 2>&1

echo "[bisect3] host $(hostname), job ${SLURM_JOB_ID:-?}, $(date +%T)"
echo "[bisect3] repo $REPO"
echo "[bisect3] ef   $EF"
echo "[bisect3] img  $SQSH"
echo "[bisect3] data $ENROOT_DATA_PATH"

run_step() {
    local t="$1" label="$2"; shift 2
    echo "[bisect3] --- $label --- $(date +%T)"
    timeout "$t" srun -n1 -N1 enroot start --root "$@" < /dev/null
    local rc=$?
    echo "[bisect3] --- $label rc=$rc --- $(date +%T)"
    return "$rc"
}

py_step() {
    local t="$1" label="$2" code="$3"
    run_step "$t" "$label" --rc "$RC" "${MS[@]}" "$CONTAINER" \
        bash -c "source /usr/lib/openfoam/openfoam2606/etc/bashrc; source /dtOO-install/bin/env.sh; $code"
}

# ------------------------------------------------------------------ mounts
run_step 90 "P1 real mounts remapped" --rc "$RC" \
    --mount "$REPO:/repo" --mount "$EF:/ef" "$CONTAINER" \
    bash -c 'date +%T; pwd; for i in 1 2 3; do echo "MARK$i"; done; ls -ld /repo /ef; ls /repo/data/dataset/sobol | head -3; exec echo P1-OK'
R1=$?

run_step 90 "P1b real mounts same-path (production form)" --rc "$RC" \
    --mount "$REPO:$REPO" --mount "$EF:$EF" "$CONTAINER" \
    bash -c "date +%T; ls -ld $REPO $EF; exec echo P1B-OK"
R1B=$?

echo "[bisect3] --- P2 start WITHOUT < /dev/null (batch form) --- $(date +%T)"
timeout 90 srun -n1 -N1 enroot start --root --rc "$RC" \
    --mount "$REPO:/repo" --mount "$EF:/ef" "$CONTAINER" bash -c 'exec echo P2-OK'
R2=$?
echo "[bisect3] --- P2 rc=$R2 --- $(date +%T)"

# ------------------------------------------------------------ python chain
if [ "$R1" -eq 0 ]; then
    MS=(--mount "$REPO:/repo" --mount "$EF:/ef"); R_IN=/repo; E_IN=/ef
elif [ "$R1B" -eq 0 ]; then
    MS=(--mount "$REPO:$REPO" --mount "$EF:$EF"); R_IN="$REPO"; E_IN="$EF"
else
    MS=(--mount "$REPO:/repo" --mount "$EF:/ef"); R_IN=/repo; E_IN=/ef
    echo "[bisect3] no mount form worked -- python steps will be skipped"
fi

R3=n/a; R4A=n/a; R4B=n/a; R4C=n/a; R4D=n/a; R5=n/a; R6=n/a
if [ "$R1" -eq 0 ] || [ "$R1B" -eq 0 ]; then
    run_step 90 "P3 sources" --rc "$RC" "${MS[@]}" "$CONTAINER" \
        bash -c 'source /usr/lib/openfoam/openfoam2606/etc/bashrc; source /dtOO-install/bin/env.sh; date +%T; exec echo P3-OK'
    R3=$?

    py_step 90 "P4a packages" \
        "python3.13 -c 'import sys; print(sys.version.split()[0]); import yaml, pydantic; print(\"P4A-PKG-OK\")'"
    R4A=$?
    py_step 120 "P4b config" \
        "python3.13 -c 'import sys; sys.path.insert(0, \"$E_IN/src\"); import eigenfrequencies.config; print(\"P4B-CFG-OK\")'"
    R4B=$?
    py_step 120 "P4c bc.builders" \
        "python3.13 -c 'import sys; sys.path.insert(0, \"$E_IN/src\"); import eigenfrequencies.bc.builders; print(\"P4C-BC-OK\")'"
    R4C=$?
    py_step 240 "P4d adapter" \
        "python3.13 -c 'import sys, time; sys.path.insert(0, \"$E_IN/src\"); t = time.time(); from eigenfrequencies.adapters.dtoo.adapter import DtooAdapter; print(\"P4D-ADAPTER-OK\", round(time.time() - t, 1))'"
    R4D=$?

    py_step 300 "P5 single export" \
        "python3.13 -u -c 'import json, os, sys, time; sys.path.insert(0, \"$E_IN/src\"); from eigenfrequencies.adapters.dtoo.adapter import DtooAdapter; a = DtooAdapter(\"$E_IN/adapters/machines/tistos.yaml\"); v = json.load(open(\"$R_IN/data/dataset/sobol/machine_0001/params.json\")); out = \"$R_IN/data/dataset/sobol/machine_0001/mesh.msh\"; t = time.time(); a.export_mesh(v, output_msh=out + \".part\"); os.replace(out + \".part\", out); print(\"P5-EXPORT-OK\", round(time.time() - t, 1))'"
    R5=$?
    ls -l "$REPO/data/dataset/sobol/machine_0001/" | sed 's/^/[bisect3] /'

    py_step 240 "P6 batch payload (generate_machines --export)" \
        "cd $R_IN; python3.13 -u experimentell/hex3d_algohex/scripts/generate_machines.py --export --count 64"
    R6=$?
fi

echo
echo "[bisect3] summary: P1=$R1 P1b=$R1B P2=$R2 P3=$R3 P4a=$R4A P4b=$R4B P4c=$R4C P4d=$R4D P5=$R5 P6=$R6"
echo "[bisect3] log $LOG"
echo "[bisect3] commit it: git add experimentell/hex3d_algohex/logs && git commit -m 'cluster: dtoo bisect3 log'"
wait
