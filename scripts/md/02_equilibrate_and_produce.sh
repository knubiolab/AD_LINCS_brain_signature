#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 02 - Equilibration and 300-ns production run
#
#   bash 02_equilibrate_and_produce.sh <system_dir> [seed]
#
# Stages (parameters are set in the .mdp files, not here):
#   NVT   100 ps, velocity-rescale thermostat, 300 K, tau_t = 0.1 ps,
#         position restraints on heavy atoms
#   NPT   1 ns, Berendsen barostat, 1 bar, tau_p = 2.0 ps, restraints retained
#   MD    300 ns, V-rescale thermostat (Protein_LIG vs solvent coupling) and
#         Parrinello-Rahman barostat, 2 fs step, LINCS on all bonds,
#         PME (1.2 nm cutoff, 0.16 nm grid), force-switch vdW 1.0-1.2 nm,
#         coordinates written every 10 ps
#
# Each system was run in triplicate from independent velocity seeds. Pass the
# seed as the second argument; it is written into gen_seed of the NVT run and
# the run is placed in <system_dir>/seed_<seed>/.
# ---------------------------------------------------------------------------
set -euo pipefail
GMX="${GMX:-gmx}"
SYS="${1:?usage: bash 02_equilibrate_and_produce.sh <system_dir> [seed]}"
SEED="${2:-}"
MDP="$(cd "$(dirname "$0")/mdp" && pwd)"
# Tune to the host; these were the settings used on a 2-GPU node.
MDRUN_OPTS="${MDRUN_OPTS:--ntmpi 2 -ntomp 32 -pin on -gpu_id 01}"

cd "$SYS"
SRC="$PWD"
if [[ -n "$SEED" ]]; then
    mkdir -p "seed_$SEED"
    for f in topol.top index.ndx posre_LIG.itp LIG.itp EM.gro *.ff; do
        [[ -e "$f" ]] && cp -r "$f" "seed_$SEED/" 2>/dev/null || true
    done
    cd "seed_$SEED"
    sed -e "s/^gen_seed.*/gen_seed = $SEED/" "$MDP/NVT.mdp" > NVT.mdp
else
    cp "$MDP/NVT.mdp" NVT.mdp
fi

# --- NVT -------------------------------------------------------------------
"$GMX" grompp -f NVT.mdp -c EM.gro -r EM.gro -p topol.top -n index.ndx -maxwarn 2 -o NVT.tpr
"$GMX" mdrun -deffnm NVT $MDRUN_OPTS

# --- NPT -------------------------------------------------------------------
"$GMX" grompp -f "$MDP/NPT.mdp" -c NVT.gro -r NVT.gro -p topol.top -n index.ndx -maxwarn 2 -o NPT.tpr
"$GMX" mdrun -deffnm NPT $MDRUN_OPTS

# --- production ------------------------------------------------------------
"$GMX" grompp -f "$MDP/MD.mdp" -c NPT.gro -t NPT.cpt -p topol.top -n index.ndx -maxwarn 2 -o MD.tpr
"$GMX" mdrun -deffnm MD $MDRUN_OPTS

# Resume an interrupted run with:
#   gmx mdrun -deffnm MD -cpi MD.cpt -append $MDRUN_OPTS

echo "02 done: $SYS${SEED:+ (seed $SEED)}"
