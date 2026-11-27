#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 06 - MM/GBSA binding free energies
#
#   bash 06_mmgbsa.sh <system_dir> [more dirs ...]
#
# Single-trajectory MM/GBSA with gmx_MMPBSA (v1.5.0.3) on the equilibrated
# phase: 401 evenly spaced frames over 200-300 ns (0.25 ns spacing).
# Generalised Born: OBC2 (igb = 5), 0.1 M salt. No entropy term is computed.
#
# Produces, in <system_dir>/mmgbsa/:
#   FINAL_RESULTS_MMPBSA.dat   per-system decomposition and total dG
#   mmgbsa.csv                 machine-readable summary
#
# gmx_MMPBSA needs the ligand parameters that were used for the simulation, so
# it is run from the directory holding topol.top and the included .itp files.
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=targets.conf
source "$HERE/targets.conf"
BEG_PS=$((WINDOW_BEGIN_NS * 1000))
END_PS=$((WINDOW_END_NS * 1000))
STRIDE_PS="${STRIDE_PS:-250}"               # 0.25 ns -> 401 frames over 100 ns
NPROC="${NPROC:-8}"

read -r -d '' MMPBSA_IN <<EOF
Sampling the equilibrated phase of the production run
&general
  sys_name     = "complex",
  startframe   = 1,
  endframe     = 401,
  interval     = 1,
  verbose      = 2,
  keep_files   = 0,
/
&gb
  igb          = 5,
  saltcon      = 0.100,
/
EOF

for SYS in "$@"; do
    [[ -d "$SYS" ]] || { echo "skip: $SYS" >&2; continue; }
    echo "=== $SYS"
    cd "$SYS" || continue
    [[ -f topol.top ]] || { echo "    topol.top not found" >&2; cd - >/dev/null; continue; }
    mkdir -p mmgbsa
    printf '%s\n' "$MMPBSA_IN" > mmgbsa/mmpbsa.in

    # frames of the analysis window at the chosen stride
    gmx trjconv -s MD.tpr -f MD_center_300.xtc -o mmgbsa/window.xtc \
        -b "$BEG_PS" -e "$END_PS" -dt "$STRIDE_PS" <<< $'0\n' >/dev/null 2>&1

    # receptor = Protein, ligand = the SwissParam residue
    ( cd mmgbsa && mpirun -np "$NPROC" gmx_MMPBSA MPI \
        -O -i mmpbsa.in -cs ../MD.tpr -ci ../index.ndx -cg 1 13 \
        -ct window.xtc -cp ../topol.top \
        -o FINAL_RESULTS_MMPBSA.dat -eo mmgbsa.csv > mmgbsa.log 2>&1 )

    if grep -qE '^ *DELTA TOTAL' mmgbsa/FINAL_RESULTS_MMPBSA.dat 2>/dev/null; then
        awk '/^ *DELTA TOTAL/ {printf "    dG = %s +/- %s kcal/mol\n", $3, $4}' \
            mmgbsa/FINAL_RESULTS_MMPBSA.dat
    else
        echo "    failed; see $SYS/mmgbsa/mmgbsa.log" >&2
    fi
    cd - >/dev/null
done
echo "06 done"
