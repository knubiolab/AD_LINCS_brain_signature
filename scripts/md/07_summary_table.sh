#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 07 - Stability summary across systems and replicates
#
#   bash 07_summary_table.sh <system_dir> [more dirs ...] > stability_summary.tsv
#
# Reports, per system, over the equilibrated phase (200-300 ns):
#   ligand heavy-atom RMSD   mean +/- SD
#   protein backbone RMSD    mean +/- SD
#   per-residue RMSF         mean / maximum
#   MM/GBSA dG               mean +/- SD across replicates
#
# This reproduces Supplementary Table 3 of the manuscript.
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=targets.conf
source "$HERE/targets.conf"
BEG="$WINDOW_BEGIN_NS"; END="$WINDOW_END_NS"

# mean and SD of column 2 of an .xvg, restricted to a time window (x in ns)
stat_window() {          # $1 = file, $2 = begin, $3 = end
    [[ -f "$1" ]] || { printf 'NA\tNA'; return; }
    awk -v b="$2" -v e="$3" '!/^[#@]/ && $1>=b && $1<=e {
            s+=$2; ss+=$2*$2; n++ }
        END { if(!n) { printf "NA\tNA" }
              else { m=s/n; v=ss/n-m*m; if(v<0) v=0; printf "%.3f\t%.3f", m, sqrt(v) } }' "$1"
}
# mean and max of column 2 over all rows (RMSF is indexed by residue, not time)
stat_all() {
    [[ -f "$1" ]] || { printf 'NA\tNA'; return; }
    awk '!/^[#@]/ { s+=$2; n++; if($2>mx) mx=$2 }
         END { if(!n) printf "NA\tNA"; else printf "%.3f\t%.3f", s/n, mx }' "$1"
}
# MM/GBSA total from gmx_MMPBSA output
mmgbsa_dg() {
    local f="$1/mmgbsa/FINAL_RESULTS_MMPBSA.dat"
    [[ -f "$f" ]] || { printf 'NA'; return; }
    awk '/^ *DELTA TOTAL/ { printf "%s", $3; exit }' "$f"
}

printf 'system\tligand_RMSD_mean\tligand_RMSD_sd\tbackbone_RMSD_mean\tbackbone_RMSD_sd\tRMSF_mean\tRMSF_max\tMMGBSA_dG\n'
for SYS in "$@"; do
    [[ -d "$SYS" ]] || continue
    printf '%s\t' "$SYS"
    stat_window "$SYS/rmsd_backbone_LIG.xvg" "$BEG" "$END"; printf '\t'
    stat_window "$SYS/rmsd_backbone.xvg"     "$BEG" "$END"; printf '\t'
    stat_all    "$SYS/rmsf.xvg";                              printf '\t'
    mmgbsa_dg   "$SYS";                                       printf '\n'
done
