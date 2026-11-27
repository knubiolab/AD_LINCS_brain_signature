#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 05 - Conformational clustering of the ligand and representative poses
#
#   bash 05_cluster_representative_pose.sh <system_dir> [more dirs ...]
#
# Clusters the ligand over the equilibrated phase (200-300 ns) with the GROMOS
# algorithm at a 0.10 nm RMSD cutoff, after superposition on the protein
# backbone, and writes the representative structure of every cluster needed to
# reach 50% cumulative occupancy. Each representative pose is then profiled
# with PLIP.
#
# Produces, in <system_dir>/rep_pose/:
#   fitted.xtc          backbone-fitted trajectory of the analysis window
#   clusters.pdb        representative structure of each cluster
#   cluster.log         gmx cluster report
#   clust_size.xvg      cluster sizes  (input for plot_cluster_fraction.R)
#   clust_id.xvg        per-frame cluster assignment
#   rep_cluster<N>.pdb  representative pose of cluster N (complex)
#   plip_rank<N>.txt/.xml/.pse
#
# Note on the gmx cluster prompts: because the trajectory is pre-fitted and
# -nofit is used, gmx cluster asks for exactly two groups -- the group the RMSD
# is computed on (ligand heavy atoms) and the group written to the output
# structures (protein + ligand).
# ---------------------------------------------------------------------------
set -uo pipefail
GMX="${GMX:-gmx}"
HERE="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=targets.conf
source "$HERE/targets.conf"
BEG_PS=$((WINDOW_BEGIN_NS * 1000))
END_PS=$((WINDOW_END_NS * 1000))
CUM_PCT="${CUM_PCT:-50}"          # cumulative occupancy to cover

for SYS in "$@"; do
    [[ -d "$SYS" ]] || { echo "skip: $SYS" >&2; continue; }
    echo "=== $SYS"
    cd "$SYS" || continue
    TPR=MD.tpr; XTC=MD_center_300.xtc; LIG="${LIGRES:-UNL}"
    D=rep_pose; mkdir -p "$D"

    # --- selection groups ---------------------------------------------------
    "$GMX" select -s "$TPR" -on "$D/rep.ndx" -select \
        "\"fit\" group \"Backbone\"; \"lig\" resname $LIG; \
         \"complex\" group \"Protein\" or resname $LIG; \"all\" all" >/dev/null 2>&1

    # --- superpose on the protein backbone, restricted to the window --------
    printf 'fit\nall\n' | "$GMX" trjconv -s "$TPR" -f "$XTC" -n "$D/rep.ndx" \
        -fit rot+trans -b "$BEG_PS" -e "$END_PS" -o "$D/fitted.xtc" >/dev/null 2>&1

    # --- cluster ------------------------------------------------------------
    printf 'lig\ncomplex\n' | "$GMX" cluster -s "$TPR" -f "$D/fitted.xtc" -nofit \
        -n "$D/rep.ndx" -method gromos -cutoff "$CLUSTER_CUTOFF_NM" \
        -cl "$D/clusters.pdb" -g "$D/cluster.log" -o "$D/rmsd_clust.xpm" \
        -dist "$D/rmsd_dist.xvg" -sz "$D/clust_size.xvg" -clid "$D/clust_id.xvg" >/dev/null 2>&1

    # --- how many clusters are needed to reach CUM_PCT ----------------------
    NEED=$(awk -v thr="$CUM_PCT" '!/^[#@]/ {n[++k]=$2; tot+=$2}
        END{c=0; for(i=1;i<=k;i++){c+=n[i]; if(100*c/tot>=thr){print i; exit}} print k}' \
        "$D/clust_size.xvg")
    TOP1=$(awk '!/^[#@]/ {n[++k]=$2; tot+=$2} END{printf "%.1f", 100*n[1]/tot}' "$D/clust_size.xvg")
    echo "    clusters to reach ${CUM_PCT}%: $NEED   (largest cluster ${TOP1}%)"

    # --- split clusters.pdb into one file per representative pose -----------
    python3 - "$D/clusters.pdb" "$D" "$NEED" <<'PY'
import sys
src, out, need = sys.argv[1], sys.argv[2], int(sys.argv[3])
models, cur = [], []
for line in open(src):
    if line.startswith("MODEL"):
        cur = []
    elif line.startswith("ENDMDL"):
        models.append(cur)
    elif line.startswith(("ATOM", "HETATM", "TER")):
        cur.append(line)
if not models and cur:
    models = [cur]
for i in range(min(need, len(models))):
    with open(f"{out}/rep_cluster{i+1}.pdb", "w") as fh:
        fh.writelines(models[i]); fh.write("END\n")
print(f"    wrote {min(need, len(models))} representative pose(s)")
PY

    # --- PLIP ---------------------------------------------------------------
    if command -v plip >/dev/null 2>&1; then
        for i in $(seq 1 "$NEED"); do
            [[ -f "$D/rep_cluster${i}.pdb" ]] || continue
            ( cd "$D" && plip -f "rep_cluster${i}.pdb" -t -x -y \
                 --name "plip_rank${i}" >/dev/null 2>&1 )
        done
    else
        echo "    note: plip not on PATH; skipping interaction profiling"
    fi

    cd - >/dev/null
done
echo "05 done"
