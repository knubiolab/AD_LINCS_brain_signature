#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 03 - PBC correction and global stability metrics
#
#   bash 03_postprocess_trajectory.sh <system_dir> [more dirs ...]
#
# Produces, in each <system_dir>:
#   MD_center_300.xtc       PBC-corrected, protein-centred, 100-ps subsampled
#   rmsd_backbone.xvg       protein backbone RMSD (LSQ fit on backbone)
#   rmsd_backbone_LIG.xvg   ligand heavy-atom RMSD after backbone fit
#   rmsf.xvg                per-residue C-alpha RMSF
#   gyrate1.xvg             radius of gyration (protein)
#   hb.xvg                  protein-ligand hydrogen bonds
#   hb_essential.xvg        hydrogen bonds to the target's key residues
#   sasa_prot/lig/complex.xvg
#
# Hydrogen bonds use the GROMACS default criteria (donor-acceptor <= 0.35 nm,
# angle >= 150 deg). SASA uses the default 0.14 nm probe.
# ---------------------------------------------------------------------------
set -uo pipefail
GMX="${GMX:-gmx}"
HERE="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=targets.conf
source "$HERE/targets.conf"
DT_PS="${DT_PS:-100}"          # subsampling interval for post-processing

run() { printf '%s\n' "$@" | "$GMX" "${CMD[@]}" ; }

for SYS in "$@"; do
    [[ -d "$SYS" ]] || { echo "skip (no such directory): $SYS" >&2; continue; }
    echo "=== $SYS"
    cd "$SYS" || continue
    TPR=MD.tpr

    # --- PBC correction + subsampling --------------------------------------
    if [[ ! -f MD_center_300.xtc ]]; then
        # prompts: centre on Protein (1), output System (0)
        printf '1\n0\n' | "$GMX" trjconv -s "$TPR" -f MD.xtc -o MD_center_300.xtc \
            -center -pbc mol -ur compact -dt "$DT_PS"
    fi
    XTC=MD_center_300.xtc

    # --- index groups: ligand heavy atoms, key residues ---------------------
    if [[ ! -f index_MD.ndx ]]; then
        LIG="$(grep -m1 -oE '^[A-Z0-9]{2,4}' <<< "${LIGRES:-UNL}")"
        printf 'r %s & !a H*\nname 21 LIG_HEAVY\nq\n' "${LIGRES:-UNL}" \
            | "$GMX" make_ndx -f "$TPR" -o index_MD.ndx
    fi

    # --- RMSD --------------------------------------------------------------
    # ligand heavy-atom RMSD, least-squares fit on the protein backbone
    printf '4\n21\n' | "$GMX" rms -s "$TPR" -f "$XTC" -n index_MD.ndx \
        -o rmsd_backbone_LIG.xvg -tu ns
    # protein backbone RMSD
    printf '4\n4\n'  | "$GMX" rms -s "$TPR" -f "$XTC" -n index_MD.ndx \
        -o rmsd_backbone.xvg -tu ns

    # --- RMSF (per residue, C-alpha) ---------------------------------------
    echo 3 | "$GMX" rmsf -s "$TPR" -f "$XTC" -o rmsf.xvg -res

    # --- radius of gyration -------------------------------------------------
    echo 1 | "$GMX" gyrate -s "$TPR" -f "$XTC" -o gyrate1.xvg

    # --- hydrogen bonds -----------------------------------------------------
    # protein (1) vs ligand (13)
    printf '1\n13\n' | "$GMX" hbond -s "$TPR" -f "$XTC" -num hb.xvg -tu ns
    # key residues vs ligand; index_hb_essential.ndx is built by 04_site_contacts.sh
    if [[ -f index_hb_essential.ndx ]]; then
        printf 'hb_essential\n%s\n' "${LIGRES:-UNL}" \
            | "$GMX" hbond -s "$TPR" -f "$XTC" -n index_hb_essential.ndx \
                     -num hb_essential.xvg -tu ns
    else
        echo "    note: index_hb_essential.ndx absent; run 04_site_contacts.sh first"
    fi

    # --- SASA ---------------------------------------------------------------
    "$GMX" sasa -s "$TPR" -f "$XTC" -surface 'group "Protein"'             -o sasa_prot.xvg    -tu ns
    "$GMX" sasa -s "$TPR" -f "$XTC" -surface "resname ${LIGRES:-UNL}"      -o sasa_lig.xvg     -tu ns
    "$GMX" sasa -s "$TPR" -f "$XTC" \
        -surface "group \"Protein\" or resname ${LIGRES:-UNL}"             -o sasa_complex.xvg -tu ns

    cd - >/dev/null
done
echo "03 done"
