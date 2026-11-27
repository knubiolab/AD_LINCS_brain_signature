#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 04 - Per-residue contacts and binding-site distances
#
#   bash 04_site_contacts.sh <target> <system_dir> [more dirs ...]
#
# <target> is one of: BCHE PPARG GR VDR HNF4A   (key residues from targets.conf)
#
# Produces, in each <system_dir>:
#   index_hb_essential.ndx     key-residue group, used by 03 for hb_essential.xvg
#   index_key.ndx              one group per key residue
#   mindist_<RES>.xvg          minimum ligand-residue heavy-atom distance
#   contacts_<RES>.xvg         number of contacts within 0.45 nm
# and for BChE additionally:
#   dist_LIG_CAS_com.xvg       ligand COM to catalytic anionic site COM
#   dist_LIG_PAS_com.xvg       ligand COM to peripheral anionic site COM
#   dist_TRP82_UNLring_com.xvg ring-centroid distance for the pi-stacking analysis
#   mindist_TRP82_UNLring.xvg
#
# A residue is counted as "in contact" when its minimum heavy-atom distance to
# the ligand falls below 0.35 nm; the 0.45 nm contact count is reported for
# reference. The occupancy figure uses the 0.35 nm criterion.
# ---------------------------------------------------------------------------
set -uo pipefail
GMX="${GMX:-gmx}"
HERE="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=targets.conf
source "$HERE/targets.conf"

TARGET="${1:?usage: bash 04_site_contacts.sh <BCHE|PPARG|GR|VDR|HNF4A> <system_dir> ...}"; shift
eval "KEYRES=\${${TARGET}_KEYRES:-}"
[[ -n "$KEYRES" ]] || { echo "unknown target: $TARGET" >&2; exit 1; }

# three-letter residue names, read once per system from the structure
resname_of() {           # $1 = reference pdb, $2 = resid
    awk -v r="$2" '/^ATOM/ && substr($0,23,4)+0==r {print substr($0,18,3); exit}' "$1"
}

for SYS in "$@"; do
    [[ -d "$SYS" ]] || { echo "skip: $SYS" >&2; continue; }
    echo "=== $SYS  ($TARGET)"
    cd "$SYS" || continue
    TPR=MD.tpr; XTC=MD_center_300.xtc
    LIG="${LIGRES:-UNL}"

    "$GMX" editconf -f "$TPR" -o _ref.pdb >/dev/null 2>&1

    # --- key residues as one group, for the hb_essential analysis ----------
    SEL="$(tr ' ' '\n' <<< "$KEYRES" | sed 's/^/r /' | paste -sd'|' -)"
    printf '%s\nname 21 hb_essential\nq\n' "$SEL" \
        | "$GMX" make_ndx -f "$TPR" -o index_hb_essential.ndx

    # --- one group per key residue -----------------------------------------
    { for r in $KEYRES; do printf 'r %s\n' "$r"; done; printf 'q\n'; } \
        | "$GMX" make_ndx -f "$TPR" -o index_key.ndx

    g=21
    for r in $KEYRES; do
        RN="$(resname_of _ref.pdb "$r")"
        printf '%d\n%s\n' "$g" "$LIG" \
            | "$GMX" mindist -s "$TPR" -f "$XTC" -n index_key.ndx -d 0.45 -group \
                     -od "mindist_${RN}${r}.xvg" -on "contacts_${RN}${r}.xvg"
        g=$((g + 1))
    done

    # --- BChE: site-level distances ----------------------------------------
    if [[ "$TARGET" == "BCHE" ]]; then
        cas="$(tr ' ' '\n' <<< "$BCHE_CAS" | sed 's/^/resid /' | paste -sd' or ' -)"
        pas="$(tr ' ' '\n' <<< "$BCHE_PAS" | sed 's/^/resid /' | paste -sd' or ' -)"
        "$GMX" distance -s "$TPR" -f "$XTC" -oall dist_LIG_CAS_com.xvg \
            -select "com of resname $LIG plus com of (group \"Protein\" and ($cas))"
        "$GMX" distance -s "$TPR" -f "$XTC" -oall dist_LIG_PAS_com.xvg \
            -select "com of resname $LIG plus com of (group \"Protein\" and ($pas))"

        # Trp82 indole ring vs the ligand aromatic ring
        printf 'r 82 & a %s\nname 21 TRP82_RING\nr %s & a C6 C7 C8 C9 C10 C11\nname 22 UNL_RING\nq\n' \
               "$BCHE_RING_ATOMS" "$LIG" | "$GMX" make_ndx -f "$TPR" -o index_pi.ndx
        "$GMX" distance -s "$TPR" -f "$XTC" -n index_pi.ndx -oall dist_TRP82_UNLring_com.xvg \
            -select 'com of group "TRP82_RING" plus com of group "UNL_RING"'
        printf 'TRP82_RING\nUNL_RING\n' \
            | "$GMX" mindist -s "$TPR" -f "$XTC" -n index_pi.ndx -od mindist_TRP82_UNLring.xvg
    fi

    rm -f _ref.pdb
    cd - >/dev/null
done
echo "04 done"
