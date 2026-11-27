#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 01 - System preparation for a protein-ligand MD simulation
#
#   bash 01_prepare_system.sh <system_dir>
#
# Expects <system_dir> to already contain, from the manual steps described in
# README.md (section "Before running 01"):
#   REC.pdb    receptor, hydrogens added, ligand and crystallographic waters removed
#   LIG.pdb    ligand pose to be simulated (the docked pose, or the co-crystallised
#              pose for a reference ligand)
#   LIG.itp    ligand topology from SwissParam (CHARMM36-compatible)
#   LIG.prm    ligand parameters from SwissParam, if produced
#   ../mdp/*.mdp
#
# Produces: box_sol_ion.gro, topol.top, EM.gro, index.ndx, posre_LIG.itp
#
# Force field / water model: CHARMM27 / TIP3P (selected at the pdb2gmx prompt).
# Box: triclinic, >=1.0 nm solvent padding. Ions: neutralised, 0.1 M NaCl.
# ---------------------------------------------------------------------------
set -euo pipefail
GMX="${GMX:-gmx}"
SYS="${1:?usage: bash 01_prepare_system.sh <system_dir>}"
MDP="$(cd "$(dirname "$0")/mdp" && pwd)"
cd "$SYS"

for f in REC.pdb LIG.pdb LIG.itp; do
    [[ -f "$f" ]] || { echo "missing: $SYS/$f  (see README.md)" >&2; exit 1; }
done

# --- receptor topology -----------------------------------------------------
# Prompts: force field (CHARMM27), then water model (TIP3P).
printf '8\n1\n' | "$GMX" pdb2gmx -f REC.pdb -o REC.gro -p topol.top -ignh

# --- merge the ligand into the complex coordinate file ---------------------
# The ligand block is appended to the receptor coordinates and the atom count
# on line 2 is updated. topol.top must also #include "LIG.itp" after the
# forcefield include, and list "LIG 1" under [ molecules ] -- see README.md.
"$GMX" editconf -f LIG.pdb -o LIG.gro
python3 - <<'PY'
rec = open("REC.gro").read().splitlines()
lig = open("LIG.gro").read().splitlines()
n_rec, n_lig = int(rec[1]), int(lig[1])
body = rec[2:2 + n_rec] + lig[2:2 + n_lig]
with open("conf.gro", "w") as fh:
    fh.write(rec[0] + "\n" + f"{n_rec + n_lig:5d}\n")
    fh.write("\n".join(body) + "\n" + rec[-1] + "\n")
print(f"conf.gro: {n_rec} receptor + {n_lig} ligand atoms")
PY

# --- box, solvent, ions ----------------------------------------------------
"$GMX" editconf -f conf.gro -o box.gro -bt triclinic -d 1.0
"$GMX" solvate  -cp box.gro -cs spc216.gro -p topol.top -o box_sol.gro
"$GMX" grompp   -f "$MDP/ions.mdp" -c box_sol.gro -p topol.top -maxwarn 2 -o ION.tpr
echo SOL | "$GMX" genion -s ION.tpr -p topol.top -conc 0.1 -neutral -o box_sol_ion.gro

# --- energy minimisation ---------------------------------------------------
# Steepest descent to Fmax < 1000 kJ mol^-1 nm^-1 (max 50000 steps; see EM.mdp).
"$GMX" grompp -f "$MDP/EM.mdp" -c box_sol_ion.gro -p topol.top -maxwarn 2 -o EM.tpr
"$GMX" mdrun  -v -deffnm EM

# --- position restraints on ligand heavy atoms -----------------------------
printf '0 & ! a H*\nq\n' | "$GMX" make_ndx -f LIG.gro -o index_LIG.ndx
echo 3 | "$GMX" genrestr -f LIG.gro -n index_LIG.ndx -o posre_LIG.itp -fc 1000 1000 1000
# topol.top must guard this with #ifdef POSRES -- see README.md.

# --- coupling groups: Protein_LIG vs solvent --------------------------------
printf '1 | 13\nq\n' | "$GMX" make_ndx -f EM.gro -o index.ndx

echo "01 done: $SYS"
