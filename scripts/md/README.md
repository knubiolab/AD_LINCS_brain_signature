# Molecular dynamics and binding free energy analyses

Scripts reproducing the molecular dynamics part of the study: 300-ns simulations
of motesanib and the corresponding reference ligands in five targets, the
trajectory analyses behind Figure 4 and Supplementary Figures 6–9, the
conformational clustering behind Supplementary Figure 6a and Supplementary
Table 3, and the MM/GBSA binding free energies.

Simulated systems (each in triplicate from independent velocity seeds):

| target | ligands | receptor |
|---|---|---|
| BChE | motesanib, co-crystallised reference ligand | 4TPK |
| PPARγ | motesanib, LBB78 | 8REJ |
| GR | motesanib, GSK866 | 3E7C |
| VDR *(negative control)* | motesanib, calcitriol, paricalcitol | 7QPP |
| HNF4α *(negative control)* | motesanib | 9IBR |

All equilibrium analyses use the final 100 ns (200–300 ns) of each trajectory,
subsampled at 100 ps.

## Requirements

| tool | version used |
|---|---|
| GROMACS | 2025.4 |
| gmx_MMPBSA | 1.5.0.3 (with AmberTools) |
| PLIP | 3.0.1 |
| R | ≥ 4.2, with `ggplot2`, `dplyr`, `tidyr`, `forcats`, `patchwork` |
| Python | ≥ 3.8 (standard library only) |

Receptor coordinates were parameterised with CHARMM27 and solvated with TIP3P;
ligand topologies were generated with [SwissParam](https://www.swissparam.ch/)
(CHARMM36-compatible) from MOL2 files.

## Layout

```text
scripts/md/
├── README.md
├── targets.conf                       # per-target residue definitions, analysis window
├── mdp/                               # GROMACS .mdp files (see mdp/README.md)
├── 01_prepare_system.sh               # topology, box, solvation, ions, minimisation
├── 02_equilibrate_and_produce.sh      # NVT -> NPT -> 300 ns production
├── 03_postprocess_trajectory.sh       # PBC correction, RMSD/RMSF/Rg/SASA/H-bonds
├── 04_site_contacts.sh                # per-residue contacts, CAS/PAS distances
├── 05_cluster_representative_pose.sh  # GROMOS clustering, representative poses, PLIP
├── 06_mmgbsa.sh                       # MM/GBSA over the equilibrated phase
├── 07_summary_table.sh                # Supplementary Table 3
└── figures/
    ├── graph_BCHE.R                   # Figure 4c–i
    ├── graph_PPARG.R                  # Supplementary Figure 7
    ├── graph_NR3C1.R                  # Supplementary Figure 8
    ├── graph_VDR.R                    # Supplementary Figure 9
    ├── plot_cluster_fraction.R        # Supplementary Figure 6a
    └── plot_stability_hist.R          # comparative RMSD/RMSF histograms
```

## Before running 01

Three steps are interactive and are not scripted:

1. **Split the docked complex.** From the selected docking pose, save the
   receptor as `REC.pdb` (ligand and crystallographic waters removed, hydrogens
   added) and the ligand as `LIG.mol2`.
2. **Parameterise the ligand.** Upload `LIG.mol2` to SwissParam and place the
   returned `unl.itp` / `unl.pdb` (and `unl.par`, if produced) in the system
   directory as `LIG.itp` / `LIG.pdb` / `LIG.prm`.
3. **Edit `topol.top`** after `01` has written it:
   - add `#include "LIG.itp"` below the force-field include;
   - add `LIG   1` to the `[ molecules ]` block, below the protein chain;
   - add the position-restraint guard below the ligand include:
     ```
     #ifdef POSRES
     #include "posre_LIG.itp"
     #endif
     ```
   - confirm that `[ moleculetype ]` in `LIG.itp` names the molecule `LIG`.

## Running the pipeline

```bash
cd scripts/md
SYS=BCHE_motesanib

bash 01_prepare_system.sh           "$SYS"
for seed in 101 202 303; do
  bash 02_equilibrate_and_produce.sh "$SYS" "$seed"
done

# analyses, per replicate directory
for R in "$SYS"/seed_*; do
  bash 04_site_contacts.sh BCHE              "$R"   # builds the key-residue index
  bash 03_postprocess_trajectory.sh          "$R"
  bash 05_cluster_representative_pose.sh     "$R"
  bash 06_mmgbsa.sh                          "$R"
done

bash 07_summary_table.sh "$SYS"/seed_* > stability_summary.tsv
```

`04` is run before `03` because `03` reads the key-residue index file that `04`
writes. The target argument to `04` selects the residue list from
`targets.conf`: `BCHE`, `PPARG`, `GR`, `VDR` or `HNF4A`.

## Figures

The R scripts read the `.xvg` files produced above and expect one subdirectory
per system, named as in `targets.conf`. The four `graph_*.R` scripts set that
location in `WORKDIR` on their first lines; the two `plot_*.R` scripts use the
current directory. Panel assignments and the 200–300 ns analysis window are set
at the top of each script (`t_split`, `t0`).

Each `graph_*.R` script ends at the `pdf()` / `dev.off()` block that writes its
figure, and creates the output directory if it does not exist.

```bash
cd <directory holding the system subdirectories>
Rscript figures/graph_BCHE.R                # set WORKDIR first
Rscript figures/plot_cluster_fraction.R     # also writes cluster_fraction.tsv
Rscript figures/plot_stability_hist.R       # also writes stability_summary.tsv
```

## Notes on specific analyses

**Clustering.** `gmx cluster` is run with `-nofit` on a trajectory that has
already been superposed on the protein backbone by `trjconv -fit rot+trans`, so
it prompts for two groups only: the group the RMSD is computed on (ligand heavy
atoms) and the group written to the output structures (protein + ligand).
Representative poses are extracted for the smallest number of clusters that
together account for ≥50% of the frames in the analysis window.

**Contact occupancy.** A residue counts as being in contact when its minimum
heavy-atom distance to the ligand is below 0.35 nm. `04` also writes the number
of contacts within 0.45 nm for reference; the occupancy panels use the 0.35 nm
criterion.

**Binding sites in BChE.** The catalytic anionic site (CAS) is taken as the
catalytic triad Ser198/His438/Glu325 for the centre-of-mass distance, and the
peripheral anionic site (PAS) as Asp70/Tyr332. Trp82 is reported separately as
the choline-binding aromatic residue.

**MM/GBSA.** Single-trajectory MM/GBSA over 401 evenly spaced frames of the
200–300 ns window, OBC2 generalised Born (`igb = 5`), 0.1 M salt, no entropy
term. Values are therefore relative and are compared only between ligands in the
same receptor.
