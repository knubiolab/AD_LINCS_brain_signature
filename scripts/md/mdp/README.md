# GROMACS run parameter files

Copy the `.mdp` files that were used for the published simulations into this
directory. The pipeline scripts expect these five names:

| file | stage |
|---|---|
| `ions.mdp` | parameter file for `grompp` before `genion` |
| `EM.mdp` | energy minimisation |
| `NVT.mdp` | 100 ps NVT equilibration |
| `NPT.mdp` | 1 ns NPT equilibration |
| `MD.mdp` | 300 ns production |

For reference, the settings reported in the Methods are:

**Energy minimisation** — steepest descent, terminated when the maximum force
fell below 1000 kJ mol⁻¹ nm⁻¹, 50,000 steps maximum.

**NVT** — 100 ps, velocity-rescale thermostat, 300 K, τ_t = 0.1 ps, position
restraints on heavy atoms.

**NPT** — 1 ns, Berendsen barostat, 1 bar, τ_p = 2.0 ps, position restraints
retained.

**Production** — 300 ns at 300 K and 1 bar; V-rescale thermostat (τ_t = 0.1 ps)
with Protein_LIG and solvent coupled separately; Parrinello-Rahman barostat
(τ_p = 2.0 ps); 2 fs time step; all bonds constrained with LINCS; Particle Mesh
Ewald electrostatics with a 1.2 nm real-space cutoff and 0.16 nm Fourier grid
spacing; force-switched van der Waals between 1.0 and 1.2 nm; coordinates
written every 10 ps.

Replicates differ only in `gen_seed` of `NVT.mdp`; `02_equilibrate_and_produce.sh`
substitutes that value when a seed is passed on the command line.
