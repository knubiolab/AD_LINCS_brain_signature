# AD LINCS brain perturbation signature analysis

This repository contains the code used to select brain-relevant LINCS/CMAP perturbation signatures, build drug-level perturbation signatures, compute Alzheimer's disease connectivity scores, and Figure 3E.

The repository is organized so that raw public data and large intermediate objects are kept outside version control. Scripts use command-line arguments and do not depend on user-specific paths.

## Repository structure

```text
AD_LINCS_brain_signature/
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
├── environment.yml
├── data/
│   └── README.md
├── scripts/
│   ├── 01_select_brain_specific_signatures.R
│   ├── 02_build_perturbation_signatures_by_drug.R
│   ├── 03_merge_perturbation_signatures.R
│   ├── 04_compute_connectivity_scores.R
|   ├── 05_cluster_tier2_compounds.R
│   └── 06_plot_cluster6_results.R
├── resources/
│   ├── deg_manifest.tsv
│   ├── PPARG_PPARgene_targets.txt
│   ├── cluster6_drugs.txt
│   └── README.md
└── example/
      └── run_example.sh
```

## Workflow

1. Select brain-relevant perturbation signatures from LINCS Phase I, LINCS Phase II, and CMAP 2020 metadata.
2. Build one perturbation signature per drug to avoid loading the full LINCS/CMAP matrices into memory.
3. Merge drug-level perturbation signatures into one matrix.
4. Compute connectivity scores for Alzheimer's disease DEG signatures across five predefined gene axes:
   - `non_PGBV`: genes outside PPARG, NR3C1, BCHE, and VEGFR-related gene sets
   - `B`: BCHE-axis genes
   - `PG`: PPARG + NR3C1 target genes
   - `PGB`: PPARG + NR3C1 + BCHE-axis genes
   - `PGBV`: PPARG + NR3C1 + BCHE-axis + VEGFR-axis genes
5. Plot cluster 6 rank and drug-specific connectivity figures.

## Installation

Create the conda environment:

```bash
conda env create -f environment.yml
conda activate ad_lincs_signature
```

Some Bioconductor packages may need to be installed inside R if they are unavailable from conda on your platform:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("PharmacoGx", "cmapR", "KEGGREST", "Biobase", "SummarizedExperiment"))
```

## Input data

Large raw files are not distributed in this repository. Download the public LINCS/CMAP files and place them under `data/raw/`.

Required small resource files are described in `resources/README.md`.

## Example run

```bash
bash example/run_example.sh "$(pwd)" "$(pwd)/results"
```

The example script assumes the following:

- Raw LINCS/CMAP files are under `data/raw/`.
- DEG files are under `resources/Prefrontal_cortex_DEG/`.
- `resources/Cmap_genes.txt` is present.
- `resources/cluster6_drugs.txt` contains the final cluster 6 drug list.

## Reproducibility notes

The original exploratory scripts used local paths and repeated scoring blocks for multiple exploratory datasets. This public version keeps only the analysis steps required for the manuscript figures and removes user-specific paths, commented exploratory analyses, and duplicated code.

Drug-level perturbation signatures are generated one drug at a time for memory efficiency. Running all drugs together would not collapse them into a single drug signature when using the higher-level PharmacoGx workflow, but it would require more memory and would make failed jobs harder to restart.

The default batch policy is `pooled`, matching the original analysis in which source-specific batch identifiers were collapsed before perturbation signature estimation. Use `--batch-policy source` for sensitivity checks.
