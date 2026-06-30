# AD LINCS brain perturbation signature analysis

This repository contains code and small input resources used for the manuscript:

**Brain cell-centric transcriptomic screening for Alzheimer’s disease prioritizes motesanib and identifies in vitro butyrylcholinesterase inhibition**

The repository supports reproduction of the main computational workflow, including selection of brain-relevant LINCS/CMap perturbation signatures, construction of drug-level perturbation signatures, Alzheimer’s disease connectivity scoring, Tier 1 candidate ranking, Tier 2 BR-focused clustering, and generation of quantitative source data and plotting outputs for the main manuscript figures.

Large public LINCS/CMap matrices and large intermediate objects are not tracked in git. Small input files required for the connectivity analysis, including AD DEG lists, CMap gene-space annotation, and target-axis gene sets, are organized under `resources/`.

## Repository structure

```text
AD_LINCS_brain_signature/
├── README.md
├── LICENSE
├── environment.yml
├── sessionInfo.txt
├── data/
│   └── README.md
├── resources/
│   ├── README.md
│   ├── Cmap_genes.txt
│   ├── Cmap_genes_L1000_BING.tsv
│   ├── deg_manifest.tsv
│   ├── AD_input_DEG_up_down_count_summary.tsv
│   ├── AD_DEG/
│   │   ├── TL/
│   │   ├── FL/
│   │   └── PL/
│   ├── target_gene_sets/
│   │   ├── README_target_gene_sets.txt
│   │   └── target_gene_set_manifest.tsv
│   └── resource_manifest.tsv
├── scripts/
│   ├── 01_select_brain_specific_signatures.R
│   ├── 02_build_perturbation_signatures_by_drug.R
│   ├── 03_merge_perturbation_signatures.R
│   ├── 04_compute_connectivity_scores.R
│   ├── 05_cluster_tier2_compounds.R
│   └── 06_plot_cluster6_results.R
└── example/
    └── run_example.sh
```

## Installation

Create the conda environment:

```bash
conda env create -f environment.yml
conda activate ad_lincs_signature
```

The `environment.yml` file pins the main R and Bioconductor package versions used in the manuscript. The exact R session used for the analyses is also documented in `sessionInfo.txt`.

## Input data

Large LINCS/CMap files are not included in this repository. Place the public LINCS Phase I, LINCS Phase II, and CMap 2020 files under `data/raw/` as described in `data/README.md`.

Expected layout:

```text
data/raw/
├── GSE92742/
│   ├── GSE92742_Broad_LINCS_sig_info.txt.gz
│   ├── GSE92742_Broad_LINCS_cell_info.txt.gz
│   └── annotated_GSE92742_Broad_LINCS_Level5_COMPZ_n473647x12328.gctx
├── GSE70138/
│   ├── GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz
│   └── GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx
└── CMAP_2020/
    ├── siginfo_beta.txt
    ├── compoundinfo_beta.txt
    ├── level5_beta_ctl_n58022x12328.gctx
    └── level5_beta_trt_cp_n720216x12328.gctx
```

Small input files required for connectivity scoring are stored under `resources/`. These include:

* `Cmap_genes.txt`: LINCS/CMap gene annotation table.
* `Cmap_genes_L1000_BING.tsv`: landmark and best-inferred gene subset used for query-gene matching.
* `deg_manifest.tsv`: metadata for AD up/down DEG input signatures.
* `AD_DEG/`: standardized AD DEG query lists.
* `AD_input_DEG_up_down_count_summary.tsv`: raw, CMap/BING-matched, and actually used up/down DEG counts for each AD query signature.
* `target_gene_sets/`: PPARG, NR3C1/GR, BChE, and VEGFR target-axis gene sets.
* `resource_manifest.tsv`: source and generation notes for resource files.

## Analysis workflow

### 1. Select brain-relevant perturbation signatures

```bash
Rscript scripts/01_select_brain_specific_signatures.R \
  --lincs-dir data/raw \
  --out-dir results/signature_selection
```

This step selects perturbation signatures from LINCS Phase I, LINCS Phase II, and CMap 2020 for the brain-relevant reference groups used in the manuscript.

### 2. Build drug-level perturbation signatures

Drug-level signatures are generated one compound at a time for memory efficiency.

```bash
Rscript scripts/02_build_perturbation_signatures_by_drug.R \
  --lincs-dir data/raw \
  --group-dir results/signature_selection/perturbations_for_NEU_NPC \
  --out-dir results/signature_selection/perturbations_for_NEU_NPC \
  --drug-index 1 \
  --batch-policy pooled \
  --nthread 1
```

### 3. Merge drug-level signatures

```bash
Rscript scripts/03_merge_perturbation_signatures.R \
  --input-dir results/signature_selection/perturbations_for_NEU_NPC \
  --output-rds results/perturbation_sig_for_NEU_NPC.rds \
  --statistic estimate
```

### 4. Compute AD connectivity scores

```bash
Rscript scripts/04_compute_connectivity_scores.R \
  --perturbation-rds results/perturbation_sig_for_NEU_NPC.rds \
  --signatures-select results/signature_selection/perturbations_for_NEU_NPC/signatures_select.txt \
  --resources-dir resources \
  --out-dir results/connectivity_scores_NEU_NPC \
  --nperm 1000 \
  --rank-direction high \
  --seed 20240601
```

Connectivity is computed for AD DEG query signatures across five predefined target-axis gene spaces:

* `non_PGBV`: genes outside PPARG, NR3C1, BChE, and VEGFR-related gene sets.
* `B`: BChE-axis genes.
* `PG`: PPARG + NR3C1/GR target genes.
* `PGB`: PPARG + NR3C1/GR + BChE-axis genes.
* `PGBV`: PPARG + NR3C1/GR + BChE-axis + VEGFR-axis genes.

The random seed used for permutation-based scoring is recorded in the output parameter file.

### 5. Reproduce Tier 1 candidate ranking

```bash
Rscript scripts/07_compute_tier1_rankings.R \
  --scores-dir results \
  --resources-dir resources \
  --out-dir results/tier1_rankings
```

Tier 1 ranking summarizes disease similarity and positive-control drug similarity across reconstructed references. Self-connectivity of positive-control drugs is excluded before calculating mean drug similarity.

### 6. Reproduce Tier 2 BR-focused clustering

```bash
Rscript scripts/05_cluster_tier2_compounds.R \
  --signature-rds results/perturbation_sig_for_NEU_NPC.rds \
  --candidate-rds results/tier1_rankings/drugs_candidate_NEU_NPC.rds \
  --out-dir results/tier2_clustering \
  --seed 20240602 \
  --k 11 \
  --cluster-id 6
```

This script performs k-means clustering of BR/NEU-NPC candidate perturbation signatures using a fixed seed and exports the full cluster assignment table, silhouette diagnostics, clustering parameters, and the regenerated Cluster 6 drug list.

## Example run

A minimal example workflow can be run with:

```bash
bash example/run_example.sh "$(pwd)" "$(pwd)/results"
```

The example assumes that public LINCS/CMap files are placed under `data/raw/` and that required small resource files are available under `resources/`.

## Reproducibility notes

* Main package versions are pinned in `environment.yml`.
* The exact R session used for the manuscript analyses is documented in `sessionInfo.txt`.
* Stochastic steps use explicit random seeds:

  * Connectivity/permutation scoring: `20240601`
  * Tier 2 k-means clustering: `20240602`
* Clustering outputs include full cluster assignments, the selected Cluster 6 members, silhouette diagnostics, and clustering parameters.
* Resource files include the actual AD DEG inputs, CMap gene-space annotation, target-axis gene sets, and manifests documenting file provenance.
* Large public LINCS/CMap matrices and large intermediate RDS objects are excluded from git because of file size.
* Quantitative source-data tables are provided where possible to support figure reproduction.

## Output files

Common output directories include:

```text
results/
├── signature_selection/
├── perturbation_sig_for_NEU_NPC.rds
├── connectivity_scores_NEU_NPC/
├── tier1_rankings/
└── tier2_clustering/
```

Key reproducibility outputs include:

```text
results/connectivity_scores_NEU_NPC/connectivity_scores_all.tsv
results/connectivity_scores_NEU_NPC/connectivity_scoring_parameters.tsv
results/tier1_rankings/tier1_rankings_all.tsv
results/tier1_rankings/tier1_candidates.tsv
results/tier2_clustering/kmeans_cluster_assignments.tsv
results/tier2_clustering/cluster6_drugs.txt
results/tier2_clustering/kmeans_clustering_parameters.tsv
```

## License

This repository is distributed under the MIT license.
