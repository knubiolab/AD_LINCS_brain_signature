# Resource files

This directory contains small input files required to reproduce the AD connectivity and target-axis analyses.

## Files

- `Cmap_genes.txt`: original LINCS/CMap gene annotation table.
- `Cmap_genes_L1000_BING.tsv`: genes retained for landmark and best inferred gene-space matching.
- `AD_DEG/`: standardized copies of AD up/down DEG query lists.
- `deg_manifest.tsv`: metadata and standardized paths for each AD DEG query signature.
- `AD_input_DEG_up_down_count_summary.tsv`: raw, CMap/BING-matched, and actually used up/down DEG counts.
- `target_gene_sets/`: target-axis input gene sets for PPARG, NR3C1/GR, BChE, and VEGFR.
- `resource_manifest.tsv`: source paths and generation notes for copied/generated files.

The files were collected using `scripts/00_collect_input_files_for_repo.R`.
