Target-axis gene set files
==========================

This directory should contain the target-axis input gene sets used for axis-restricted connectivity:

  PPARG_PPARgene_targets.txt
  PPARG_ChEA_targets.txt
  PPARG_JASPAR_2025_targets.txt
  NR3C1_ChEA_targets.txt
  NR3C1_ENCODE_targets.txt
  NR3C1_JASPAR_2025_targets.txt
  BCHE_KEGG_hsa04725_cholinergic_synapse.txt
  VEGFR_KEGG_hsa04370_VEGF_signaling.txt

PPARG/NR3C1 source files are copied automatically only when plausible local files are found
under --target-source-dir or --source-root. If these files are not found, place the curated
gene-list files here manually and commit them with resources/target_gene_sets/target_gene_set_manifest.tsv.

BChE/VEGFR KEGG pathway files are regenerated when KEGGREST is installed and internet access is available.

Recommended format: one gene symbol per line, no header.
