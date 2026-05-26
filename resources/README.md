# Resources

Required files:

- `deg_manifest.tsv`: Manifest for Alzheimer's disease DEG signatures.
- `Prefrontal_cortex_DEG/*_up.txt` and `Prefrontal_cortex_DEG/*_down.txt`: DEG gene lists. Each file should contain one gene symbol per line.
- `Cmap_genes.txt`: CMAP gene annotation table with at least `Entrez.ID`, `Symbol`, and `Type` columns.
- `PPARG_PPARgene_targets.txt`: PPARG target symbols collected from PPARgene.
- `cluster6_drugs.txt`: Final cluster 6 drug names or perturbagen IDs, one per line.

Optional target gene files:

If available, place Enrichr-style JSON files under `resources/target_gene_sets/`:

```text
resources/target_gene_sets/
├── PPARG_ChEA_Transcription_Factor_Targets.txt
├── PPARG_JASPAR_Predicted_Human_Transcription_Factor_Targets_2025.txt
├── NR3C1_ChEA_Transcription_Factor_Targets.txt
├── NR3C1_ENCODE_Transcription_Factor_Targets.txt
└── NR3C1_JASPAR_Predicted_Human_Transcription_Factor_Targets_2025.txt
```

If local BCHE or VEGFR target files are available, they can be supplied as:

```text
resources/BCHE_targets.txt
resources/VEGFR_targets.txt
```

Otherwise, the connectivity script attempts to retrieve BCHE and VEGFR pathway genes from KEGG using `KEGGREST`.
