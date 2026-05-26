#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="${1:-$(pwd)}"
OUT_DIR="${2:-${PROJECT_DIR}/results}"
DATA_DIR="${DATA_DIR:-${PROJECT_DIR}/data/raw}"
RESOURCES_DIR="${RESOURCES_DIR:-${PROJECT_DIR}/resources}"
GROUP="${GROUP:-NEU_NPC}"
BATCH_POLICY="${BATCH_POLICY:-pooled}"
NTHREAD="${NTHREAD:-1}"
NPERM="${NPERM:-100}"

cd "${PROJECT_DIR}"

mkdir -p "${OUT_DIR}"

Rscript scripts/01_select_brain_specific_signatures.R \
  --lincs-dir "${DATA_DIR}" \
  --out-dir "${OUT_DIR}/signature_selection"

GROUP_DIR="${OUT_DIR}/signature_selection/perturbations_for_${GROUP}"

N_DRUGS="$(awk -F'\t' '
  NR == 1 {
    for (i = 1; i <= NF; i++) {
      if ($i == "drugid") drugid_col = i
      if ($i == "xptype") xptype_col = i
    }
    next
  }
  $xptype_col == "perturbation" {print $drugid_col}
' "${GROUP_DIR}/signatures_select.txt" | sort -u | wc -l | tr -d ' ')"

echo "Processing ${N_DRUGS} drugs for group ${GROUP}"

for IDX in $(seq 1 "${N_DRUGS}"); do
  Rscript scripts/02_build_perturbation_signatures_by_drug.R \
    --lincs-dir "${DATA_DIR}" \
    --group-dir "${GROUP_DIR}" \
    --out-dir "${GROUP_DIR}" \
    --drug-index "${IDX}" \
    --batch-policy "${BATCH_POLICY}" \
    --nthread "${NTHREAD}"
done

Rscript scripts/03_merge_perturbation_signatures.R \
  --input-dir "${GROUP_DIR}" \
  --output-rds "${OUT_DIR}/perturbation_sig_for_${GROUP}.rds" \
  --statistic estimate

Rscript scripts/04_compute_connectivity_scores.R \
  --perturbation-rds "${OUT_DIR}/perturbation_sig_for_${GROUP}.rds" \
  --signatures-select "${GROUP_DIR}/signatures_select.txt" \
  --resources-dir "${RESOURCES_DIR}" \
  --out-dir "${OUT_DIR}/connectivity_scores" \
  --nperm "${NPERM}" \
  --rank-direction high

Rscript scripts/05_plot_cluster6_results.R \
  --scores "${OUT_DIR}/connectivity_scores/connectivity_scores_all.tsv" \
  --cluster6-drugs "${RESOURCES_DIR}/cluster6_drugs.txt" \
  --out-dir "${OUT_DIR}/figures" \
  --highlight-drugs "donepezil,motesanib" \
  --selected-drug "motesanib"
