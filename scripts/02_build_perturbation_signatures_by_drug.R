#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(cmapR)
  library(PharmacoGx)
})

option_list <- list(
  make_option("--lincs-dir", type = "character", default = file.path("data", "raw"),
              help = "Directory containing GSE92742, GSE70138, and CMAP_2020 subdirectories."),
  make_option("--group-dir", type = "character", default = file.path("results", "signature_selection", "perturbations_for_NEU_NPC"),
              help = "Directory containing signatures_select.txt."),
  make_option("--drug-index", type = "integer", default = NA,
              help = "One-based index of the perturbagen to process."),
  make_option("--drug-id", type = "character", default = NA,
              help = "Perturbagen ID to process. Overrides --drug-index if provided."),
  make_option("--out-dir", type = "character", default = NA,
              help = "Output directory. Defaults to --group-dir."),
  make_option("--batch-policy", type = "character", default = "pooled",
              help = "Batch policy: pooled or source."),
  make_option("--nthread", type = "integer", default = 1,
              help = "Number of threads passed to PharmacoGx."),
  make_option("--overwrite", action = "store_true", default = FALSE,
              help = "Overwrite existing perturbation RDS file.")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.na(opt$out_dir)) {
  opt$out_dir <- opt$group_dir
}
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

if (!opt$batch_policy %in% c("pooled", "source")) {
  stop("--batch-policy must be either pooled or source.")
}

signatures_path <- file.path(opt$group_dir, "signatures_select.txt")
if (!file.exists(signatures_path)) {
  stop("Missing signatures_select.txt: ", signatures_path)
}

sig <- data.table::fread(signatures_path, sep = "\t", data.table = FALSE)
required <- c("sig_id", "drugid", "cellid", "xptype", "concentration", "duration", "batchid", "source")
missing_cols <- setdiff(required, colnames(sig))
if (length(missing_cols) > 0) {
  stop("signatures_select.txt is missing columns: ", paste(missing_cols, collapse = ", "))
}

drug_ids <- sort(unique(sig$drugid[sig$xptype == "perturbation" & !is.na(sig$drugid) & nzchar(sig$drugid)]))
if (length(drug_ids) == 0) {
  stop("No perturbation drugs found in signatures_select.txt.")
}

if (!is.na(opt$drug_id) && nzchar(opt$drug_id)) {
  target_drug <- opt$drug_id
  if (!target_drug %in% drug_ids) {
    stop("--drug-id not found among selected perturbagens: ", target_drug)
  }
} else {
  if (is.na(opt$drug_index) || opt$drug_index < 1 || opt$drug_index > length(drug_ids)) {
    stop("--drug-index must be between 1 and ", length(drug_ids), ".")
  }
  target_drug <- drug_ids[[opt$drug_index]]
}

safe_name <- gsub("[^A-Za-z0-9_.-]+", "_", target_drug)
out_rds <- file.path(opt$out_dir, paste0("perturbation_", safe_name, ".rds"))
lock_file <- paste0(out_rds, ".lock")

if (file.exists(out_rds) && !opt$overwrite) {
  message("Output already exists, skipping: ", out_rds)
  quit(save = "no", status = 0)
}
if (file.exists(lock_file) && !opt$overwrite) {
  stop("Lock file exists. Another job may be running: ", lock_file)
}
file.create(lock_file)
on.exit(unlink(lock_file), add = TRUE)

read_feature_metadata <- function(lincs_dir) {
  gctx <- file.path(lincs_dir, "GSE92742", "annotated_GSE92742_Broad_LINCS_Level5_COMPZ_n473647x12328.gctx")
  if (!file.exists(gctx)) {
    stop("Missing GSE92742 annotated GCTX file used for feature metadata: ", gctx)
  }
  fdata <- cmapR::read_gctx_meta(gctx, dim = "row")
  rownames(fdata) <- fdata$id
  fdata
}

read_gctx_matrix <- function(path, cids, gene_order) {
  cids <- unique(cids[!is.na(cids) & nzchar(cids)])
  if (length(cids) == 0) {
    return(NULL)
  }
  if (!file.exists(path)) {
    stop("Missing GCTX file: ", path)
  }
  g <- cmapR::parse_gctx(fname = path, cid = cids, matrix_only = TRUE)
  mat <- as.matrix(g@mat)
  mat <- mat[match(gene_order, rownames(mat)), , drop = FALSE]
  missing_genes <- sum(is.na(match(gene_order, rownames(g@mat))))
  if (missing_genes > 0) {
    warning("Missing ", missing_genes, " features in ", basename(path))
  }
  mat
}

target_cells <- unique(sig$cellid[sig$drugid == target_drug & sig$xptype == "perturbation"])
if (length(target_cells) == 0) {
  stop("No target cells found for drug: ", target_drug)
}

sample_rows <- sig[
  sig$cellid %in% target_cells &
    (sig$xptype == "control" | (sig$xptype == "perturbation" & sig$drugid == target_drug)),
  ,
  drop = FALSE
]
sample_rows <- sample_rows[!is.na(sample_rows$sig_id) & nzchar(sample_rows$sig_id), , drop = FALSE]

if (sum(sample_rows$xptype == "perturbation") == 0 || sum(sample_rows$xptype == "control") == 0) {
  stop("Drug ", target_drug, " does not have both perturbation and control samples in the selected cells.")
}

fdata <- read_feature_metadata(opt$lincs_dir)
gene_order <- rownames(fdata)

mat_list <- list()

phase1_rows <- sample_rows[sample_rows$batchid == 1, , drop = FALSE]
if (nrow(phase1_rows) > 0) {
  mat_list[["GSE92742"]] <- read_gctx_matrix(
    file.path(opt$lincs_dir, "GSE92742", "annotated_GSE92742_Broad_LINCS_Level5_COMPZ_n473647x12328.gctx"),
    phase1_rows$sig_id,
    gene_order
  )
}

phase2_rows <- sample_rows[sample_rows$batchid == 2, , drop = FALSE]
if (nrow(phase2_rows) > 0) {
  mat_list[["GSE70138"]] <- read_gctx_matrix(
    file.path(opt$lincs_dir, "GSE70138", "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx"),
    phase2_rows$sig_id,
    gene_order
  )
}

cmap_rows <- sample_rows[sample_rows$batchid == 3, , drop = FALSE]
if (nrow(cmap_rows) > 0) {
  cmap_ctrl <- cmap_rows[cmap_rows$xptype == "control", , drop = FALSE]
  cmap_trt <- cmap_rows[cmap_rows$xptype == "perturbation", , drop = FALSE]

  ctrl_mat <- read_gctx_matrix(
    file.path(opt$lincs_dir, "CMAP_2020", "level5_beta_ctl_n58022x12328.gctx"),
    cmap_ctrl$sig_id,
    gene_order
  )
  trt_mat <- read_gctx_matrix(
    file.path(opt$lincs_dir, "CMAP_2020", "level5_beta_trt_cp_n720216x12328.gctx"),
    cmap_trt$sig_id,
    gene_order
  )

  if (!is.null(ctrl_mat)) mat_list[["CMAP_2020_control"]] <- ctrl_mat
  if (!is.null(trt_mat)) mat_list[["CMAP_2020_treatment"]] <- trt_mat
}

mat_list <- mat_list[!vapply(mat_list, is.null, logical(1))]
if (length(mat_list) == 0) {
  stop("No expression matrices were loaded for drug: ", target_drug)
}

assay_mat <- do.call(cbind, mat_list)
assay_mat <- assay_mat[, !duplicated(colnames(assay_mat)), drop = FALSE]

pdat <- sample_rows[match(colnames(assay_mat), sample_rows$sig_id), , drop = FALSE]
rownames(pdat) <- pdat$sig_id

if (any(is.na(pdat$sig_id))) {
  stop("Failed to match expression columns to sample metadata.")
}

pdat$drugid[pdat$xptype == "control"] <- "control"
pdat$concentration[pdat$xptype == "control" | is.na(pdat$concentration)] <- 0
pdat$duration[is.na(pdat$duration)] <- 0

if (opt$batch_policy == "pooled") {
  pdat$batchid_model <- "pooled"
} else {
  pdat$batchid_model <- paste0("source", pdat$batchid)
}

rank_fun <- getFromNamespace("rankGeneDrugPerturbation", "PharmacoGx")

rank_res <- rank_fun(
  data = t(assay_mat),
  drug = target_drug,
  drug.id = as.character(pdat$drugid),
  drug.concentration = as.numeric(pdat$concentration),
  type = as.character(pdat$cellid),
  xp = as.character(pdat$xptype),
  batch = as.character(pdat$batchid_model),
  duration = as.character(pdat$duration),
  single.type = FALSE,
  nthread = opt$nthread,
  verbose = FALSE
)

if (is.null(rank_res$all)) {
  stop("PharmacoGx returned no all-type perturbation signature for drug: ", target_drug)
}

result <- rank_res$all
result <- result[, intersect(c("estimate", "se", "n", "tstat", "fstat", "pvalue", "fdr",
                               "type.fstat", "type.pvalue"), colnames(result)), drop = FALSE]

obj <- list(
  drug_id = target_drug,
  pert_iname = unique(sample_rows$pert_iname[sample_rows$drugid == target_drug])[1],
  result = result,
  sample_info = pdat,
  batch_policy = opt$batch_policy,
  source = "PharmacoGx::rankGeneDrugPerturbation",
  created = as.character(Sys.time()),
  session_info = utils::sessionInfo()
)
class(obj) <- "LINCSDrugPerturbationSignature"

saveRDS(obj, out_rds)
message("Saved: ", out_rds)
