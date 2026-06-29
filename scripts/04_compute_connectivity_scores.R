#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(jsonlite)
  library(PharmacoGx)
})

option_list <- list(
  make_option("--perturbation-rds", type = "character", default = file.path("results", "perturbation_sig_for_NEU_NPC.rds"),
              help = "Merged perturbation signature matrix RDS."),
  make_option("--signatures-select", type = "character", default = file.path("results", "signature_selection", "perturbations_for_NEU_NPC", "signatures_select.txt"),
              help = "signatures_select.txt used to map perturbagen IDs to names."),
  make_option("--resources-dir", type = "character", default = "resources",
              help = "Directory containing deg_manifest.tsv, Cmap_genes.txt, target gene files, and DEG files."),
  make_option("--out-dir", type = "character", default = file.path("results", "connectivity_scores"),
              help = "Output directory."),
  make_option("--nperm", type = "integer", default = 100,
              help = "Number of permutations used by PharmacoGx connectivityScore."),
  make_option("--rank-direction", type = "character", default = "high",
              help = "Use high if larger connectivity means stronger disease similarity; use low otherwise.")
  make_option("--seed", type = "integer", default = 20240601,
            help = "Random seed for fgsea permutation-based connectivity scoring.")
)
opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)

if (!opt$rank_direction %in% c("high", "low")) {
  stop("--rank-direction must be high or low.")
}

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

read_symbols <- function(path) {
  if (!file.exists(path)) {
    warning("Missing gene-symbol file: ", path)
    return(character(0))
  }
  x <- readLines(path, warn = FALSE)
  x <- trimws(x)
  x <- x[!is.na(x) & nzchar(x) & !grepl("^#", x)]
  unique(x)
}

read_enrichr_json_symbols <- function(path) {
  if (!file.exists(path)) {
    return(character(0))
  }
  x <- jsonlite::fromJSON(path)
  symbols <- tryCatch(x$associations$gene$symbol, error = function(e) character(0))
  unique(symbols[!is.na(symbols) & nzchar(symbols)])
}

read_kegg_symbols <- function(pathway_id) {
  if (!requireNamespace("KEGGREST", quietly = TRUE)) {
    warning("KEGGREST is not installed; returning an empty set for ", pathway_id)
    return(character(0))
  }
  out <- tryCatch({
    x <- KEGGREST::keggGet(pathway_id)[[1]]$GENE
    x <- x[grep(";", x)]
    sym <- sub(";.*$", "", x)
    unique(unname(sym[!is.na(sym) & nzchar(sym)]))
  }, error = function(e) {
    warning("Failed to retrieve KEGG pathway ", pathway_id, ": ", conditionMessage(e))
    character(0)
  })
  out
}

load_cmap_genes <- function(path) {
  if (!file.exists(path)) {
    stop("Missing Cmap_genes.txt: ", path)
  }
  x <- data.table::fread(path, sep = "\t", data.table = FALSE)
  required <- c("Entrez.ID", "Symbol", "Type")
  missing_cols <- setdiff(required, colnames(x))
  if (length(missing_cols) > 0) {
    stop("Cmap_genes.txt is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  x <- x[x$Type %in% c("best inferred", "landmark"), c("Entrez.ID", "Symbol")]
  x <- x[!is.na(x$Entrez.ID) & !is.na(x$Symbol), , drop = FALSE]
  x$Entrez.ID <- as.character(x$Entrez.ID)
  x
}

symbols_to_entrez <- function(symbols, cmap_genes) {
  symbols <- unique(trimws(symbols))
  ids <- cmap_genes$Entrez.ID[match(symbols, cmap_genes$Symbol)]
  unique(as.character(ids[!is.na(ids) & nzchar(ids)]))
}

make_deg_signature <- function(up_path, down_path, cmap_genes, down_top_n = NA) {
  up_symbols <- read_symbols(up_path)
  down_symbols <- read_symbols(down_path)

  if (!is.na(down_top_n) && is.finite(as.numeric(down_top_n))) {
    down_symbols <- head(down_symbols, as.integer(down_top_n))
  }

  up_ids <- symbols_to_entrez(up_symbols, cmap_genes)
  down_ids <- symbols_to_entrez(down_symbols, cmap_genes)

  deg <- data.frame(
    feature = c(up_ids, down_ids),
    direction = c(rep(1, length(up_ids)), rep(-1, length(down_ids))),
    stringsAsFactors = FALSE
  )
  deg <- deg[!duplicated(deg$feature), , drop = FALSE]
  rownames(deg) <- deg$feature
  deg[, "direction", drop = FALSE]
}

compute_connectivity <- function(sig_mat, deg, nperm) {
  if (!is.null(seed)) set.seed(seed)
  common <- intersect(rownames(sig_mat), rownames(deg))
  if (length(common) < 10) {
    warning("Fewer than 10 overlapping genes. Returning NA scores.")
    return(data.frame(
      pert_id = colnames(sig_mat),
      connectivity = NA_real_,
      p_value = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  sig_mat <- sig_mat[common, , drop = FALSE]
  deg <- deg[common, , drop = FALSE]

  res <- apply(sig_mat, 2, function(x) {
    ans <- tryCatch(
      PharmacoGx::connectivityScore(x = x, y = deg, method = "fgsea", nperm = nperm),
      error = function(e) c(Connectivity = NA_real_, `P Value` = NA_real_)
    )
    c(connectivity = as.numeric(ans[1]), p_value = as.numeric(ans[2]))
  })
  res <- as.data.frame(t(res))
  res$pert_id <- rownames(res)
  rownames(res) <- NULL
  res[, c("pert_id", "connectivity", "p_value")]
}

perturbation_sig <- readRDS(opt$perturbation_rds)
perturbation_sig <- as.matrix(perturbation_sig)
rownames(perturbation_sig) <- as.character(rownames(perturbation_sig))

signatures_select <- data.table::fread(opt$signatures_select, sep = "\t", data.table = FALSE)
name_map <- unique(signatures_select[, intersect(c("drugid", "pert_iname"), colnames(signatures_select)), drop = FALSE])
if (!all(c("drugid", "pert_iname") %in% colnames(name_map))) {
  name_map <- data.frame(drugid = colnames(perturbation_sig), pert_iname = colnames(perturbation_sig))
}

cmap_genes <- load_cmap_genes(file.path(opt$resources_dir, "Cmap_genes.txt"))

target_dir <- file.path(opt$resources_dir, "target_gene_sets")
pparg_symbols <- unique(c(
  read_symbols(file.path(opt$resources_dir, "PPARG_PPARgene_targets.txt")),
  read_enrichr_json_symbols(file.path(target_dir, "PPARG_ChEA_Transcription_Factor_Targets.txt")),
  read_enrichr_json_symbols(file.path(target_dir, "PPARG_JASPAR_Predicted_Human_Transcription_Factor_Targets_2025.txt"))
))
nr3c1_symbols <- unique(c(
  read_enrichr_json_symbols(file.path(target_dir, "NR3C1_ChEA_Transcription_Factor_Targets.txt")),
  read_enrichr_json_symbols(file.path(target_dir, "NR3C1_ENCODE_Transcription_Factor_Targets.txt")),
  read_enrichr_json_symbols(file.path(target_dir, "NR3C1_JASPAR_Predicted_Human_Transcription_Factor_Targets_2025.txt"))
))

bche_symbols <- read_symbols(file.path(opt$resources_dir, "BCHE_targets.txt"))
if (length(bche_symbols) == 0) {
  bche_symbols <- setdiff(read_kegg_symbols("hsa04725"), "ACHE")
}

vegfr_symbols <- read_symbols(file.path(opt$resources_dir, "VEGFR_targets.txt"))
if (length(vegfr_symbols) == 0) {
  vegfr_symbols <- setdiff(read_kegg_symbols("hsa04370"), c("KDR", "VEGFA"))
}

pparg_ids <- symbols_to_entrez(pparg_symbols, cmap_genes)
nr3c1_ids <- symbols_to_entrez(nr3c1_symbols, cmap_genes)
bche_ids <- symbols_to_entrez(bche_symbols, cmap_genes)
vegfr_ids <- symbols_to_entrez(vegfr_symbols, cmap_genes)

pg_ids <- unique(c(pparg_ids, nr3c1_ids))
pgb_ids <- unique(c(pg_ids, bche_ids))
pgbv_ids <- unique(c(pgb_ids, vegfr_ids))
all_axis_ids <- unique(c(pg_ids, bche_ids, vegfr_ids))

axes <- list(
  non_PGBV = setdiff(rownames(perturbation_sig), all_axis_ids),
  B = bche_ids,
  PG = pg_ids,
  PGB = pgb_ids,
  PGBV = pgbv_ids
)

axis_summary <- data.frame(
  axis = names(axes),
  n_genes_total = vapply(axes, length, integer(1)),
  n_genes_in_signature = vapply(axes, function(g) length(intersect(g, rownames(perturbation_sig))), integer(1)),
  stringsAsFactors = FALSE
)
data.table::fwrite(axis_summary, file.path(opt$out_dir, "axis_gene_counts.tsv"), sep = "\t")

manifest_path <- file.path(opt$resources_dir, "deg_manifest.tsv")
if (!file.exists(manifest_path)) {
  stop("Missing DEG manifest: ", manifest_path)
}
manifest <- data.table::fread(manifest_path, sep = "\t", data.table = FALSE)
manifest <- manifest[toupper(as.character(manifest$enabled)) %in% c("TRUE", "T", "1", "YES"), , drop = FALSE]

all_scores <- list()

for (i in seq_len(nrow(manifest))) {
  signature_id <- manifest$signature_id[i]
  up_path <- file.path(opt$resources_dir, manifest$up_file[i])
  down_path <- file.path(opt$resources_dir, manifest$down_file[i])
  down_top_n <- suppressWarnings(as.numeric(manifest$down_top_n[i]))

  deg <- make_deg_signature(up_path, down_path, cmap_genes, down_top_n = down_top_n)

  for (axis_name in names(axes)) {
    genes <- intersect(axes[[axis_name]], rownames(perturbation_sig))
    axis_mat <- perturbation_sig[genes, , drop = FALSE]
    axis_seed <- opt$seed + i * 100L + match(axis_name, names(axes))
    score <- compute_connectivity(axis_mat, deg, nperm = opt$nperm, seed = axis_seed)
    score$signature_id <- signature_id
    score$axis <- axis_name
    score$n_axis_genes <- nrow(axis_mat)
    score$n_deg_genes <- nrow(deg)

    all_scores[[paste(signature_id, axis_name, sep = "__")]] <- score
  }
}

scores <- data.table::rbindlist(all_scores, use.names = TRUE, fill = TRUE)
scores$pert_name <- name_map$pert_iname[match(scores$pert_id, name_map$drugid)]
scores$pert_name[is.na(scores$pert_name)] <- scores$pert_id

scores <- as.data.frame(scores)
scores$rank <- NA_real_

for (key in unique(paste(scores$signature_id, scores$axis, sep = "__"))) {
  idx <- paste(scores$signature_id, scores$axis, sep = "__") == key
  vals <- scores$connectivity[idx]
  scores$rank[idx] <- if (opt$rank_direction == "high") {
    rank(-vals, ties.method = "average", na.last = "keep")
  } else {
    rank(vals, ties.method = "average", na.last = "keep")
  }
}

scores <- scores[, c("signature_id", "axis", "pert_id", "pert_name", "connectivity",
                     "p_value", "rank", "n_axis_genes", "n_deg_genes")]

data.table::fwrite(scores, file.path(opt$out_dir, "connectivity_scores_all.tsv"), sep = "\t")
message("Saved: ", file.path(opt$out_dir, "connectivity_scores_all.tsv"))
