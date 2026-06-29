#!/usr/bin/env Rscript

# 05_cluster_tier2_compounds.R
#
# Reproduce Tier 2 BR-focused k-means clustering used to define the 11
# perturbational-signature clusters and the cluster of interest.
#
# This script extracts the clustering-relevant steps from the original
# Fisher_exact_test.R workflow:
#   1) read BR/NEU-NPC candidate compounds,
#   2) subset the BR perturbation signature matrix,
#   3) compute signature-correlation distance,
#   4) run k-means clustering with a fixed random seed,
#   5) export cluster assignments, cluster-of-interest members, parameters,
#      and diagnostic plots/tables.
#
# Example:
# Rscript scripts/05_cluster_tier2_compounds.R \
#   --signature-rds perturbation_sig_for_NEU_NPC_only.rds \
#   --candidate-rds drugs_candidate_NEU_NPC.rds \
#   --out-dir results/tier2_clustering \
#   --seed 20240602 \
#   --k 11 \
#   --cluster-id 6

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(factoextra)
  library(cluster)
})

option_list <- list(
  make_option(
    "--signature-rds",
    type = "character",
    default = "perturbation_sig_for_NEU_NPC_only.rds",
    help = "RDS file containing the BR/NEU-NPC perturbation signature matrix. Rows = genes, columns = pert_id. [default: %default]"
  ),
  make_option(
    "--candidate-rds",
    type = "character",
    default = "drugs_candidate_NEU_NPC.rds",
    help = "RDS file containing candidate pert_ids to cluster. Ignored if --candidate-file is supplied. [default: %default]"
  ),
  make_option(
    "--candidate-file",
    type = "character",
    default = NA,
    help = "Optional text/TSV file containing candidate pert_ids. The script uses a column named pert_id if present, otherwise the first column. [default: %default]"
  ),
  make_option(
    "--out-dir",
    type = "character",
    default = file.path("results", "tier2_clustering"),
    help = "Output directory. [default: %default]"
  ),
  make_option(
    "--seed",
    type = "integer",
    default = 20240602,
    help = "Random seed for k-means clustering and k-selection diagnostics. [default: %default]"
  ),
  make_option(
    "--k",
    type = "integer",
    default = 11,
    help = "Number of k-means clusters. [default: %default]"
  ),
  make_option(
    "--nstart",
    type = "integer",
    default = 100,
    help = "Number of random starts for k-means. [default: %default]"
  ),
  make_option(
    "--cluster-id",
    type = "integer",
    default = 6,
    help = "Cluster ID to export as the cluster of interest. [default: %default]"
  ),
  make_option(
    "--outlier-ids",
    type = "character",
    default = "BRD-A36471396,BRD-K76205745",
    help = "Comma-separated pert_ids flagged as outliers in the original plotting workflow. They are not removed from clustering, but are flagged in outputs. [default: %default]"
  ),
  make_option(
    "--kmax",
    type = "integer",
    default = 20,
    help = "Maximum k for silhouette diagnostic table. Set to 0 to skip. [default: %default]"
  ),
  make_option(
    "--make-plots",
    action = "store_true",
    default = TRUE,
    help = "Write PCA-based cluster diagnostic plots. [default: %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

read_candidate_ids <- function(candidate_rds, candidate_file) {
  if (!is.na(candidate_file) && nzchar(candidate_file)) {
    dat <- data.table::fread(candidate_file)
    if ("pert_id" %in% names(dat)) {
      ids <- dat[["pert_id"]]
    } else {
      ids <- dat[[1]]
    }
  } else {
    ids <- readRDS(candidate_rds)
  }
  ids <- unique(as.character(ids))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  ids
}

clean_outlier_ids <- function(x) {
  ids <- unlist(strsplit(x, ",", fixed = TRUE))
  ids <- trimws(ids)
  ids[nzchar(ids)]
}

message("Reading perturbation signature matrix: ", opt$signature_rds)
sig <- readRDS(opt$signature_rds)
sig <- as.matrix(sig)

message("Reading candidate perturbagens.")
candidate_ids <- read_candidate_ids(opt$candidate_rds, opt$candidate_file)

candidate_ids_found <- intersect(candidate_ids, colnames(sig))
candidate_ids_missing <- setdiff(candidate_ids, colnames(sig))

if (length(candidate_ids_found) < opt$k) {
  stop(
    "Fewer candidate signatures were found than the requested number of clusters. ",
    "Found ", length(candidate_ids_found), " candidate columns, but k = ", opt$k, "."
  )
}

if (length(candidate_ids_missing) > 0) {
  data.table::fwrite(
    data.table(pert_id = candidate_ids_missing),
    file.path(opt$out_dir, "candidate_ids_missing_from_signature_matrix.tsv"),
    sep = "\t"
  )
}

sig <- sig[, candidate_ids_found, drop = FALSE]

# Remove genes with no finite values and compounds with zero variance, if present.
keep_gene <- apply(sig, 1, function(z) any(is.finite(z)))
sig <- sig[keep_gene, , drop = FALSE]

compound_sd <- apply(sig, 2, sd, na.rm = TRUE)
zero_var_compounds <- names(compound_sd)[is.na(compound_sd) | compound_sd == 0]

if (length(zero_var_compounds) > 0) {
  warning("Removing ", length(zero_var_compounds), " zero-variance compound signatures.")
  data.table::fwrite(
    data.table(pert_id = zero_var_compounds),
    file.path(opt$out_dir, "zero_variance_compounds_removed.tsv"),
    sep = "\t"
  )
  sig <- sig[, !colnames(sig) %in% zero_var_compounds, drop = FALSE]
}

message("Final clustering matrix: ", nrow(sig), " genes x ", ncol(sig), " compounds.")

# Original workflow used correlation among perturbation signatures and pearson-based
# clustering diagnostics.
cmap_corr <- cor(sig, use = "pairwise.complete.obs")
cmap_dist <- as.dist(1 - cmap_corr)

# Main k-means clustering.
set.seed(opt$seed)
clusters_cmap <- factoextra::eclust(
  t(sig),
  FUNcluster = "kmeans",
  k = opt$k,
  nstart = opt$nstart,
  hc_metric = "pearson",
  graph = FALSE
)

outlier_ids <- clean_outlier_ids(opt$outlier_ids)

cluster_assignments <- data.table(
  pert_id = names(clusters_cmap$cluster),
  cluster = as.integer(clusters_cmap$cluster),
  flagged_outlier_for_plotting = names(clusters_cmap$cluster) %in% outlier_ids
)

setorder(cluster_assignments, cluster, pert_id)

data.table::fwrite(
  cluster_assignments,
  file.path(opt$out_dir, "kmeans_cluster_assignments.tsv"),
  sep = "\t"
)

clusters_list <- split(cluster_assignments$pert_id, cluster_assignments$cluster)
saveRDS(clusters_list, file.path(opt$out_dir, "kmeans_clusters.rds"))

cluster_of_interest <- cluster_assignments[cluster == opt$cluster_id]
data.table::fwrite(
  cluster_of_interest,
  file.path(opt$out_dir, paste0("cluster", opt$cluster_id, "_pert_ids.tsv")),
  sep = "\t"
)

writeLines(
  cluster_of_interest$pert_id,
  file.path(opt$out_dir, paste0("cluster", opt$cluster_id, "_drugs.txt"))
)

# Silhouette for the selected k.
sil <- cluster::silhouette(clusters_cmap$cluster, cmap_dist)
sil_dt <- data.table(
  pert_id = rownames(sil),
  cluster = as.integer(sil[, "cluster"]),
  neighbor_cluster = as.integer(sil[, "neighbor"]),
  silhouette_width = as.numeric(sil[, "sil_width"])
)

data.table::fwrite(
  sil_dt,
  file.path(opt$out_dir, "kmeans_silhouette_widths.tsv"),
  sep = "\t"
)

# k-selection diagnostics, using the same distance matrix and fixed seeds.
if (!is.null(opt$kmax) && opt$kmax >= 2) {
  k_values <- 2:min(opt$kmax, ncol(sig) - 1)
  if (length(k_values) > 0) {
    k_diag <- rbindlist(lapply(k_values, function(kk) {
      set.seed(opt$seed + kk)
      km_kk <- factoextra::eclust(
        t(sig),
        FUNcluster = "kmeans",
        k = kk,
        nstart = opt$nstart,
        hc_metric = "pearson",
        graph = FALSE
      )
      ss <- cluster::silhouette(km_kk$cluster, cmap_dist)
      sil_per_cluster <- sapply(seq_len(kk), function(j) {
        mean(ss[ss[, "cluster"] == j, "sil_width"], na.rm = TRUE)
      })
      data.table(
        k = kk,
        mean_silhouette_width = mean(ss[, "sil_width"], na.rm = TRUE),
        variance_of_cluster_mean_silhouette = var(sil_per_cluster, na.rm = TRUE)
      )
    }))
    
    data.table::fwrite(
      k_diag,
      file.path(opt$out_dir, "k_selection_silhouette_diagnostics.tsv"),
      sep = "\t"
    )
  }
}

parameters <- data.table(
  parameter = c(
    "signature_rds",
    "candidate_rds",
    "candidate_file",
    "seed",
    "k",
    "nstart",
    "cluster_id",
    "outlier_ids",
    "n_genes_used",
    "n_candidate_ids_input",
    "n_candidate_ids_found",
    "n_candidate_ids_missing",
    "n_compounds_clustered"
  ),
  value = c(
    opt$signature_rds,
    opt$candidate_rds,
    ifelse(is.na(opt$candidate_file), "", opt$candidate_file),
    opt$seed,
    opt$k,
    opt$nstart,
    opt$cluster_id,
    opt$outlier_ids,
    nrow(sig),
    length(candidate_ids),
    length(candidate_ids_found),
    length(candidate_ids_missing),
    ncol(sig)
  )
)

data.table::fwrite(
  parameters,
  file.path(opt$out_dir, "kmeans_clustering_parameters.tsv"),
  sep = "\t"
)

# Diagnostic plots. These do not define clusters; they visualize the fixed
# cluster assignments from the k-means result above.
if (isTRUE(opt$make_plots)) {
  pca <- prcomp(t(sig), center = TRUE, scale. = TRUE)
  plot_dt <- data.table(
    pert_id = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]
  )
  plot_dt <- merge(plot_dt, cluster_assignments, by = "pert_id", all.x = TRUE)
  plot_dt[, cluster := factor(cluster)]
  plot_dt[, highlight_cluster := as.integer(as.character(cluster)) == opt$cluster_id]
  
  p_all <- ggplot(plot_dt, aes(x = PC1, y = PC2, color = cluster)) +
    geom_point(aes(shape = flagged_outlier_for_plotting), size = 2, alpha = 0.9) +
    theme_classic() +
    labs(
      title = "BR perturbation signature-based k-means clusters",
      subtitle = paste0("k = ", opt$k, ", seed = ", opt$seed),
      color = "Cluster",
      shape = "Flagged outlier"
    )
  
  ggsave(
    filename = file.path(opt$out_dir, "kmeans_cluster_pca_all.pdf"),
    plot = p_all,
    width = 5,
    height = 4
  )
  
  p_highlight <- ggplot(plot_dt, aes(x = PC1, y = PC2)) +
    geom_point(data = plot_dt[highlight_cluster == FALSE],
               color = "grey80", size = 1.8, alpha = 0.8) +
    geom_point(data = plot_dt[highlight_cluster == TRUE],
               aes(shape = flagged_outlier_for_plotting),
               size = 2.5, alpha = 0.95) +
    theme_classic() +
    labs(
      title = paste0("Cluster ", opt$cluster_id, " highlighted"),
      subtitle = paste0("k = ", opt$k, ", seed = ", opt$seed),
      shape = "Flagged outlier"
    )
  
  ggsave(
    filename = file.path(opt$out_dir, paste0("cluster", opt$cluster_id, "_highlight_pca.pdf")),
    plot = p_highlight,
    width = 4,
    height = 3.5
  )
}

message("Done. Outputs written to: ", opt$out_dir)
