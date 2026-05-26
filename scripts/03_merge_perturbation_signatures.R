#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--input-dir", type = "character", default = file.path("results", "signature_selection", "perturbations_for_NEU_NPC"),
              help = "Directory containing drug-level perturbation_*.rds files."),
  make_option("--output-rds", type = "character", default = file.path("results", "perturbation_sig_for_NEU_NPC.rds"),
              help = "Output RDS file for the merged gene-by-drug matrix."),
  make_option("--statistic", type = "character", default = "estimate",
              help = "Statistic to extract from each drug-level signature, e.g. estimate or tstat.")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(dirname(opt$output_rds), recursive = TRUE, showWarnings = FALSE)

files <- list.files(opt$input_dir, pattern = "^perturbation_.*\\.rds$", full.names = TRUE)
files <- files[!grepl("\\.lock$|fake", files)]
if (length(files) == 0) {
  stop("No perturbation_*.rds files found in: ", opt$input_dir)
}

extract_signature <- function(obj, statistic) {
  if (is.list(obj) && "result" %in% names(obj)) {
    res <- obj$result
    if (!statistic %in% colnames(res)) {
      stop("Statistic ", statistic, " not found in result for drug ", obj$drug_id)
    }
    v <- res[, statistic]
    names(v) <- rownames(res)
    drug_id <- obj$drug_id
    return(list(drug_id = drug_id, values = v))
  }

  if (is.array(obj) && length(dim(obj)) == 3) {
    stat_names <- dimnames(obj)[[3]]
    if (!statistic %in% stat_names) {
      stop("Statistic ", statistic, " not found in PharmacoSig array.")
    }
    mat <- obj[, , statistic, drop = TRUE]
    if (is.null(dim(mat))) {
      drug_id <- dimnames(obj)[[2]][1]
      names(mat) <- dimnames(obj)[[1]]
      return(list(drug_id = drug_id, values = mat))
    }
    out <- lapply(seq_len(ncol(mat)), function(j) {
      v <- mat[, j]
      names(v) <- rownames(mat)
      list(drug_id = colnames(mat)[j], values = v)
    })
    return(out)
  }

  stop("Unsupported RDS object type: ", paste(class(obj), collapse = ", "))
}

items <- list()
for (path in files) {
  obj <- readRDS(path)
  item <- extract_signature(obj, opt$statistic)
  if (is.list(item) && length(item) > 0 && is.list(item[[1]]) && "values" %in% names(item[[1]])) {
    items <- c(items, item)
  } else {
    items <- c(items, list(item))
  }
}

all_genes <- Reduce(union, lapply(items, function(x) names(x$values)))
merged <- matrix(NA_real_, nrow = length(all_genes), ncol = length(items),
                 dimnames = list(all_genes, vapply(items, function(x) x$drug_id, character(1))))

for (i in seq_along(items)) {
  v <- items[[i]]$values
  merged[names(v), i] <- v
}

if (any(duplicated(colnames(merged)))) {
  warning("Duplicated drug IDs were found. Keeping the first occurrence for each ID.")
  merged <- merged[, !duplicated(colnames(merged)), drop = FALSE]
}

saveRDS(merged, opt$output_rds)

summary_path <- sub("\\.rds$", "_summary.tsv", opt$output_rds)
summary <- data.frame(
  statistic = opt$statistic,
  n_genes = nrow(merged),
  n_drugs = ncol(merged),
  output_rds = opt$output_rds,
  stringsAsFactors = FALSE
)
data.table::fwrite(summary, summary_path, sep = "\t")

message("Saved merged signature matrix: ", opt$output_rds)
message("Dimensions: ", nrow(merged), " genes x ", ncol(merged), " drugs")
