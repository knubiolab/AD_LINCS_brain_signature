#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--lincs-dir", type = "character", default = file.path("data", "raw"),
              help = "Directory containing GSE92742, GSE70138, and CMAP_2020 subdirectories."),
  make_option("--out-dir", type = "character", default = file.path("results", "signature_selection"),
              help = "Output directory.")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  if (!file.exists(path)) {
    stop("Missing input file: ", path)
  }
  data.table::fread(path, sep = "\t", header = TRUE, data.table = FALSE,
                    quote = "", fill = TRUE, encoding = "UTF-8")
}

drop_pert_types <- function(x, patterns) {
  if (!"pert_type" %in% colnames(x)) {
    stop("Input metadata does not contain a pert_type column.")
  }
  keep <- !grepl(paste(patterns, collapse = "|"), x$pert_type)
  x[keep, , drop = FALSE]
}

parse_dose_um <- function(x) {
  z <- tolower(trimws(as.character(x)))
  z[z %in% c("", "na", "nan", "-666")] <- NA_character_
  numeric_part <- suppressWarnings(as.numeric(sub("^\\s*([-+0-9.eE]+).*$", "\\1", z)))
  unit <- rep(NA_real_, length(z))
  unit[grepl("nm", z)] <- 1e-3
  unit[grepl("um|µm|μm", z)] <- 1
  unit[grepl("mm", z)] <- 1e3
  unit[is.na(unit) & !is.na(numeric_part)] <- 1
  numeric_part * unit
}

parse_duration_h <- function(x) {
  z <- tolower(trimws(as.character(x)))
  z[z %in% c("", "na", "nan", "-666")] <- NA_character_
  numeric_part <- suppressWarnings(as.numeric(sub("^\\s*([-+0-9.eE]+).*$", "\\1", z)))
  out <- numeric_part
  out[grepl("min", z) & !is.na(numeric_part)] <- numeric_part[grepl("min", z) & !is.na(numeric_part)] / 60
  out[grepl("day|d$", z) & !is.na(numeric_part)] <- numeric_part[grepl("day|d$", z) & !is.na(numeric_part)] * 24
  out
}

standardize_signature_table <- function(x, source_name, batchid) {
  rename_map <- c(
    cmap_name = "pert_iname",
    cell_iname = "cell_id"
  )
  for (old in names(rename_map)) {
    if (old %in% colnames(x) && !(rename_map[[old]] %in% colnames(x))) {
      colnames(x)[match(old, colnames(x))] <- rename_map[[old]]
    }
  }

  required <- c("sig_id", "pert_iname", "pert_id", "cell_id", "pert_idose", "pert_itime", "pert_type")
  missing_cols <- setdiff(required, colnames(x))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", source_name, ": ", paste(missing_cols, collapse = ", "))
  }

  dose_um <- parse_dose_um(x$pert_idose)
  duration_h <- parse_duration_h(x$pert_itime)
  xptype <- ifelse(grepl("^ctl|control|vehicle", tolower(x$pert_type)), "control", "perturbation")
  dose_um[xptype == "control" & is.na(dose_um)] <- 0
  dose_um[xptype == "control"] <- 0

  out <- data.frame(
    sig_id = as.character(x$sig_id),
    pert_iname = as.character(x$pert_iname),
    pert_id = as.character(x$pert_id),
    drugid = as.character(x$pert_id),
    cell_id = as.character(x$cell_id),
    cellid = as.character(x$cell_id),
    pert_idose = as.character(x$pert_idose),
    concentration_um = dose_um,
    concentration = dose_um * 1e-6,
    pert_itime = as.character(x$pert_itime),
    duration = duration_h,
    batchid = batchid,
    source = source_name,
    pert_type = as.character(x$pert_type),
    xptype = xptype,
    stringsAsFactors = FALSE
  )

  out <- out[!is.na(out$sig_id) & nzchar(out$sig_id), , drop = FALSE]
  out
}

phase1 <- read_tsv(file.path(opt$lincs_dir, "GSE92742", "GSE92742_Broad_LINCS_sig_info.txt.gz"))
phase1 <- drop_pert_types(phase1, c("vector", "trt_lig", "trt_sh", "trt_oe"))

phase2 <- read_tsv(file.path(opt$lincs_dir, "GSE70138", "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz"))
phase2 <- drop_pert_types(phase2, c("vector", "trt_xpr"))

cmap2020 <- read_tsv(file.path(opt$lincs_dir, "CMAP_2020", "siginfo_beta.txt"))
cmap2020 <- drop_pert_types(cmap2020, c("vector", "cgs", "cns", "ctl_x", "trt_lig",
                                        "trt_sh", "trt_oe", "trt_si", "trt_xpr"))

all_signatures <- rbind(
  standardize_signature_table(phase1, "GSE92742", 1),
  standardize_signature_table(phase2, "GSE70138", 2),
  standardize_signature_table(cmap2020, "CMAP_2020", 3)
)

groups <- list(
  NEU_NPC = c("NEU", "NPC"),
  BT = c("GI1", "LN229", "U251MG", "YH13"),
  NDD = c("MICROGLIA-PSEN1", "SHSY5Y")
)

data.table::fwrite(all_signatures,
                   file.path(opt$out_dir, "all_selected_candidate_signatures.tsv"),
                   sep = "\t")

for (group_name in names(groups)) {
  group_dir <- file.path(opt$out_dir, paste0("perturbations_for_", group_name))
  dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)

  selected <- all_signatures[all_signatures$cellid %in% groups[[group_name]], , drop = FALSE]
  selected <- selected[order(selected$batchid, selected$cellid, selected$pert_id, selected$sig_id), , drop = FALSE]

  data.table::fwrite(selected, file.path(group_dir, "signatures_select.txt"), sep = "\t")
  message(group_name, ": ", nrow(selected), " signatures; ",
          length(unique(selected$pert_id[selected$xptype == "perturbation"])), " perturbagens")
}
