#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})

option_list <- list(
  make_option("--scores", type = "character", default = file.path("results", "connectivity_scores", "connectivity_scores_all.tsv"),
              help = "Combined connectivity score table from 04_compute_connectivity_scores.R."),
  make_option("--cluster6-drugs", type = "character", default = file.path("resources", "cluster6_drugs.txt"),
              help = "Text file containing cluster 6 drug names or IDs, one per line."),
  make_option("--out-dir", type = "character", default = file.path("results", "figures"),
              help = "Figure output directory."),
  make_option("--highlight-drugs", type = "character", default = "donepezil,motesanib",
              help = "Comma-separated drug names or IDs to overlay on the cluster rank plot."),
  make_option("--selected-drug", type = "character", default = "motesanib",
              help = "Drug shown in the rank/connectivity violin figure.")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

axis_levels <- c("non_PGBV", "B", "PG", "PGB", "PGBV")
axis_labels <- c(non_PGBV = "non-PGBV", B = "B", PG = "PG", PGB = "PGB", PGBV = "PGBV")

scores <- data.table::fread(opt$scores, sep = "\t", data.table = FALSE)
scores <- scores[scores$axis %in% axis_levels, , drop = FALSE]
scores$axis <- factor(scores$axis, levels = axis_levels, labels = unname(axis_labels[axis_levels]))

normalize <- function(x) tolower(trimws(as.character(x)))

read_drug_list <- function(path) {
  if (!file.exists(path)) {
    stop("Missing drug list: ", path)
  }
  x <- readLines(path, warn = FALSE)
  x <- trimws(x)
  x <- x[!is.na(x) & nzchar(x) & !grepl("^#", x)]
  unique(x)
}

matches_drug <- function(df, drug_terms) {
  terms <- normalize(drug_terms)
  normalize(df$pert_id) %in% terms | normalize(df$pert_name) %in% terms
}

cluster6 <- read_drug_list(opt$cluster6_drugs)
highlight <- trimws(unlist(strsplit(opt$highlight_drugs, ",")))
highlight <- highlight[nzchar(highlight)]

scores$display_group <- NA_character_
scores$display_group[matches_drug(scores, cluster6)] <- "cluster 6 drugs"

for (drug in highlight) {
  scores$display_group[matches_drug(scores, drug)] <- drug
}

plot_groups <- c("cluster 6 drugs", highlight)
rank_df <- scores[!is.na(scores$display_group), , drop = FALSE]
rank_df$display_group <- factor(rank_df$display_group, levels = plot_groups)

rank_summary <- aggregate(
  rank ~ display_group + axis,
  data = rank_df,
  FUN = function(x) mean(x, na.rm = TRUE)
)
rank_summary <- rank_summary[!is.na(rank_summary$rank), , drop = FALSE]

rank_plot <- ggplot(rank_summary, aes(x = axis, y = rank, group = display_group,
                                      color = display_group)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  scale_y_reverse() +
  labs(x = NULL, y = "Rank in disease similarity", color = "group") +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

ggsave(file.path(opt$out_dir, "cluster6_rank.pdf"), rank_plot, width = 4.0, height = 3.2)
ggsave(file.path(opt$out_dir, "cluster6_rank.png"), rank_plot, width = 4.0, height = 3.2, dpi = 300)

selected_rows <- scores[matches_drug(scores, opt$selected_drug), , drop = FALSE]
if (nrow(selected_rows) == 0) {
  stop("Selected drug not found in score table: ", opt$selected_drug)
}

selected_rank <- aggregate(rank ~ axis, data = selected_rows, FUN = function(x) mean(x, na.rm = TRUE))
selected_conn <- selected_rows[!is.na(selected_rows$connectivity), , drop = FALSE]

top_panel <- ggplot(selected_rank, aes(x = axis, y = rank, group = 1)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  scale_y_reverse() +
  labs(x = NULL, y = "Rank in disease similarity") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

bottom_panel <- ggplot(selected_conn, aes(x = axis, y = connectivity)) +
  geom_violin(width = 0.8, trim = FALSE) +
  geom_jitter(width = 0.08, size = 1.4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  labs(x = "Signature", y = "AD DEGs connectivity") +
  theme_classic(base_size = 10)

if (requireNamespace("ggpubr", quietly = TRUE)) {
  comparisons <- lapply(axis_labels[axis_levels[-1]], function(x) c(axis_labels[["non_PGBV"]], x))
  bottom_panel <- bottom_panel +
    ggpubr::stat_compare_means(
      comparisons = comparisons,
      method = "wilcox.test",
      label = "p.signif",
      hide.ns = TRUE,
      size = 3
    )
}

if (requireNamespace("ggpubr", quietly = TRUE)) {
  combined <- ggpubr::ggarrange(top_panel, bottom_panel, ncol = 1, heights = c(1, 2), align = "v")
  ggsave(file.path(opt$out_dir, paste0(opt$selected_drug, "_rank_connectivity.pdf")),
         combined, width = 3.6, height = 5.0)
  ggsave(file.path(opt$out_dir, paste0(opt$selected_drug, "_rank_connectivity.png")),
         combined, width = 3.6, height = 5.0, dpi = 300)
} else {
  ggsave(file.path(opt$out_dir, paste0(opt$selected_drug, "_rank.pdf")),
         top_panel, width = 3.6, height = 1.8)
  ggsave(file.path(opt$out_dir, paste0(opt$selected_drug, "_connectivity.pdf")),
         bottom_panel, width = 3.6, height = 3.2)
}

data.table::fwrite(rank_summary, file.path(opt$out_dir, "cluster6_rank_summary.tsv"), sep = "\t")
data.table::fwrite(selected_rows, file.path(opt$out_dir, paste0(opt$selected_drug, "_scores.tsv")), sep = "\t")

message("Saved figures to: ", opt$out_dir)
