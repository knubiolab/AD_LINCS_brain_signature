# Run this script from the directory that holds one subdirectory per system
# (names as in targets.conf). Set WORKDIR to that path.
WORKDIR <- "~/"
setwd(WORKDIR)

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

# =========================
# shared setup: three datasets
# =========================
datasets3 <- c("GROMACS_VDR_7qpp", "GROMACS_7QPP", "GROMACS_VDR_Paricalcitol")

cols3 <- c(
  "GROMACS_VDR_7qpp"          = "#FAD21C",
  "GROMACS_7QPP"              = "#76BAFA",
  "GROMACS_VDR_Paricalcitol"  = "#2A9D8F"
)

read_xvg <- function(path) {
  lines <- readLines(path, warn = FALSE)
  dat_lines <- lines[!grepl("^\\s*[@#]", lines)]
  if (length(dat_lines) == 0) stop("No data lines in: ", path)
  
  dat <- read.table(text = dat_lines, header = FALSE)
  if (ncol(dat) < 2) stop("Need >=2 columns in: ", path)
  
  tibble(x = dat[[1]], y = dat[[2]])
}

convert_time_if_ps <- function(df) {
  if (max(df$x, na.rm = TRUE) > 1000) {
    df$x <- df$x / 1000
    attr(df, "xlab") <- "Time (ns)"
  } else {
    attr(df, "xlab") <- "Time"
  }
  df
}

load_one <- function(dir, fname, label, do_time_convert = TRUE) {
  d <- read_xvg(file.path(dir, fname))
  if (do_time_convert) d <- convert_time_if_ps(d)
  d %>% dplyr::mutate(dataset = label)
}

# =========================
# 1) time-series plots (three datasets)
# =========================
make_ts_plot <- function(fname, ylab, title, front = NULL) {
  d1 <- load_one("GROMACS_VDR_7qpp",         fname, "GROMACS_VDR_7qpp",         do_time_convert = TRUE)
  xlab <- attr(d1, "xlab")
  d2 <- load_one("GROMACS_7QPP",             fname, "GROMACS_7QPP",             do_time_convert = TRUE)
  d3 <- load_one("GROMACS_VDR_Paricalcitol", fname, "GROMACS_VDR_Paricalcitol", do_time_convert = TRUE)
  
  dat <- dplyr::bind_rows(d1, d2, d3)
  
  if (!is.null(front)) {
    dat <- dat %>%
      dplyr::mutate(dataset = factor(dataset, levels = c(setdiff(datasets3, front), front))) %>%
      dplyr::arrange(dataset)
  }
  
  ggplot2::ggplot(dat, ggplot2::aes(x = x, y = y, color = dataset)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::scale_color_manual(values = cols3, breaks = datasets3) +
    ggplot2::labs(title = title, x = xlab, y = ylab, color = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 11),
      legend.position = "bottom"
    )
}

# =========================
# 2) RMSF window facets (three datasets)
# =========================
make_rmsf_plot <- function(
    fname   = "rmsf.xvg",
    res_min = 120,
    res_max = 430,
    datasets = datasets3,
    cols     = cols3,
    front    = NULL
) {
  # read RMSF from each dataset directory
  dat_list <- lapply(datasets, function(ds) {
    read_xvg(file.path(ds, fname)) |>
      dplyr::mutate(dataset = ds) |>
      dplyr::filter(x >= res_min, x <= res_max)
  })
  
  dat <- dplyr::bind_rows(dat_list)
  
  # choose which dataset is drawn on top of the overlay
  if (!is.null(front)) {
    dat <- dat |>
      dplyr::mutate(dataset = factor(dataset, levels = c(setdiff(datasets, front), front))) |>
      dplyr::arrange(dataset)
  } else {
    dat <- dat |>
      dplyr::mutate(dataset = factor(dataset, levels = datasets))
  }
  
  ggplot2::ggplot(dat, ggplot2::aes(x = x, y = y, color = dataset)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::scale_color_manual(values = cols, breaks = datasets) +
    ggplot2::labs(title = "RMSF", x = "Residue", y = "RMSF (nm)", color = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 11),
      legend.position = "bottom"
    )
}

# =========================
# 3) ligand SASA (three datasets)
# =========================
make_lig_sasa_plot <- function(front = NULL) {
  d1 <- load_one("GROMACS_VDR_7qpp",         "sasa_lig.xvg", "GROMACS_VDR_7qpp",         do_time_convert = TRUE)
  xlab <- attr(d1, "xlab")
  d2 <- load_one("GROMACS_7QPP",             "sasa_lig.xvg", "GROMACS_7QPP",             do_time_convert = TRUE)
  d3 <- load_one("GROMACS_VDR_Paricalcitol", "sasa_lig.xvg", "GROMACS_VDR_Paricalcitol", do_time_convert = TRUE)
  
  dat <- dplyr::bind_rows(d1, d2, d3)
  
  if (!is.null(front)) {
    dat <- dat %>%
      dplyr::mutate(dataset = factor(dataset, levels = c(setdiff(datasets3, front), front))) %>%
      dplyr::arrange(dataset)
  }
  
  ggplot2::ggplot(dat, ggplot2::aes(x = x, y = y, color = dataset)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::scale_color_manual(values = cols3, breaks = datasets3) +
    ggplot2::labs(title = "Ligand SASA", x = xlab, y = expression("SASA (nm"^2*")"), color = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 11),
      legend.position = "bottom"
    )
}

# =========================
# 4) 6 panels
# =========================
p1 <- make_ts_plot("rmsd_backbone_LIG.xvg", "RMSD (nm)", "Ligand RMSD", front = "GROMACS_VDR_7qpp")
p2 <- make_ts_plot("rmsd_backbone.xvg",     "RMSD (nm)", "Backbone RMSD")
p3 <- make_rmsf_plot(res_min = 118, res_max = 425, front = "GROMACS_VDR_Paricalcitol")
p4 <- make_ts_plot("hb.xvg",          "count", "H-bonds", front = "GROMACS_VDR_7qpp")
p5 <- make_ts_plot("gyrate1.xvg",     "Rg (nm)",         "Radius of Gyration", front = "GROMACS_VDR_7qpp")
p6 <- make_lig_sasa_plot(front = "GROMACS_7QPP")
p7 <- make_ts_plot("hb_essential.xvg",              "count", "H-bonds to key residues", front = "GROMACS_VDR_7qpp")

inner_gap <- 22  # pt: gap between the two columns
outer_gap <- 6   # pt: outer margin

left_col_margin  <- theme(plot.margin = margin(outer_gap, inner_gap, outer_gap, outer_gap, unit = "pt"))
right_col_margin <- theme(plot.margin = margin(outer_gap, outer_gap, outer_gap, inner_gap, unit = "pt"))

# left column: p1, p3, p5 (wider right margin)
p1 <- p1 + left_col_margin
p3 <- p3 + left_col_margin
p5 <- p5 + left_col_margin

# right column: p2, p4, p6 (wider left margin)
p2 <- p2 + right_col_margin
p4 <- p4 + right_col_margin
p6 <- p6 + right_col_margin
blank_right <- plot_spacer() & theme_void() & right_col_margin

combined <- (p1 + p2 + p7 + p3 + p4 + blank_right + p5 + p6 + blank_right) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

combined

dir.create("R/Figures.Positive_control_drug_revised", recursive = TRUE, showWarnings = FALSE)
pdf("R/Figures.Positive_control_drug_revised/VDR.pdf", height = 7.96, width = 9)
combined
dev.off()
