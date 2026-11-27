# Run this script from the directory that holds one subdirectory per system
# (names as in targets.conf). Set WORKDIR to that path.
WORKDIR <- "~/"
setwd(WORKDIR)

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

# =========================
# 1) helpers
# =========================
read_xvg <- function(path) {
  lines <- readLines(path, warn = FALSE)
  dat_lines <- lines[!grepl("^\\s*[@#]", lines)]
  if (length(dat_lines) == 0) stop("No data lines in: ", path)
  
  dat <- read.table(text = dat_lines, header = FALSE)
  if (ncol(dat) < 2) stop("Need >=2 columns in: ", path)
  
  tibble(x = dat[[1]], y = dat[[2]])
}

convert_time_if_ps <- function(df) {
  # GROMACS time is often ps. If max looks like ps-scale, convert to ns.
  if (max(df$x, na.rm = TRUE) > 1000) {
    df$x <- df$x / 1000
    attr(df, "xlab") <- "Time (ns)"
  } else {
    attr(df, "xlab") <- "Time"
  }
  df
}

make_ts_plot <- function(fname, ylab, title, front = NULL) {
  d1 <- read_xvg(file.path("GROMACS", fname)) |>
    convert_time_if_ps() |>
    mutate(dataset = "GROMACS")
  xlab <- attr(d1, "xlab")
  
  d2 <- read_xvg(file.path("GROMACS_8REJ", fname)) |>
    convert_time_if_ps() |>
    mutate(dataset = "GROMACS_8REJ")
  
  dat <- bind_rows(d1, d2)
  
  # --- key: control draw order so the chosen system is plotted on top ---
  if (!is.null(front)) {
    dat <- dat |>
      mutate(dataset = factor(dataset, levels = c(setdiff(c("GROMACS","GROMACS_8REJ"), front), front))) |>
      arrange(dataset)
  }
  
  dat |>
    ggplot(aes(x = x, y = y, color = dataset)) +
    geom_line(linewidth = 0.7) +
    scale_color_manual(
      values = c("GROMACS" = "#FAD21C", "GROMACS_8REJ" = "#76BAFA"),
      breaks = c("GROMACS", "GROMACS_8REJ")  # fixed legend order
    ) +
    labs(title = title, x = xlab, y = ylab, color = NULL) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 11),
      legend.position = "bottom"
    )
}


make_rmsf_plot <- function(fname = "rmsf.xvg", res_min = 526, res_max = 776, front = NULL) {
  d1 <- read_xvg(file.path("GROMACS", fname)) |>
    mutate(dataset = "GROMACS") |>
    filter(x >= res_min, x <= res_max)
  
  d2 <- read_xvg(file.path("GROMACS_8REJ", fname)) |>
    mutate(dataset = "GROMACS_8REJ") |>
    filter(x >= res_min, x <= res_max)
  
  dat <- bind_rows(d1, d2)
  
  if (!is.null(front)) {
    dat <- dat |>
      mutate(dataset = factor(dataset, levels = c(setdiff(c("GROMACS","GROMACS_8REJ"), front), front))) |>
      arrange(dataset)
  }
  
  dat |>
    ggplot(aes(x = x, y = y, color = dataset)) +
    geom_line(linewidth = 0.7) +
    scale_color_manual(
      values = c("GROMACS" = "#FAD21C", "GROMACS_8REJ" = "#76BAFA"),
      breaks = c("GROMACS", "GROMACS_8REJ")
    ) +
    labs(title = "RMSF", x = "Residue", y = "RMSF (nm)", color = NULL) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 11),
      legend.position = "bottom"
    )
}


make_lig_sasa_plot <- function(front = NULL) {
  d1 <- read_xvg(file.path("GROMACS", "sasa_lig.xvg")) |>
    convert_time_if_ps() |>
    mutate(dataset = "GROMACS")
  xlab <- attr(d1, "xlab")
  
  d2 <- read_xvg(file.path("GROMACS_8REJ", "sasa_lig.xvg")) |>
    convert_time_if_ps() |>
    mutate(dataset = "GROMACS_8REJ")
  
  dat <- bind_rows(d1, d2)
  
  # put the chosen system in front
  if (!is.null(front)) {
    dat <- dat |>
      mutate(dataset = factor(dataset, levels = c(setdiff(c("GROMACS","GROMACS_8REJ"), front), front))) |>
      arrange(dataset)
  }
  
  dat |>
    ggplot(aes(x = x, y = y, color = dataset)) +
    geom_line(linewidth = 0.7) +
    scale_color_manual(
      values = c("GROMACS" = "#FAD21C", "GROMACS_8REJ" = "#76BAFA"),
      breaks = c("GROMACS", "GROMACS_8REJ")
    ) +
    labs(title = "Ligand SASA", x = xlab, y = expression("SASA (nm"^2*")"), color = NULL) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 11),
      legend.position = "bottom"
    )
}


# =========================
# 2) build 6 panels (order fixed)
# =========================
p1 <- make_ts_plot("rmsd_backbone_LIG.xvg", "RMSD (nm)", "Ligand RMSD", front = "GROMACS")
p2 <- make_ts_plot("rmsd_backbone.xvg",     "RMSD (nm)", "Backbone RMSD")
p3 <- make_rmsf_plot("rmsf.xvg", 200, 480, front = "GROMACS")
p4 <- make_ts_plot("hb.xvg",              "count", "H-bonds")
p5 <- make_ts_plot("gyrate1.xvg",         "Rg (nm)",         "Radius of Gyration", front = "GROMACS")
p6 <- make_lig_sasa_plot(front = "GROMACS")
p7 <- make_ts_plot("hb_essential.xvg",              "count", "H-bonds to key residues", front = "GROMACS")

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
pdf("R/Figures.Positive_control_drug_revised/PPARG.pdf", height = 7.96, width = 9)
combined
dev.off()
