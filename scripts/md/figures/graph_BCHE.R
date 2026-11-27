# Run this script from the directory that holds one subdirectory per system
# (names as in targets.conf). Set WORKDIR to that path.
WORKDIR <- "~/"
setwd(WORKDIR)

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

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
  d1 <- read_xvg(file.path("GROMACS_BCHE", fname)) |>
    convert_time_if_ps() |>
    mutate(dataset = "GROMACS_BCHE")
  xlab <- attr(d1, "xlab")
  
  d2 <- read_xvg(file.path("GROMACS_4TPK", fname)) |>
    convert_time_if_ps() |>
    mutate(dataset = "GROMACS_4TPK")
  
  dat <- bind_rows(d1, d2)
  
  # --- key: control draw order so the chosen system is plotted on top ---
  if (!is.null(front)) {
    dat <- dat |>
      mutate(dataset = factor(dataset, levels = c(setdiff(c("GROMACS_BCHE","GROMACS_4TPK"), front), front))) |>
      arrange(dataset)
  }
  
  dat |>
    ggplot(aes(x = x, y = y, color = dataset)) +
    geom_line(linewidth = 0.7) +
    scale_color_manual(
      values = c("GROMACS_BCHE" = "#FAD21C", "GROMACS_4TPK" = "#76BAFA"),
      breaks = c("GROMACS_BCHE", "GROMACS_4TPK")  # fixed legend order
    ) +
    labs(title = title, x = xlab, y = ylab, color = NULL) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 11),
      legend.position = "bottom"
    )
}


make_rmsf_plot <- function(fname = "rmsf.xvg", res_min = 526, res_max = 776, front = NULL) {
  d1 <- read_xvg(file.path("GROMACS_BCHE", fname)) |>
    mutate(dataset = "GROMACS_BCHE") |>
    filter(x >= res_min, x <= res_max)
  
  d2 <- read_xvg(file.path("GROMACS_4TPK", fname)) |>
    mutate(dataset = "GROMACS_4TPK") |>
    filter(x >= res_min, x <= res_max)
  
  dat <- bind_rows(d1, d2)
  
  if (!is.null(front)) {
    dat <- dat |>
      mutate(dataset = factor(dataset, levels = c(setdiff(c("GROMACS_BCHE","GROMACS_4TPK"), front), front))) |>
      arrange(dataset)
  }
  
  dat |>
    ggplot(aes(x = x, y = y, color = dataset)) +
    geom_line(linewidth = 0.7) +
    scale_color_manual(
      values = c("GROMACS_BCHE" = "#FAD21C", "GROMACS_4TPK" = "#76BAFA"),
      breaks = c("GROMACS_BCHE", "GROMACS_4TPK")
    ) +
    labs(title = "RMSF", x = "Residue", y = "RMSF (nm)", color = NULL) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 11),
      legend.position = "bottom"
    )
}


make_overlap_bin2d <- function(df, bins = 60,
                               col_bche = "#FAD21C",
                               col_4tpk = "#76BAFA",
                               alpha_each = 1) {
  
  # common bin edges so both systems fall on the same grid
  xbreaks <- seq(min(df$d_cas, na.rm=TRUE), max(df$d_cas, na.rm=TRUE), length.out = bins + 1)
  ybreaks <- seq(min(df$d_pas, na.rm=TRUE), max(df$d_pas, na.rm=TRUE), length.out = bins + 1)
  min <- min(c(df$d_cas, df$d_pas))
  max <- max(c(df$d_cas, df$d_pas))
  
  # 2D bin counts, per system
  binned <- df %>%
    filter(system %in% c("GROMACS_BCHE", "GROMACS_4TPK")) %>%
    mutate(
      xb = cut(d_cas, breaks = xbreaks, include.lowest = TRUE, right = FALSE),
      yb = cut(d_pas, breaks = ybreaks, include.lowest = TRUE, right = FALSE)
    ) %>%
    count(system, xb, yb, name = "n") %>%
    group_by(system) %>%            # scale within each system
    mutate(n_scaled = n / max(n)) %>%
    ungroup()
  
  # bin centre coordinates, for geom_tile()
  xmid <- (xbreaks[-1] + xbreaks[-length(xbreaks)]) / 2
  ymid <- (ybreaks[-1] + ybreaks[-length(ybreaks)]) / 2
  x_levels <- levels(cut(xmid, breaks = xbreaks, include.lowest = TRUE, right = FALSE))
  y_levels <- levels(cut(ymid, breaks = ybreaks, include.lowest = TRUE, right = FALSE))
  
  binned <- binned %>%
    mutate(
      x = xmid[match(as.character(xb), x_levels)],
      y = ymid[match(as.character(yb), y_levels)]
    )
  
  p <- ggplot() +
    # 4TPK layer
    geom_tile(
      data = filter(binned, system == "GROMACS_4TPK"),
      aes(x = x, y = y, alpha = n_scaled),
      fill = col_4tpk, color = NA
    ) +
    # BCHE layer
    geom_tile(
      data = filter(binned, system == "GROMACS_BCHE"),
      aes(x = x, y = y, alpha = n_scaled),
      fill = col_bche, color = NA
    ) +
    scale_alpha(range = c(0.5, alpha_each), guide = "none") +
    labs(
      title = "State landscape",
      x = "d(LIG–CAS) (nm)",
      y = "d(LIG–PAS) (nm)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 11),
      legend.position = "bottom"
    ) + lims(x = c(min, max), y = c(min, max))
  
  p
}

make_lig_sasa_plot <- function(front = NULL) {
  d1 <- read_xvg(file.path("GROMACS_BCHE", "sasa_lig.xvg")) |>
    convert_time_if_ps() |>
    mutate(dataset = "GROMACS_BCHE")
  xlab <- attr(d1, "xlab")
  
  d2 <- read_xvg(file.path("GROMACS_4TPK", "sasa_lig.xvg")) |>
    convert_time_if_ps() |>
    mutate(dataset = "GROMACS_4TPK")
  
  dat <- bind_rows(d1, d2)
  
  # put the chosen system in front
  if (!is.null(front)) {
    dat <- dat |>
      mutate(dataset = factor(dataset, levels = c(setdiff(c("GROMACS_BCHE","GROMACS_4TPK"), front), front))) |>
      arrange(dataset)
  }
  
  dat |>
    ggplot(aes(x = x, y = y, color = dataset)) +
    geom_line(linewidth = 0.7) +
    scale_color_manual(
      values = c("GROMACS_BCHE" = "#FAD21C", "GROMACS_4TPK" = "#76BAFA"),
      breaks = c("GROMACS_BCHE", "GROMACS_4TPK")
    ) +
    labs(title = "Ligand SASA", x = xlab, y = expression("SASA (nm"^2*")"), color = NULL) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 11),
      legend.position = "bottom"
    )
}

read_cas_pas <- function(dir, system,
                         cas_file = "dist_LIG_CAS_com.xvg",
                         pas_file = "dist_LIG_PAS_com.xvg") {
  
  cas <- read_xvg(file.path(dir, cas_file)) %>%
    rename(time = x, d_cas = y)
  
  pas <- read_xvg(file.path(dir, pas_file)) %>%
    rename(time = x, d_pas = y)
  
  df <- inner_join(cas, pas, by = "time") %>%
    arrange(time)
  
  # convert time once (based on CAS)
  tmp <- tibble(x = df$time, y = df$d_cas) %>% convert_time_if_ps()
  xlab <- attr(tmp, "xlab")
  
  out <- df %>%
    mutate(time = tmp$x, system = system)
  
  attr(out, "xlab") <- xlab
  out
}

cols <- c(
  "GROMACS_BCHE" = "#FAD21C",
  "GROMACS_4TPK" = "#76BAFA"
)

# residue groups
group_def <- list(
  "TRP82" = c("TRP82"),
  "Oxyanion-hole" = c("GLY116","GLY117"),
  "Catalytic triad" = c("SER198","HIS438","GLU325"),
  "PAS" = c("ASP70","TYR332")
)

# -----------------------------
# per dataset: residue mindist -> group mindist (minimum over the group, per frame)
# -----------------------------
load_dataset_group_mindist <- function(dir, dataset, group_def, warmup_ns = 0) {
  
  residues <- unique(unlist(group_def))
  
  # join the per-residue mindist series into wide form, keyed on time
  wide <- residues %>%
    set_names() %>%
    map(~{
      f <- file.path(dir, paste0("mindist_", .x, ".xvg"))
      read_xvg(f) %>%
        convert_time_if_ps() %>%
        rename(time = x, dist = y) %>%
        dplyr::filter(time >= warmup_ns) %>%
        transmute(time, !!.x := dist)
    }) %>%
    reduce(full_join, by = "time") %>%
    arrange(time)
  
  # minimum distance within each group, frame by frame
  group_long <- imap_dfr(group_def, function(res_list, gname) {
    # pmin() must tolerate NA in any of the res_list columns
    cols_exist <- intersect(res_list, colnames(wide))
    if (length(cols_exist) == 0) stop("No residue columns found for group: ", gname)
    
    dmin <- do.call(pmin, c(wide[cols_exist], na.rm = TRUE))
    
    tibble(
      dataset = dataset,
      group   = gname,
      dist    = dmin
    ) %>% filter(is.finite(dist))
  })
  
  group_long
}

dir_bche <- "GROMACS_BCHE"
dir_4tpk <- "GROMACS_4TPK"

df_bche <- read_cas_pas("GROMACS_BCHE", "GROMACS_BCHE")
df_4tpk <- read_cas_pas("GROMACS_4TPK",  "GROMACS_4TPK")

df <- bind_rows(df_bche, df_4tpk)

dat <- bind_rows(
  load_dataset_group_mindist(dir_bche, "GROMACS_BCHE", group_def),
  load_dataset_group_mindist(dir_4tpk, "GROMACS_4TPK", group_def)
)

dat <- dat %>%
  mutate(group = factor(group, levels = names(group_def)))

# -----------------------------
# Plot
# -----------------------------

p7 <- ggplot(dat, aes(x = dist, y = after_stat(scaled), color = dataset, fill = dataset)) +
  geom_density(alpha = 0.25, linewidth = 0.8, adjust = 1) +
  facet_grid(group ~ ., scales = "free_y", switch = "y") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  labs(x = "Distance (nm)", color = NULL, fill = NULL) +
  theme_classic() +
  ggtitle("Minimum distance") +
  theme(legend.position = "none",
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 360), 
        plot.title = element_text(size = 11), 
        axis.title.y = element_blank()) +
  geom_vline(xintercept = 0.35, color = "red")

read_mindist_dir <- function(dir, system_name) {
  files <- list.files(dir, pattern = "^mindist_.*\\.xvg$", full.names = TRUE)
  if (length(files) == 0) stop("No mindist_*.xvg in: ", dir)
  
  map_dfr(files, function(f) {
    residue <- gsub("^mindist_|\\.xvg$", "", basename(f))  # file name encodes the residue
    read_xvg(f) |>
      convert_time_if_ps() |>
      transmute(
        time = x,
        dist = y,
        residue = residue,
        system = system_name
      )
  })
}

# 1) read both systems
d_bche <- read_mindist_dir(dir_bche, "BCHE")
d_4tpk <- read_mindist_dir(dir_4tpk, "4TPK")

# 2) keep only residues present in both systems
common_res <- intersect(unique(d_bche$residue), unique(d_4tpk$residue))
mindist_all <- bind_rows(d_bche, d_4tpk) |>
  filter(residue %in% common_res)

# 3) order residues by residue number
mindist_all <- mindist_all |>
  mutate(
    resnum = readr::parse_number(residue),
    residue = fct_reorder(residue, resnum)
  )

cutoff <- 0.35   # nm
t_split <- 200   # ns: before/after boundary (must be defined before occ_res)
t0 <- 200         # ns: start of the analysis window
occ_res <- mindist_all |>
  mutate(window = case_when(
    time < t_split ~ "before",
    TRUE ~ "after"
  )) |>
  group_by(system, residue, window) |>
  summarise(
    n = n(),
    occ = mean(dist < cutoff, na.rm = TRUE) * 100,
    mean_dist = mean(dist, na.rm = TRUE),
    median_dist = median(dist, na.rm = TRUE),
    .groups = "drop"
  ) |>
  # pivot before/after to wide form
  pivot_wider(
    names_from  = window,
    values_from = c(n, occ, mean_dist, median_dist)
  ) |>
  # add the full-trajectory value alongside
  left_join(
    mindist_all |>
      group_by(system, residue) |>
      summarise(
        n_all = n(),
        occ_all = mean(dist < cutoff, na.rm = TRUE) * 100,
        mean_all = mean(dist, na.rm = TRUE),
        median_all = median(dist, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c("system", "residue")
  ) |>
  arrange(desc(occ_after), desc(occ_before), residue)

occ_res %>% arrange(residue)

# residue -> group mapping
res2group <- tribble(
  ~residue,  ~group,
  "ASP70",   "PAS",
  "TYR332",  "PAS",
  "TRP82",   "TRP82",
  "GLY116",  "OXYANION",
  "GLY117",  "OXYANION",
  "SER198",  "TRIAD",
  "HIS438",  "TRIAD",
  "GLU325",  "TRIAD",
  "TRP231",  "CAS",
  "LEU286",  "CAS",
  "PHE329",  "CAS",
  "PHE398",  "CAS"
)

occ_res$residue <- sub("\\..*$", "", as.character(occ_res$residue))
# 1) occ_res -> long
occ_long <- occ_res %>%
  dplyr::filter(!residue %in% c("045", "GLY118", "GLU325")) %>%
  dplyr::select(system, residue, occ_before, occ_after) %>%
  pivot_longer(
    cols = c(occ_before, occ_after),
    names_to = "window",
    values_to = "occupancy"
  ) %>%
  mutate(
    window = recode(window, occ_before = "before", occ_after = "after"),
    window = factor(window, levels = c("before", "after")),
    residue = toupper(residue)
  ) %>%
  left_join(res2group, by = "residue") %>%
  mutate(
    group = replace_na(group, "Other"),
    group = factor(group, levels = c("TRP82","OXYANION","TRIAD","PAS","CAS","Other"))
  )

cols_bar <- c(
  "BCHE" = "#FAD21C",
  "4TPK" = "#76BAFA"
)

# 2) bar plot: x=before/after, fill=system, facet=group~residue
library(ggh4x)

occ_long$window <- factor(occ_long$window, levels = c("after", "before"))
p8 <- ggplot(occ_long, aes(y = occupancy, x = window, fill = system)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = cols_bar) +
  facet_nested(~ group + residue + system, scales = "free", space = "free", nest_line = element_line()) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.05))) +
  labs(y = NULL, x = "Occupancy (<0.35 nm, %)", fill = NULL) +
  theme_classic() +
  geom_hline(yintercept = 70, color = "red") +
  ggtitle("Contact Occupancy (Before vs After 200 ns)") +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    strip.background = element_blank(), 
    title = element_text(size = 11), 
    axis.text.x = element_text(angle = 90, color = "black"), 
    axis.text.y = element_text(angle = 90, color = "black")
  )
p8

# =========================
# 2) build 6 panels (order fixed)
# =========================
p1 <- make_ts_plot("rmsd_backbone_LIG.xvg", "RMSD (nm)", "Ligand RMSD", front = "GROMACS_BCHE")
p2 <- make_ts_plot("rmsd_backbone.xvg",     "RMSD (nm)", "Backbone RMSD")                 # both systems, no reordering
p3 <- make_rmsf_plot("rmsf.xvg", 2, 530, front = "GROMACS_BCHE")

p4 <- make_overlap_bin2d(df, bins = 60,
                         col_bche = "#FAD21C",
                         col_4tpk = "#76BAFA",
                         alpha_each = 1)
p5 <- make_ts_plot("gyrate1.xvg",         "Rg (nm)",         "Radius of Gyration", front = "GROMACS_BCHE")
p6 <- make_lig_sasa_plot(front = "GROMACS_BCHE")

inner_gap <- 22  # pt: gap between the two columns
outer_gap <- 6   # pt: outer margin

left_col_margin  <- theme(plot.margin = margin(outer_gap, inner_gap, outer_gap, outer_gap, unit = "pt"))
right_col_margin <- theme(plot.margin = margin(outer_gap, outer_gap, outer_gap, inner_gap, unit = "pt"))

p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")
p4 <- p4 + theme(legend.position = "none")
p5 <- p5 + theme(legend.position = "none")
p6 <- p6 + theme(legend.position = "none")

# left column: p1, p3, p5 (wider right margin)
p1 <- p1 + left_col_margin
p3 <- p3 + left_col_margin
# p3 <- patchwork::wrap_elements(full = ggplot2::ggplotGrob(p3))
p5 <- p5 + left_col_margin

# right column: p2, p4, p6 (wider left margin)
p2 <- p2 + right_col_margin
p4 <- p4 + right_col_margin
p6 <- p6 + right_col_margin
p7 <- p7 + right_col_margin
blank_right <- plot_spacer() & theme_void() & right_col_margin

combined <- (p1 + p2 + p7 + p3 + p4 + blank_right + p5 + p6 + blank_right) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")


dir.create("R/Figures.Positive_control_drug_revised", recursive = TRUE, showWarnings = FALSE)
pdf("R/Figures.Positive_control_drug_revised/BCHE.pdf", height = 7.96, width = 9)
combined
dev.off()

pdf("R/Figures.Positive_control_drug_revised/BCHE_time.pdf", height = 3, width = 9)
p8
dev.off()
