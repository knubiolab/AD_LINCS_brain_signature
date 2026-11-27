#!/usr/bin/env Rscript
# Comparative ligand RMSD / RMSF histograms (reviewer request)
#   Rscript plot_stability_hist.R
# Output: stability_hist.png, stability_hist.pdf, stability_summary.tsv
#   top row = ligand RMSD per target (200-300 ns); bottom row = per-residue RMSF for the same targets
suppressWarnings({

# Targets to plot, in order; the five named by the reviewer are the default.
KEEP <- c("BChE", "PPARg", "GR", "VDR", "HNF4a")

cfg <- list(
  list(t="BChE",   d="GROMACS_BCHE",                              l="Motesanib (rank 3)", c="#C0392B"),
  list(t="BChE",   d="GROMACS_BCHE_docked1/MD_rank1_300ns_v2/MD", l="Motesanib (rank 1)", c="#E67E22"),
  list(t="BChE",   d="GROMACS_4TPK",                              l="Reference ligand",   c="#2C7BB6"),
  list(t="PPARg",  d="GROMACS",                                   l="Motesanib",          c="#C0392B"),
  list(t="PPARg",  d="GROMACS_8REJ",                              l="Reference (LBB78)",  c="#2C7BB6"),
  list(t="GR",     d="GROMACS_NR3C1",                             l="Motesanib",          c="#C0392B"),
  list(t="GR",     d="GROMACS_3E7C",                              l="Reference (GSK866)", c="#2C7BB6"),
  list(t="VDR",    d="GROMACS_VDR_7qpp",                          l="Motesanib",          c="#C0392B"),
  list(t="VDR",    d="GROMACS_7QPP",                              l="Calcitriol",         c="#2C7BB6"),
  list(t="VDR",    d="GROMACS_VDR_Paricalcitol",                  l="Paricalcitol",       c="#4DAF4A"),
  list(t="HNF4a",  d="GROMACS_HNF4a",                             l="Motesanib",          c="#C0392B"),
  list(t="VEGFR2", d="GROMACS_VEGFR",                             l="Motesanib",          c="#C0392B")
)
BEG <- 200; END <- 300

read_xvg <- function(f) {
  if (!file.exists(f)) return(NULL)
  x <- readLines(f, warn = FALSE); x <- x[!grepl("^[#@]", x)]
  if (!length(x)) return(NULL)
  m <- do.call(rbind, lapply(strsplit(trimws(x), "[[:space:]]+"), function(z) as.numeric(z[1:2])))
  m <- m[stats::complete.cases(m), , drop = FALSE]
  if (!nrow(m)) return(NULL)
  data.frame(x = m[, 1], y = m[, 2])
}

# ---- collect data ----
rmsd <- list(); rmsf <- list(); rows <- list()
for (k in cfg) {
  a <- read_xvg(file.path(k$d, "rmsd_backbone_LIG.xvg"))
  if (!is.null(a)) { a <- a[a$x >= BEG & a$x <= END, ]
                     if (nrow(a) > 1) rmsd[[length(rmsd)+1]] <- c(k, list(v = a$y)) }
  b <- read_xvg(file.path(k$d, "rmsf.xvg"))
  if (!is.null(b) && nrow(b) > 1) rmsf[[length(rmsf)+1]] <- c(k, list(v = b$y))
  rows[[length(rows)+1]] <- data.frame(
    target = k$t, ligand = k$l, dir = k$d,
    rmsd_mean = if (!is.null(a) && nrow(a) > 1) round(mean(a$y), 3) else NA,
    rmsd_sd   = if (!is.null(a) && nrow(a) > 1) round(sd(a$y), 3)   else NA,
    rmsd_n    = if (!is.null(a)) nrow(a) else 0,
    rmsf_mean = if (!is.null(b) && nrow(b) > 1) round(mean(b$y), 3) else NA,
    rmsf_max  = if (!is.null(b) && nrow(b) > 1) round(max(b$y), 3)  else NA,
    stringsAsFactors = FALSE)
}
summ <- do.call(rbind, rows)
write.table(summ, "stability_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote: stability_summary.tsv\n"); print(summ, row.names = FALSE)

# Column order: only KEEP targets that actually have RMSD data (a target without a 300 ns run, e.g. VEGFR2, drops out)
have <- unique(sapply(rmsd, `[[`, "t"))
TG <- KEEP[KEEP %in% have]
if (!length(TG)) stop("no target has ligand RMSD data for 200-300 ns")
dropped <- setdiff(have, TG)
if (length(dropped)) cat("skipped (not in KEEP):", paste(dropped, collapse = ", "), "\n")

# ---- panels ----
panel <- function(items, xlab, main, nb = 30) {
  if (!length(items)) { plot.new(); title(main = paste(main, "(no data)")); return(invisible()) }
  rng <- range(unlist(lapply(items, `[[`, "v")))
  if (diff(rng) == 0) rng <- rng + c(-0.01, 0.01)
  br <- seq(rng[1], rng[2], length.out = nb + 1)
  hs <- lapply(items, function(it) hist(it$v, breaks = br, plot = FALSE))
  ymax <- max(sapply(hs, function(h) max(h$density)))
  plot(NA, xlim = rng, ylim = c(0, ymax * 1.10), xlab = xlab, ylab = "Density",
       main = main, las = 1, cex.main = 1.05, cex.axis = 0.85)
  for (i in seq_along(hs)) {
    cc <- items[[i]]$c
    rect(hs[[i]]$breaks[-length(hs[[i]]$breaks)], 0, hs[[i]]$breaks[-1], hs[[i]]$density,
         col = adjustcolor(cc, alpha.f = 0.38), border = adjustcolor(cc, alpha.f = 0.8))
  }
  legend("topright", legend = sapply(items, `[[`, "l"),
         fill = adjustcolor(sapply(items, `[[`, "c"), alpha.f = 0.5), bty = "n", cex = 0.62)
}

# ---- 2 rows x (number of targets) columns ----
draw <- function() {
  op <- par(mfrow = c(2, length(TG)), mar = c(4.2, 4.2, 2.8, 0.8), oma = c(0, 0, 2.2, 0))
  for (g in TG) panel(Filter(function(z) z$t == g, rmsd), "Ligand RMSD (nm)", g)
  for (g in TG) panel(Filter(function(z) z$t == g, rmsf), "Residue RMSF (nm)", "")
  mtext(sprintf("Ligand RMSD over %d-%d ns (top) and per-residue RMSF (bottom)", BEG, END),
        outer = TRUE, line = 0.4, cex = 1.05, font = 2)
  par(op)
}

W <- 330 * length(TG) + 130; H <- 950
png("stability_hist.png", width = W, height = H, res = 130); draw(); invisible(dev.off())
if (capabilities("cairo")) {
  cairo_pdf("stability_hist.pdf", width = W/130, height = H/130)
} else {
  pdf("stability_hist.pdf", width = W/130, height = H/130)
}
draw(); invisible(dev.off())
cat("wrote: stability_hist.png / stability_hist.pdf\n")

})
