#!/usr/bin/env Rscript
# Supplementary Fig. 6a: ligand cluster occupancy per target-ligand pair, with the 50% cumulative mark
#   Rscript plot_cluster_fraction.R
# Output: cluster_fraction.png, cluster_fraction.pdf, cluster_fraction.tsv
suppressWarnings({
OUTDIR <- "rep_poses"; NSHOW <- 6; THR <- 50
cfg <- list(
  c("BChE","Motesanib (contact-selected)","GROMACS_BCHE"),
  c("BChE","Motesanib (top-scoring)","GROMACS_BCHE_docked1/MD_rank1_300ns_v2/MD"),
  c("BChE","Reference ligand","GROMACS_4TPK"),
  c("PPARg","Motesanib","GROMACS"),
  c("PPARg","Reference ligand","GROMACS_8REJ"),
  c("GR","Motesanib","GROMACS_NR3C1"),
  c("GR","Reference ligand","GROMACS_3E7C"),
  c("VDR","Motesanib","GROMACS_VDR_7qpp"),
  c("VDR","Reference ligand","GROMACS_7QPP"),
  c("VDR","Paricalcitol","GROMACS_VDR_Paricalcitol"),
  c("HNF4a","Motesanib","GROMACS_HNF4a")
)
rows <- list(); labs <- character(); mat <- NULL
for (k in cfg) {
  f <- file.path(OUTDIR, k[3], "clust_size.xvg")
  if (!file.exists(f)) { message("missing: ", f); next }
  x <- readLines(f, warn = FALSE); x <- x[!grepl("^[#@]", x)]
  v <- suppressWarnings(as.numeric(sapply(strsplit(trimws(x), "[[:space:]]+"), `[`, 2)))
  v <- v[!is.na(v)]; if (!length(v)) next
  p <- sort(100 * v / sum(v), decreasing = TRUE)
  need <- which(cumsum(p) >= THR)[1]
  rows[[length(rows)+1]] <- data.frame(target = k[1], ligand = k[2],
      n_clusters = length(p), top1 = round(p[1],1),
      top2 = round(ifelse(length(p)>1, p[2], NA),1),
      n_to_50 = need, stringsAsFactors = FALSE)
  pp <- rep(0, NSHOW); m <- min(NSHOW, length(p)); pp[1:m] <- p[1:m]
  if (length(p) > NSHOW) pp[NSHOW] <- pp[NSHOW] + sum(p[(NSHOW+1):length(p)])
  mat <- rbind(mat, pp); labs <- c(labs, paste0(k[1], " – ", k[2]))
}
summ <- do.call(rbind, rows)
write.table(summ, "cluster_fraction.tsv", sep="\t", quote=FALSE, row.names=FALSE)
print(summ, row.names = FALSE)

cols <- c("#2C7BB6","#7FB8DA","#BFD9EC","#E8E8E8","#D9D9D9","#BDBDBD")

# the drawing code is wrapped in a function so the png and pdf devices get identical output
draw <- function() {
  op <- par(mar = c(4.5, 15, 5.5, 2))
  bp <- barplot(t(mat[nrow(mat):1, , drop=FALSE]), horiz = TRUE, beside = FALSE,
          names.arg = rev(labs), las = 1, xlim = c(0, 100), col = cols,
          border = "white", xlab = "Fraction of analysed frames (%)",
          main = "", cex.names = 0.8)
  title(main = "Ligand cluster occupancy, 200–300 ns", line = 3.8, cex.main = 1.15)
  legend(x = 50, y = par("usr")[4] + 0.05 * diff(par("usr")[3:4]),
         xjust = 0.5, yjust = 0, horiz = TRUE, xpd = NA, bty = "n", cex = 0.75,
         legend = c("cluster 1","cluster 2","cluster 3","4","5","≥6"), fill = cols)
  abline(v = THR, col = "red", lwd = 2, lty = 2)
  text(THR, par("usr")[4], " 50%", col = "red", adj = c(0, 1.3), cex = 0.85, xpd = NA)
  text(mat[nrow(mat):1, 1] / 2, bp, rev(sprintf("%.0f%%", mat[, 1])),
       cex = 0.72, col = "white")
  par(op)
}

png("cluster_fraction.png", width = 1500, height = 950, res = 150)
draw(); invisible(dev.off())

# cairo_pdf handles UTF-8 glyphs (en dash, >=) and font embedding; fall back to pdf() if unavailable.
if (capabilities("cairo")) {
  cairo_pdf("cluster_fraction.pdf", width = 10, height = 6.33)
} else {
  pdf("cluster_fraction.pdf", width = 10, height = 6.33)
}
draw(); invisible(dev.off())

cat("wrote: cluster_fraction.png / cluster_fraction.pdf / cluster_fraction.tsv\n")
})
