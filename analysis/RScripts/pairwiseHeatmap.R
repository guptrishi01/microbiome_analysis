#Author: Rishi Gupta / Claude
#Date: 2026
#Description: Pairwise study-timepoint heatmap from a run's Genus_PairwiseComparison.txt.
#Cell colour = Spearman r between signed -log10(p) vectors; cell text = r and the p-value.
#The p-value is the BH-adjusted one compareStudies_16S.R writes (adjusted across all pairs
#in that run), not the raw cor.test p.
#Usage: Rscript pairwiseHeatmap.R <Genus_PairwiseComparison.txt> <out_prefix>

args <- commandArgs(trailingOnly = TRUE)
pairs <- read.table(args[1], sep = "\t", header = TRUE, quote = "\"", stringsAsFactors = FALSE)

#studyPairs look like "RYGB_SaoPaulo-T2_Lazaro_GB-pPost": study names contain underscores,
#but the study/timepoint separator is a hyphen and timepoint suffixes contain no underscores
parts <- regmatches(pairs$studyPairs, regexec("^([^-]+-[^_]+)_(.+-[^_]+)$", pairs$studyPairs))
pairs$a <- sapply(parts, `[`, 2)
pairs$b <- sapply(parts, `[`, 3)
labels <- unique(c(pairs$a, pairs$b))
n <- length(labels)

R <- matrix(NA_real_, n, n, dimnames = list(labels, labels))
P <- R
for (i in seq_len(nrow(pairs))) {
  ia <- match(pairs$a[i], labels)
  ib <- match(pairs$b[i], labels)
  R[ia, ib] <- R[ib, ia] <- pairs$r[i]
  P[ia, ib] <- P[ib, ia] <- pairs$pval[i]
}

fmtP <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
study <- sub("-.*", "", labels)
sameStudy <- outer(study, study, "==")

lo <- min(R, na.rm = TRUE)
hi <- max(R, na.rm = TRUE)
pal <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)
colIdx <- function(v) pmin(101, pmax(1, round((v - lo) / (hi - lo) * 100) + 1))

pdf(paste0(args[2], "_heatmap.pdf"), width = 11, height = 11)
par(mar = c(9, 9, 4, 2))
plot(NA, xlim = c(0.5, n + 0.5), ylim = c(0.5, n + 0.5), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", asp = 1, bty = "n")
for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    x <- j
    y <- n - i + 1
    if (i == j) {
      rect(x - 0.5, y - 0.5, x + 0.5, y + 0.5, col = "grey90", border = "white")
      next
    }
    fill <- if (sameStudy[i, j]) "grey85" else pal[colIdx(R[i, j])]
    rect(x - 0.5, y - 0.5, x + 0.5, y + 0.5, col = fill, border = "white")
    text(x, y, sprintf("r=%.2f\np=%s", R[i, j], fmtP(P[i, j])), cex = 0.55)
  }
}
axis(1, at = seq_len(n), labels = labels, las = 2, cex.axis = 0.7)
axis(2, at = seq_len(n), labels = rev(labels), las = 2, cex.axis = 0.7)
title(main = "Pairwise Spearman r (signed -log10 p) between study-timepoints", cex.main = 1)
mtext(paste0("Cell text: r and BH-adjusted p. Grey = same study (excluded from headline r). ",
             "Colour: ", sprintf("%.2f", lo), " (blue) to ", sprintf("%.2f", hi), " (red)"),
      side = 3, line = 0.2, cex = 0.7)
dev.off()

write.csv(R, paste0(args[2], "_r.csv"))
write.csv(P, paste0(args[2], "_p_BHadjusted.csv"))
cat("wrote", paste0(args[2], "_heatmap.pdf"), "and r/p csv tables for", n, "study-timepoints\n")
