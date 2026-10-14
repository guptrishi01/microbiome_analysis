#Author: Rishi Gupta / Claude
#Date: 2026
#Description: One-page PDF section divider for the combined comparison report.
#Usage: Rscript titlePage.R <out.pdf> <title> <subtitle_line1> [<subtitle_line2> ...]

args <- commandArgs(trailingOnly = TRUE)
outFile <- args[1]
title <- args[2]
lines <- args[-c(1,2)]

pdf(outFile, width = 11, height = 8.5)
plot.new()
text(0.5, 0.65, title, cex = 1.8, font = 2)
yy <- 0.5
for (ln in lines) {
  text(0.5, yy, ln, cex = 1.05)
  yy <- yy - 0.06
}
dev.off()
cat("wrote", outFile, "\n")
