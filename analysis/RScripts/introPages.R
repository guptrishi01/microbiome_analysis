#Author: Rishi Gupta / Claude
#Date: 2026
#Description: Report cover pages -- title, the full cross-phase dataset table, and a
#plain-language explanation of how the heatmaps/scatterplots are built. No statistics
#computed here; purely descriptive front matter for the combined report.
#Usage: Rscript introPages.R <out.pdf>

args <- commandArgs(trailingOnly = TRUE)
outFile <- args[1]

pdf(outFile, width = 11, height = 8.5)

# ---- Page 1: title ----
plot.new()
text(0.5, 0.6, "Gut Microbiome Signature Correlation Analysis", cex=1.8, font=2)
text(0.5, 0.52, "Pairwise Comparison Report", cex=1.3)
text(0.5, 0.40, "Reproduction (Phase 1) + extension to non-bariatric GI surgery (Phase 2)", cex=0.9)
text(0.5, 0.35, "+ bariatric positive control, RYGB-specific and non-RYGB sub-phases (Phase 3)", cex=0.9)

# ---- Page 2: all datasets table ----
plot.new()
par(family = "mono")
text(0.5, 0.97, "All Datasets Used In This Analysis", cex=1.5, font=2, family="sans")
text(0.5, 0.925, "12 cohorts across 3 phases, all 16S rRNA, all repeat-sampled pre/post surgery", cex=0.85, family="sans")

tbl <- rbind(
 c("Cohort","Phase","Accession","Procedure","Subjects","Runs"),
 c("BS","1 (OG)","PRJNA668472","RYGB","30","61"),
 c("Assal","1 (OG)","PRJNA395660","RYGB","25","71"),
 c("Afshar","1 (OG)","PRJEB14082","RYGB","31","50"),
 c("Ilhan","1 (OG)","PRJNA612042","RYGB","45","82"),
 c("Cholecystectomy","2","PRJNA1188648","Gallbladder removal","10","30"),
 c("IleocecalResection","2","PRJNA514452","Ileocecal resection (Crohn's)","54","189"),
 c("Ileostomy","2","PRJNA1480144 (IL)","Temporary ileostomy","44","65"),
 c("SDT","2","PRJNA1480144 (SD)","Stent-based fecal diversion","34","45"),
 c("RYGB_SaoPaulo","3b","PRJNA1364303","RYGB +/- exercise","31","87"),
 c("Lazaro_GB","3b","PRJNA859272 (GB)","RYGB","10","20"),
 c("Lazaro_SG","3a","PRJNA859272 (SG)","Sleeve gastrectomy","14","28"),
 c("Chaudhari_SG","3a","PRJNA635168","Sleeve gastrectomy","17","34")
)
cex_tbl <- 0.72
gap <- 0.015
ncol_tbl <- ncol(tbl)
colx <- numeric(ncol_tbl)
colx[1] <- 0.04
for (j in 2:ncol_tbl) {
  widest <- max(strwidth(tbl[,j-1], cex=cex_tbl, family="mono", font=2))
  colx[j] <- colx[j-1] + widest + gap
}
yy <- 0.85
for (r in 1:nrow(tbl)) {
  fw <- if (r==1) 2 else 1
  text(colx, yy, tbl[r,], adj=0, cex=cex_tbl, font=fw, family="mono")
  yy <- yy - 0.045
}
text(0.5, 0.10, "Phase 1 = original paper's 4 cohorts (\"OG\"). Phase 2 = unrelated GI surgeries, tests whether", cex=0.72, family="sans")
text(0.5, 0.07, "the signature generalizes to any gut surgery. Phase 3 = more bariatric cohorts as a positive", cex=0.72, family="sans")
text(0.5, 0.04, "control, split into RYGB-specific (3b) and non-RYGB bariatric / sleeve gastrectomy (3a).", cex=0.72, family="sans")

# ---- Page 3: methods in simple terms ----
plot.new()
text(0.5, 0.97, "How the Heatmaps and Scatterplots Were Built", cex=1.5, font=2, family="sans")
text(0.5, 0.94, "In plain terms, no statistics background required", cex=0.85, family="sans")

lines <- c(
"1. For every cohort, each subject's gut bacteria were measured before surgery and at one or",
"   more points after. Each distinct type of bacteria (a \"taxon\", e.g. the genus Veillonella)",
"   gets its own abundance reading per sample.",
"",
"2. For each taxon in each cohort, a statistical test asks: did this bacterium's abundance",
"   change significantly after surgery? That test returns two things -- a p-value (how",
"   confident we are it's a real change, not noise) and a direction (did it go up or down).",
"",
"3. Those two numbers are combined into one score per taxon: a large positive score means",
"   \"significantly increased after surgery\"; a large negative score means \"significantly",
"   decreased\"; near zero means \"no real change detected.\"",
"",
"4. SCATTERPLOT: to compare two cohorts, take every taxon they both measured and plot one",
"   cohort's scores on the x-axis against the other's on the y-axis. Each dot is one taxon.",
"   If the same bacteria rise and fall together in both cohorts, the dots trend along a",
"   diagonal line. If the two cohorts have nothing in common, the dots scatter randomly.",
"",
"5. The tightness of that diagonal trend is summarized as a single number, the Spearman",
"   correlation (r), printed above each scatterplot along with a p-value for r itself.",
"   r near 1 = strong agreement between cohorts; r near 0 = no agreement.",
"",
"6. HEATMAP: rather than looking at one scatterplot pair at a time, the heatmap puts every",
"   cohort/timepoint against every other in one grid. Each cell is just that pair's r and p",
"   value from step 5 -- the heatmap doesn't compute anything new, it's a compact way to see",
"   every pairwise comparison at once. Warmer colors (red) = stronger agreement, cooler",
"   colors (blue) = weaker. Grey cells compare a cohort's own two timepoints to each other",
"   (e.g. 1-month vs. 6-month within the same study) and are excluded from headline numbers,",
"   since the real test is agreement BETWEEN independent cohorts, not within one.",
"",
"The headline number reported for each phase (e.g. \"r = 0.40\") is just the average of all",
"the DIFFERENT-cohort cells in that phase's heatmap."
)
n_lines <- length(lines)
yy_top <- 0.89
yy_bot <- 0.04
step <- (yy_top - yy_bot) / (n_lines - 1)
yy <- yy_top
for (ln in lines) {
  text(0.05, yy, ln, adj=0, cex=0.74, family="mono")
  yy <- yy - step
}
dev.off()
cat("wrote", outFile, "\n")
