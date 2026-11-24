#Author: Rishi Gupta / Claude
#Date: 2026
#Description: Report cover pages for the WGS/shotgun metagenomics taxonomy track -- title,
#the full dataset table (all 6 identified cohorts, resolved and pending), and a
#plain-language explanation of how the species-level pairwise scatterplots are built.
#Mirrors introPages.R's structure and tone (same report style as
#results/phase3/report/pairwise_comparison_report.pdf), adapted for this track's
#different pipeline (kneaddata + Kraken2/Bracken, not DADA2 + 16S Kraken2) and its
#partial-completion state (2 of 6 cohorts still blocked on metadata).
#No statistics computed here; purely descriptive front matter for the combined report.
#Usage: Rscript introPagesMetagenomics.R <out.pdf>

args <- commandArgs(trailingOnly = TRUE)
outFile <- args[1]

pdf(outFile, width = 11, height = 8.5)

# ---- Page 1: title ----
plot.new()
text(0.5, 0.6, "Gut Microbiome Signature Correlation Analysis", cex=1.8, font=2)
text(0.5, 0.52, "WGS / Shotgun Metagenomics Taxonomy Track", cex=1.3)
text(0.5, 0.40, "Phase 4: extension of the 16S reproduction (see the main report) to whole-genome", cex=0.9)
text(0.5, 0.35, "shotgun sequencing cohorts, classified with Kraken2 + Bracken at species level", cex=0.9)

# ---- Page 2: all datasets table ----
plot.new()
par(family = "mono")
text(0.5, 0.97, "All WGS Datasets Identified For This Track", cex=1.5, font=2, family="sans")
text(0.5, 0.925, "4 of 6 cohorts fully resolved and analyzed below; 2 remain blocked on metadata", cex=0.85, family="sans")

tbl <- rbind(
 c("Cohort","Accession","Procedure","Subjects","Runs","Status"),
 c("Penney_RYGB","PRJNA473348","RYGB","10","25","Analyzed"),
 c("Penney_SG","PRJNA473348","Sleeve gastrectomy","14","31","Analyzed"),
 c("Oseberg_RYGB","PRJEB59861","RYGB","39","78","Analyzed"),
 c("Oseberg_SG","PRJEB59861","Sleeve gastrectomy","38","76","Analyzed"),
 c("COCKTAIL_GBP","PRJNA904651","Gastric bypass","unresolved","408*","PENDING (Arm unknown)"),
 c("Hayashi_LSG","PRJNA1329344","Lap. sleeve gastrectomy","19","38*","PENDING (Arm unknown)")
)
cex_tbl <- 0.72
gap <- 0.02
ncol_tbl <- ncol(tbl)
colx <- numeric(ncol_tbl)
colx[1] <- 0.04
for (j in 2:ncol_tbl) {
  widest <- max(strwidth(tbl[,j-1], cex=cex_tbl, family="mono", font=2))
  colx[j] <- colx[j-1] + widest + gap
}
yy <- 0.84
for (r in 1:nrow(tbl)) {
  fw <- if (r==1) 2 else 1
  text(colx, yy, tbl[r,], adj=0, cex=cex_tbl, font=fw, family="mono")
  yy <- yy - 0.05
}
text(0.5, 0.17, "210 samples analyzed across the 4 resolved cohorts (Penney_RYGB/SG, Oseberg_RYGB/SG).", cex=0.72, family="sans")
text(0.5, 0.14, "* COCKTAIL_GBP and Hayashi_LSG run counts are the FULL BioProject, including an explicitly", cex=0.72, family="sans")
text(0.5, 0.11, "excluded second arm (VLED diet / IGBP) not distinguishable from public metadata alone --", cex=0.72, family="sans")
text(0.5, 0.08, "see cohorts/<Study>/NEEDS_METADATA.md. Their per-cohort scripts stop() rather than run on", cex=0.72, family="sans")
text(0.5, 0.05, "placeholder metadata, so nothing below includes them until that's resolved.", cex=0.72, family="sans")

# ---- Page 3: methods in simple terms ----
plot.new()
text(0.5, 0.97, "How the Species-Level Scatterplots Were Built", cex=1.5, font=2, family="sans")
text(0.5, 0.94, "In plain terms, no statistics background required", cex=0.85, family="sans")

lines <- c(
"1. Unlike the 16S track (which sequences one marker gene), WGS/shotgun sequencing reads",
"   essentially all DNA in a sample. Each sample's reads are first cleaned (human DNA and",
"   low-quality bases removed, via kneaddata), then every remaining read is assigned to a",
"   taxon (a specific species, e.g. Veillonella parvula) by Kraken2, a reference-database",
"   classifier. Bracken then re-estimates each species' true abundance from those",
"   read-level calls -- more accurate than reading percentages straight off Kraken2's own",
"   report for WGS data.",
"",
"2. For each taxon in each cohort, a statistical test asks: did this species' abundance",
"   change significantly after surgery? That test returns a p-value (how confident we are",
"   it's a real change, not noise) and a direction (did it go up or down).",
"",
"3. Those two numbers are combined into one score per taxon: a large positive score means",
"   \"significantly increased after surgery\"; a large negative score means \"significantly",
"   decreased\"; near zero means \"no real change detected.\"",
"",
"4. SCATTERPLOT: to compare two cohorts (or two timepoints within one cohort), take every",
"   species they both measured and plot one cohort's scores on the x-axis against the",
"   other's on the y-axis. Each dot is one species. If the same bacteria rise and fall",
"   together in both, the dots trend along a diagonal line; if the two have nothing in",
"   common, the dots scatter randomly. The scatterplots on the following pages cover every",
"   one of the 15 possible pairs across the 4 analyzed cohorts' timepoints.",
"",
"5. The tightness of that diagonal trend is summarized as a single number, the Spearman",
"   correlation (r), printed above each scatterplot along with a p-value for r itself.",
"   r near 1 = strong agreement between cohorts; r near 0 = no agreement. Standout species",
"   (significant in the SAME direction in both cohorts) are labeled by name on each plot.",
"",
"6. HEATMAP: the same 15 pairwise r/p values, laid out as one grid (every cohort-timepoint",
"   against every other) instead of 15 separate plots -- a compact way to see every",
"   comparison at once. Grey cells compare a cohort's own two timepoints to each other and",
"   are excluded from the headline numbers, since the real test is agreement BETWEEN",
"   independent cohorts, not within one.",
"",
"7. TOP MOVERS: beyond pairwise agreement, which specific species show up as significant",
"   (BH-adjusted p<0.05) most consistently across the 6 cohort-timepoints analyzed here?",
"   See the final pages for that list and an important caveat about reading it.",
"",
"This is the same method the main 16S report uses (see its own methods page) -- the only",
"difference is the unit being compared is a WGS-classified species rather than a 16S genus."
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
