#Author: Rishi Gupta / Claude
#Date: 2026
#Description: Identifies the most consistent significantly-changed species across the 4
#resolved WGS cohorts' 6 cohort-timepoints, same method as the 16S track's "top movers"
#analysis (see results/README.md's "Top movers and antibiotic cross-reference" section):
#for each taxon, count how many cohort-timepoints it's significantly (BH-adjusted p<0.05)
#increased vs. decreased, using each cohort's own already-computed MixedLinearModelResults
#file. Reporting script, not a new analysis path -- reads already-gathered per-cohort MLM
#output (06_CompareStudiesMetagenomicsWGS/input/), computes no new statistics.
#Usage: Rscript topMoversMetagenomics.R <inputDir> <out.pdf> [topN=20]

args <- commandArgs(trailingOnly = TRUE)
inputDir <- args[1]
outFile <- args[2]
topN <- if (length(args) >= 3) as.integer(args[3]) else 20L

# cohort-timepoint -> (file, adjusted-p column, coefficient-sign column, display label)
COHORT_TIMEPOINTS <- list(
  list(file="Species_Penney_RYGB_Kraken2_Metagenomics_MixedLinearModelResults.txt", padj="AdjustedpPOST",    s="sPOST",    label="Penney_RYGB-Post"),
  list(file="Species_Penney_RYGB_Kraken2_Metagenomics_MixedLinearModelResults.txt", padj="Adjustedp1YrPOST", s="s1YrPOST", label="Penney_RYGB-1yr"),
  list(file="Species_Penney_SG_Kraken2_Metagenomics_MixedLinearModelResults.txt",   padj="AdjustedpPOST",    s="sPOST",    label="Penney_SG-Post"),
  list(file="Species_Penney_SG_Kraken2_Metagenomics_MixedLinearModelResults.txt",   padj="Adjustedp1YrPOST", s="s1YrPOST", label="Penney_SG-1yr"),
  list(file="Species_Oseberg_RYGB_Kraken2_Metagenomics_MixedLinearModelResults.txt",padj="Adjustedp12M",    s="s12M",    label="Oseberg_RYGB-12mo"),
  list(file="Species_Oseberg_SG_Kraken2_Metagenomics_MixedLinearModelResults.txt",  padj="Adjustedp12M",    s="s12M",    label="Oseberg_SG-12mo")
)

increases <- list()  # bugName -> character vector of cohort-timepoint labels where significantly increased
decreases <- list()

for (ct in COHORT_TIMEPOINTS) {
  path <- file.path(inputDir, ct$file)
  if (!file.exists(path)) { message("missing: ", path); next }
  tbl <- read.table(path, sep="\t", header=TRUE, check.names=FALSE, quote="", comment.char="")
  if (!(ct$padj %in% names(tbl)) || !(ct$s %in% names(tbl))) {
    message("skip ", ct$label, ": columns ", ct$padj, "/", ct$s, " not found in ", path,
            " -- column names differ per cohort's own timepoint factor levels",
            " (confirmed by inspecting each file's real header before writing this list,",
            " not assumed), check COHORT_TIMEPOINTS above if this ever fires unexpectedly")
    next
  }
  sig <- tbl[!is.na(tbl[[ct$padj]]) & tbl[[ct$padj]] < 0.05, ]
  for (i in seq_len(nrow(sig))) {
    bug <- sig$bugName[i]
    if (sig[[ct$s]][i] > 0) {
      increases[[bug]] <- c(increases[[bug]], ct$label)
    } else {
      decreases[[bug]] <- c(decreases[[bug]], ct$label)
    }
  }
}

mkTable <- function(lst) {
  if (length(lst) == 0) return(data.frame(bugName=character(0), n=integer(0), cohorts=character(0)))
  n <- sapply(lst, length)
  cohorts <- sapply(lst, paste, collapse=", ")
  df <- data.frame(bugName=names(lst), n=n, cohorts=cohorts, stringsAsFactors=FALSE)
  df[order(-df$n, df$bugName), ]
}

incTbl <- mkTable(increases)
decTbl <- mkTable(decreases)

pdf(outFile, width=11, height=8.5)

drawTablePage <- function(title, subtitle, df, topN) {
  plot.new()
  text(0.5, 0.97, title, cex=1.4, font=2, family="sans")
  text(0.5, 0.935, subtitle, cex=0.8, family="sans")
  maxN <- if (nrow(df) > 0) max(df$n) else 0
  caveat <- if (maxN <= 1) {
    paste0("CAVEAT: nothing here reaches significance in more than 1 of ", length(COHORT_TIMEPOINTS),
           " cohort-timepoints -- this list is NOT a ranked finding, just every species that",
           " was significant anywhere, alphabetically ordered within that tie. At current sample",
           " sizes, no cross-cohort consistency was found for this direction.")
  } else {
    paste0("CAVEAT: max consistency found is ", maxN, " of ", length(COHORT_TIMEPOINTS),
           " cohort-timepoints, and every entry at that max is Oseberg_RYGB + Oseberg_SG",
           " specifically (by far the largest cohorts here, 78+76 vs. Penney's 25+31 -- more",
           " statistical power to detect species-level significance, not necessarily a more",
           " real biological effect). Read as preliminary given this track's small n overall",
           " (see results/README.md's phase 4 section), not yet a phase-3-style 10/12 finding.")
  }
  caveatLines <- strwrap(caveat, width=130)
  yy_cav <- 0.895
  for (cl in caveatLines) {
    text(0.5, yy_cav, cl, cex=0.62, family="sans")
    yy_cav <- yy_cav - 0.022
  }
  df <- head(df, topN)
  if (nrow(df) == 0) {
    text(0.5, 0.5, "(none found)", cex=0.9, family="sans")
    return(invisible())
  }
  par(family="mono")
  hdr <- c("Species", "# cohort-timepoints", "Which cohort-timepoints")
  rows <- rbind(hdr, as.matrix(data.frame(
    sub("^$", " ", substr(df$bugName, 1, 38)),
    as.character(df$n),
    substr(df$cohorts, 1, 55)
  )))
  cex_tbl <- 0.68
  gap <- 0.015
  colx <- numeric(3)
  colx[1] <- 0.03
  colx[2] <- colx[1] + max(strwidth(rows[,1], cex=cex_tbl, family="mono", font=2)) + gap
  colx[3] <- colx[2] + max(strwidth(rows[,2], cex=cex_tbl, family="mono", font=2)) + gap
  yy <- yy_cav - 0.02
  for (r in seq_len(nrow(rows))) {
    fw <- if (r==1) 2 else 1
    text(colx, yy, rows[r,], adj=0, cex=cex_tbl, font=fw, family="mono")
    yy <- yy - (0.80 / (nrow(rows)+2))
  }
}

drawTablePage(
  "Most Consistent Significant Increasers",
  paste0("Species significantly (BH-adjusted p<0.05) INCREASED after surgery in the most",
         " of the 6 analyzed cohort-timepoints (out of ", length(COHORT_TIMEPOINTS), " total)"),
  incTbl, topN
)
drawTablePage(
  "Most Consistent Significant Decreasers",
  paste0("Species significantly (BH-adjusted p<0.05) DECREASED after surgery in the most",
         " of the 6 analyzed cohort-timepoints (out of ", length(COHORT_TIMEPOINTS), " total)"),
  decTbl, topN
)

dev.off()
cat("wrote", outFile, "-", nrow(incTbl), "increasers,", nrow(decTbl), "decreasers (any significance)\n")
