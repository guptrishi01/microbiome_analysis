#Author: Rishi Gupta / Claude
#Date: 2026
#Description: Re-plots the WGS track's pairwise scatterplots WITH taxon labels, matching
#the 16S track's plotPairwiseStudies() (functions.R) -- it labels taxa that are
#significantly concordant (BH p<0.05 in the SAME direction in both studies) via
#geom_text_repel, whereas plotPairwiseStudiesMetagenomics() (also functions.R, used by the
#real pipeline module compareStudies_Metagenomics.R) deliberately has no labeling at all.
#This is a REPORTING script, not a new analysis path (CLAUDE.md's "Reporting scripts"
#convention): it reads the already-gathered MixedLinearModelResults files
#(06_CompareStudiesMetagenomicsWGS's own input/ dir, populated by the real pipeline run)
#and re-plots a filtered/labeled view, computing nothing new -- same r/p/taxon data as
#Species_scatterPlots_MetagenomicsAllStudies.pdf, just with names on the standout points.
#
#Caps labels at the top N (by whichever study's significance is weaker, so a taxon has to
#be robustly significant in BOTH, not just one) per panel when the concordant set is too
#large to read -- confirmed necessary empirically: most of the 15 real pairs have 30-110
#concordant taxa (fine for geom_text_repel, same scale as 16S), but Oseberg_RYGB vs.
#Oseberg_SG has 860 (both large cohorts, high power to detect small differences as
#"significant" in both) -- unreadable unlabeled, let alone labeled. The 16S track's
#original df1 filter has no such cap since it never hits this scale.
#
#Usage: Rscript labeledScatterplotsMetagenomics.R <inputDir> <funcScript> <out.pdf> [cap=25]

suppressMessages({
  library(ggplot2)
  library(gridExtra)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
inputDir <- args[1]
funcScript <- args[2]
outFile <- args[3]
CAP <- if (length(args) >= 4) as.integer(args[4]) else 25L
source(funcScript)

taxaLevel <- "Species"

# Same 4 resolved cohorts, same timepoint labels as the real pipeline's
# compareStudies_Metagenomics.R STUDIES config -- kept in sync by hand, mirroring how
# 03c/04c already mirror 03b/04b by convention in this repo (no single source of truth).
STUDIES <- list(
  Penney_RYGB  = list(timepoints=list(c("pPOST","Penney RYGB-Post"), c("p1YrPOST","Penney RYGB-1 year"))),
  Penney_SG    = list(timepoints=list(c("pPOST","Penney SG-Post"), c("p1YrPOST","Penney SG-1 year"))),
  Oseberg_RYGB = list(timepoints=list(c("p12M","Oseberg RYGB-12 months"))),
  Oseberg_SG   = list(timepoints=list(c("p12M","Oseberg SG-12 months")))
)

studies<-character(0); studyNames<-character(0)
for (study in names(STUDIES)) {
  for (tp in STUDIES[[study]]$timepoints) {
    studies<-c(studies, paste0(study,"_Kraken2_Metagenomics-",tp[1]))
    studyNames<-c(studyNames, tp[2])
  }
}

# plotPairwiseStudiesMetagenomicsLabeled: same shape as functions.R's
# plotPairwiseStudiesMetagenomics(), with plotPairwiseStudies()'s (16S) concordant-taxa
# labeling added, plus the readability cap described above.
plotLabeled <- function(df, xlab, ylab, coeficient, p, cap) {
  df1 <- df[(df[,"pval1"]<log10(0.05) & df[,"pval2"]<log10(0.05)) |
            (df[,"pval1"]>-log10(0.05) & df[,"pval2"]>-log10(0.05)), ]
  if (nrow(df1) > cap) {
    strength <- pmin(abs(df1$pval1), abs(df1$pval2))
    df1 <- df1[order(-strength), ][seq_len(cap), ]
  }
  if (p == 0) p1 <- "< 2.2e-16" else p1 <- paste("=", format(p, digits=3))
  theme_set(theme_classic())
  ggplot(data=df, aes(x=pval1, y=pval2)) +
    geom_point(size=0.5, color=alpha("black",0.4)) +
    geom_hline(yintercept=0, linetype="dashed", color="red") +
    geom_vline(xintercept=0, linetype="dashed", color="red") +
    labs(x=xlab, y=ylab, title=paste0("Spearman coefficient = ", format(coeficient, digits=3),
                                       "\nAdjusted p ", p1,
                                       "  (", nrow(df1), " labeled)")) +
    geom_text_repel(data=df1, aes(x=pval1, y=pval2, label=bugName), segment.colour="red",
                     size=2, min.segment.length=0, segment.color="grey", segment.size=0.2,
                     max.overlaps=30)
}

# First pass: correlation + BH-adjusted p per pair (same order/method as the real pipeline)
r<-vector(); pval<-vector(); index<-1
for (s in seq_along(studies)) {
  if (s==length(studies)) next
  for (s1 in (s+1):length(studies)) {
    df<-compareStudies(inputDir,taxaLevel,strsplit(studies[s],"-")[[1]][1],strsplit(studies[s1],"-")[[1]][1],strsplit(studies[s],"-")[[1]][2],strsplit(studies[s1],"-")[[1]][2])
    r[index]<-correlationBetweenStudies(df)[[1]]
    pval[index]<-correlationBetweenStudies(df)[[2]]
    index<-index+1
  }
}
pval<-p.adjust(pval, method="BH")

plotList<-list(); count<-0
for (s in seq_along(studies)) {
  if (s==length(studies)) next
  for (s1 in (s+1):length(studies)) {
    count<-count+1
    df<-compareStudies(inputDir,taxaLevel,strsplit(studies[s],"-")[[1]][1],strsplit(studies[s1],"-")[[1]][1],strsplit(studies[s],"-")[[1]][2],strsplit(studies[s1],"-")[[1]][2])
    xlab<-paste0(studyNames[s]," vs. baseline")
    ylab<-paste0(studyNames[s1]," vs. baseline")
    plotList[[count]]<-plotLabeled(df, xlab, ylab, r[count], pval[count], CAP)
  }
}

pdf(outFile, width=10, height=10)
theme_set(theme_classic(base_size=9))
n<-length(plotList)
for (i in seq(1,n,by=9)) {
  pageEnd<-min(i+8,n)
  do.call(grid.arrange, c(plotList[i:pageEnd], ncol=3, nrow=3))
}
dev.off()
cat("wrote", outFile, "(", n, "pairs)\n")
