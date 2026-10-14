#Author: Rishi Gupta / Claude
#Date: 2026
#Description: Cross-group-only scatter plots for the lab-meeting report: every study-timepoint
#in groupA plotted against every study-timepoint in groupB (e.g. "new cohort vs. OG"), skipping
#within-group pairs entirely. Reuses the already-gathered MixedLinearModelResults files and
#compareStudies()/plotPairwiseStudies()/correlationBetweenStudies() from a completed
#CompareStudies16S module -- no new R statistics, just a filtered re-plot of the same pairwise
#comparison the pipeline already computed.
#Usage: Rscript crossGroupScatter.R <CompareStudies16S module dir> <groupA CSV> <groupB CSV> <title> <out_prefix>
#  <CompareStudies16S module dir> must contain input/ (gathered MixedLinearModelResults files)
#  and resources/functions.R, exactly as 06_analysis.sbatch leaves it after a run.

library(ggplot2)
library(gridExtra)
library(ggsignif)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
moduleDir <- args[1]
groupA <- trimws(strsplit(args[2], ",")[[1]])
groupB <- trimws(strsplit(args[3], ",")[[1]])
titleText <- args[4]
outPrefix <- args[5]

inputDir <- file.path(moduleDir, "input")
source(file.path(moduleDir, "resources", "functions.R"))
taxaLevel <- "Genus"

#Same STUDIES config as analysis/RScripts/compareStudies/compareStudies_16S.R -- duplicated
#rather than sourced so this script has no dependency on the pipeline's module layout beyond
#the one directory it's pointed at.
STUDIES <- list(
  BS     = list(timepoints=list(c("p1M","BS-1 month"), c("p6M","BS-6 months"))),
  Assal  = list(timepoints=list(c("p3M","Assal-3 months"), c("p1Y","Assal-1 year"), c("p2Y","Assal-2 years"))),
  Ilhan  = list(timepoints=list(c("p6M","Ilhan-6 months"), c("p1Y","Ilhan-1 year"))),
  Afshar = list(timepoints=list(c("p6M","Afshar-Post surgery"))),
  Cholecystectomy    = list(timepoints=list(c("p6M","Cholecystectomy-6 months"), c("p12M","Cholecystectomy-12 months"))),
  IleocecalResection = list(timepoints=list(c("p1","IleocecalResection-1 month"), c("p3","IleocecalResection-3 months"), c("p6","IleocecalResection-6 months"))),
  Ileostomy          = list(timepoints=list(c("pPost","Ileostomy-Post surgery"))),
  SDT                = list(timepoints=list(c("pPost","SDT-Post surgery"))),
  RYGB_SaoPaulo = list(timepoints=list(c("pT2","RYGB_SaoPaulo-T2"), c("pT3","RYGB_SaoPaulo-T3"))),
  Lazaro_GB     = list(timepoints=list(c("pPost","Lazaro_GB-Post surgery"))),
  Lazaro_SG     = list(timepoints=list(c("pPost","Lazaro_SG-Post surgery"))),
  Chaudhari_SG  = list(timepoints=list(c("pPost","Chaudhari_SG-Post surgery")))
)

studies<-character(0); studyNames<-character(0); studyOf<-character(0)
for (study in names(STUDIES)) {
  for (tp in STUDIES[[study]]$timepoints) {
    studies<-c(studies, paste0(study,"-",tp[1]))
    studyNames<-c(studyNames, tp[2])
    studyOf<-c(studyOf, study)
  }
}

r<-vector(); pval<-vector(); plotList<-list(); xlabs<-character(0); ylabs<-character(0)
index<-1
for (s in seq_along(studies)) {
  for (s1 in seq_along(studies)) {
    if (s1 <= s) next
    cross <- (studyOf[s] %in% groupA && studyOf[s1] %in% groupB) ||
             (studyOf[s] %in% groupB && studyOf[s1] %in% groupA)
    if (!cross) next

    df<-compareStudies(inputDir,taxaLevel,strsplit(studies[s],"-")[[1]][1],strsplit(studies[s1],"-")[[1]][1],
                        strsplit(studies[s],"-")[[1]][2],strsplit(studies[s1],"-")[[1]][2])
    rr<-correlationBetweenStudies(df)
    r[index]<-rr[[1]]; pval[index]<-rr[[2]]
    plotList[[index]]<-plotPairwiseStudies(df, paste0(studyNames[s]," vs. baseline"),
                                            paste0(studyNames[s1]," vs. baseline"), r[index], pval[index])
    xlabs[index]<-studyNames[s]; ylabs[index]<-studyNames[s1]
    index<-index+1
  }
}
pval<-p.adjust(pval,method="BH")

n <- length(plotList)
pdf(paste0(outPrefix,"_crossGroupScatter.pdf"), width = 10, height = 10)
theme_set(theme_classic(base_size = 9))
plot.new()
text(0.5, 0.75, titleText, cex = 1.4, font = 2)
text(0.5, 0.60, paste0("Group A: ", paste(groupA, collapse=", ")), cex = 0.9)
text(0.5, 0.53, paste0("Group B: ", paste(groupB, collapse=", ")), cex = 0.9)
text(0.5, 0.42, paste0(n, " cross-group pairs (within-group pairs excluded)"), cex = 0.9)
text(0.5, 0.35, "Each point: one taxon's signed -log10(p). r = Spearman correlation, p = BH-adjusted.", cex = 0.75)
if (n > 0) for (i in seq(1, n, by = 9)) {
  pageEnd <- min(i+8, n)
  do.call(gridExtra::grid.arrange, c(plotList[i:pageEnd], ncol=3, nrow=3))
}
dev.off()

if (n > 0) {
  write.table(data.frame(groupA_vs=xlabs, groupB_vs=ylabs, r=r, pval_BH=pval),
              paste0(outPrefix,"_crossGroupPairwise.txt"), sep="\t", row.names=FALSE, quote=FALSE)
}
cat("wrote", paste0(outPrefix,"_crossGroupScatter.pdf"), "with", n, "cross-group pairs\n")
