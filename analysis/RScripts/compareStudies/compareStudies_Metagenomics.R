#Author: WGS-taxonomy extension (Rishi Gupta / Claude), rewritten from Farnaz Fouladi's
#        compareStudies_Metagenomics.R
#Date: 2026
#Description: p-value versus p-value plots for WGS/shotgun metagenomic datasets
#
#The original hardcodes exactly 2 studies (BS/Palleja) as separate file-copy blocks, a
#fixed studies/studyNames vector of 5 study-timepoint entries, and a plotting loop
#assuming exactly 10 pairs. This is a generic rewrite over a STUDIES config, same pattern
#already used for compareStudies_16S.R (see that file's own header comment and CLAUDE.md's
#"Adding a new cohort" section) — a new WGS cohort is one new list entry here, not a
#rewrite of this script.
#
#compareStudies(), correlationBetweenStudies(), and plotPairwiseStudiesMetagenomics()
#(functions.R) already take plain study/timepoint-suffix strings and needed NO changes —
#same as when compareStudies_16S.R was generalized. Only this script's per-study config
#and plotting scaffolding are new.
#
#Unlike compareStudies_16S.R, COCKTAIL_GBP and Hayashi_LSG are NOT yet runnable end to
#end — see each cohort's cohorts/<Study>/NEEDS_METADATA.md. This script defensively
#SKIPS (with a clear message) any study whose module directory doesn't exist in the
#pipeline yet, rather than failing hard, since "not all configured WGS cohorts have
#resolved metadata yet" is expected to be true for a while, not a bug to crash on.
#
#Oseberg_RYGB/Oseberg_SG are fully resolved and wired (2026-11-14) — the original unsplit
#"Oseberg" entry was replaced with these two once Arm was recovered from the paper's own
#Zenodo source-data deposit; see cohorts/Oseberg_RYGB/ and cohorts/Oseberg_SG/ and
#Oseberg_RYGB_kraken2_Metagenomics.R / Oseberg_SG_kraken2_Metagenomics.R.
#
#Penney_RYGB/Penney_SG are PARTIALLY resolved and wired (2026-11-14) — the original
#unsplit "Penney" entry (n=120 runs) was replaced with these two once a confirmed-but-
#incomplete Arm mapping was recovered from the paper's own Mendeley deposit (only 25+31=56
#of the 120 runs have a known Arm); see cohorts/Penney_RYGB/ and cohorts/Penney_SG/
#NEEDS_METADATA.md for exactly what's missing and why.

library(ggplot2)
library(gridExtra)
library(ggsignif)
library(ggrepel)
library(stringr)

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
dir.create(paste0(moduleDir,"/input/"), showWarnings = FALSE)
inputDir = file.path(moduleDir,"input/")
outDir = file.path(moduleDir, "output/")

taxaLevel <- "Species"
funcScript <- paste0(moduleDir,"/resources/functions.R")
source(funcScript)

#study -> (MLM module name fragment, list of (column suffix, display label) pairs)
#BS/Palleja display labels kept identical to the original script's studyNames text, since
#that text ends up in output plots/tables that predate this generalization.
STUDIES <- list(
  BS      = list(module="BSMetagenomics",      timepoints=list(c("p1M","BS-1 month"), c("p6M","BS-6 months"), c("p1Y","BS-1 year"))),
  Palleja = list(module="PallejaMetagenomics", timepoints=list(c("p3M","Palleja-3 months"), c("p1Y","Palleja-1 year"))),
  # New WGS cohorts, identified this session — see cohorts/<Study>/NEEDS_METADATA.md for
  # each one's metadata status. Penney_RYGB/SG and Oseberg_RYGB/SG are fully resolved and
  # have run successfully (see results/phase4/); COCKTAIL_GBP and Hayashi_LSG remain
  # blocked on unresolved Arm metadata and correctly SKIP below (module directory exists,
  # but its per-cohort script stop()s before producing output — see each
  # kraken2_Metagenomics/*.R). All 6 are listed unconditionally (not gated behind
  # PHASE4_WGS_STUDIES) since the gather loop below already skips whichever aren't ready.
  Penney_RYGB      = list(module="Penney_RYGBMetagenomics",      timepoints=list(c("pPOST","Penney RYGB-Post"), c("p1YrPOST","Penney RYGB-1 year"))),
  Penney_SG        = list(module="Penney_SGMetagenomics",        timepoints=list(c("pPOST","Penney SG-Post"), c("p1YrPOST","Penney SG-1 year"))),
  Oseberg_RYGB     = list(module="Oseberg_RYGBMetagenomics",     timepoints=list(c("p12M","Oseberg RYGB-12 months"))),
  Oseberg_SG       = list(module="Oseberg_SGMetagenomics",       timepoints=list(c("p12M","Oseberg SG-12 months"))),
  COCKTAIL_GBP     = list(module="COCKTAIL_GBPMetagenomics",     timepoints=list(c("pWeek0","COCKTAIL-Week 0 (pre-surgery)"), c("pWeek2","COCKTAIL-Week 2"), c("pWeek6","COCKTAIL-Week 6"), c("pMonth4","COCKTAIL-Month 4"), c("pYear1","COCKTAIL-Year 1"), c("pYear2","COCKTAIL-Year 2"))),
  Hayashi_LSG      = list(module="Hayashi_LSGMetagenomics",      timepoints=list(c("pPost","Hayashi LSG-6 months")))
)
active <- Sys.getenv("PHASE4_WGS_STUDIES", unset = "")
if (nzchar(active)) {
  keep <- trimws(strsplit(active, ",")[[1]])
  STUDIES <- STUDIES[intersect(names(STUDIES), keep)]
}

#Gather each active study's MixedLinearModelResults file into inputDir, skipping (with a
#message, not an error) any whose module hasn't actually produced output yet.
present <- character(0)
for (study in names(STUDIES)) {
  cfg <- STUDIES[[study]]
  moduleMatch <- dir(pipeRoot, cfg$module, full.names = TRUE)
  if (length(moduleMatch) == 0) {
    message("SKIP ", study, ": no pipeline module matching '", cfg$module,
            "' found under ", pipeRoot, " — see cohorts/", study, "/NEEDS_METADATA.md")
    next
  }
  fileName <- paste0(taxaLevel,"_",study,"_Kraken2_Metagenomics_MixedLinearModelResults.txt")
  srcFile <- file.path(moduleMatch[1], "output", "MixedLinearModels", "Kraken2", fileName)
  if (!file.exists(srcFile)) {
    message("SKIP ", study, ": module directory exists but ", fileName, " was not found in its output")
    next
  }
  message("Grabbing input file: ", fileName, " from module: ", moduleMatch[1])
  file.copy(from = srcFile, to = inputDir, overwrite = TRUE)
  present <- c(present, study)
}
STUDIES <- STUDIES[present]
if (length(STUDIES) < 2) {
  stop("Fewer than 2 WGS cohorts have usable output (", length(STUDIES),
       " found: ", paste(present, collapse=", "), ") — nothing to cross-compare yet. ",
       "This is expected until the metadata gaps in cohorts/*/NEEDS_METADATA.md are resolved.")
}

#Flatten STUDIES into the same studies/studyNames vectors the original built by hand
studies<-character(0); studyNames<-character(0)
for (study in names(STUDIES)) {
  for (tp in STUDIES[[study]]$timepoints) {
    # bugName file suffix convention matches each *_kraken2_Metagenomics.R script's own
    # naming (e.g. "BS" not "BS_Kraken2_Metagenomics" — compareStudies() builds the
    # filename as "<taxa>_<study>_MixedLinearModelResults.txt", but our per-cohort
    # scripts write "<taxa>_<Study>_Kraken2_Metagenomics_MixedLinearModelResults.txt" —
    # see the file.copy block above, which renames nothing, so this must match exactly
    # what was grabbed. Kept as "<Study>_Kraken2_Metagenomics" to match verbatim.
    studies<-c(studies, paste0(study,"_Kraken2_Metagenomics-",tp[1]))
    studyNames<-c(studyNames, tp[2])
  }
}

r<-vector()
pval<-vector()
plotList<-list()
studyPairs<-vector()
study1Key<-vector()
study2Key<-vector()
compariosn<-vector()
compariosnTime<-vector()
index<-1

favs<-vector(length=2)

for (s in 1:length(studies)){

  if (s!=length(studies)){

    otherStudies<-c((s+1):length(studies))

    for (s1 in otherStudies){

      df<-compareStudies(inputDir,taxaLevel,strsplit(studies[s],"-")[[1]][1],strsplit(studies[s1],"-")[[1]][1],strsplit(studies[s],"-")[[1]][2],strsplit(studies[s1],"-")[[1]][2])

      r[index]<-correlationBetweenStudies(df)[[1]]
      pval[index]<-correlationBetweenStudies(df)[[2]]
      index<-index+1

    }
  }
}

pval<-p.adjust(pval,method = "BH")
count=0

for (s in 1:length(studies)){

  if (s!=length(studies)){
    otherStudies<-c((s+1):length(studies))

    for (s1 in otherStudies){

      count<-count+1
      df<-compareStudies(inputDir,taxaLevel,strsplit(studies[s],"-")[[1]][1],strsplit(studies[s1],"-")[[1]][1],strsplit(studies[s],"-")[[1]][2],strsplit(studies[s1],"-")[[1]][2])
      xlab=paste0(studyNames[s]," vs. baseline")
      ylab=paste0(studyNames[s1]," vs. baseline")
      plot<-plotPairwiseStudiesMetagenomics(df,xlab,ylab,r[count],pval[count])
      plotList[[count]]<-plot
      studyPairs[count]<-paste0(studies[s],"_",studies[s1])
      # Captured directly here (not re-parsed from studyPairs later) because several of
      # this script's study keys contain underscores themselves (e.g. "BS_Kraken2_
      # Metagenomics", "Penney_RYGB_Kraken2_Metagenomics") — the original compareStudies_
      # 16S.R's later strsplit(x,"_")[[1]][2] trick only works when the study-name half of
      # the pair has exactly one underscore-free token, which isn't true here.
      study1Key[count]<-strsplit(studies[s],"-")[[1]][1]
      study2Key[count]<-strsplit(studies[s1],"-")[[1]][1]
      if(strsplit(studies[s],"-")[[1]][1] == strsplit(studies[s1],"-")[[1]][1])
        compariosn[[count]]<-"Same study"
      else
        compariosn[[count]]<-"Different study"

      if(strsplit(studies[s],"-")[[1]][2] == strsplit(studies[s1],"-")[[1]][2])
        compariosnTime[[count]]<-"Same timepoint"
      else
        compariosnTime[[count]]<-"Different timepoint"

      # the original paper's Figure 3A/3B pair — only meaningful while BS/Palleja are
      # both active; harmless no-op (favs stays 0/0, nothing drawn for it) otherwise.
      if (xlab == paste0("BS-1 month"," vs. baseline") & ylab == paste0("Palleja-3 months"," vs. baseline")){
        favs[1] <- count
      }
      if (xlab == paste0("BS-1 month"," vs. baseline") & ylab == paste0("BS-6 months"," vs. baseline")){
        favs[2] <- count
      }
    }
  }

}

df<-data.frame(studyPairs,study1Key,study2Key,compariosn,compariosnTime,pval,r)
write.table(df, paste0(outDir,taxaLevel,"_PairwiseComaprison_MetagenomicsData.txt"),sep="\t",row.names = FALSE)

#Wilcoxon test for comparing coefficients
p1<-wilcox.test(df$r[df$compariosn=="Different study" & df$compariosnTime=="Different timepoint"],
                df$r[df$compariosn=="Different study" & df$compariosnTime=="Same timepoint"])

p2<-wilcox.test(df$r[df$compariosn=="Different study" & df$compariosnTime=="Different timepoint"],
                df$r[df$compariosn=="Same study" & df$compariosnTime=="Different timepoint"])

p3<-wilcox.test(df$r[df$compariosn=="Different study" & df$compariosnTime=="Same timepoint"],
                df$r[df$compariosn=="Same study" & df$compariosnTime=="Different timepoint"])

p4<-wilcox.test(df$r[df$compariosn=="Different study"],
                df$r[df$compariosn=="Same study" & df$compariosnTime=="Different timepoint"])

adjustedp<-p.adjust(c(p1$p.value,p2$p.value,p3$p.value,p4$p.value),method = "BH")
capture.output(p1, p2, p3, p4, adjustedp, file = paste0(outDir,taxaLevel,"_wilcoxResults.txt"))

plot1<-ggplot(data=df,aes(x=compariosn,y=r))+
  geom_boxplot(aes(color=compariosnTime),position = position_dodge(0.6), width = 0.5, size = 0.4,outlier.shape = NA)+
  geom_jitter(aes(color=compariosnTime),size=0.8,position = position_dodge(0.6))+labs(x="",y="Spearman Coefficient",color="")+
  lims(y=c(-1,1))

#Generic pagination: 9 plots per page (3x3 grid), any remainder on its own final page —
#same approach already used in compareStudies_16S.R, replacing the original's two
#hardcoded grid.arrange() calls assuming exactly 10 pairs (9 + 1 leftover).
scatPlotOut <- paste0(outDir,taxaLevel,"_scatterPlots_MetagenomicsAllStudies.pdf")
pdf(scatPlotOut, width = 10, height = 10)
theme_set(theme_classic(base_size = 9))
n <- length(plotList)
for (i in seq(1, n, by = 9)) {
  pageEnd <- min(i+8, n)
  do.call(grid.arrange, c(plotList[i:pageEnd], ncol=3, nrow=3))
}
dev.off()

if (favs[1] > 0 && favs[2] > 0) {
  pdf(paste0(outDir,taxaLevel,"_mainFigure3ab.pdf"), width = 7, height = 3)
  theme_set(theme_classic(base_size = 9))
  grid.arrange(plotList[[favs[1]]], plotList[[favs[2]]],ncol=2,nrow=1)
  dev.off()
}

pdf(paste0(outDir,taxaLevel,"_coefficientsFromScatterPlots_MetagenomicsAllStudies.pdf"), width = 5,height = 5)
theme_set(theme_classic(base_size = 14))
print(plot1)
dev.off()

#plot coefficients for each study — same generic per-study extraction already used in
#compareStudies_16S.R, but keyed off study1Key/study2Key (captured directly during the
#pairwise loop above) rather than re-parsed from studyPairs: unlike the 16S script's
#single-token study names (BS, Assal, ...), these WGS module keys contain underscores
#themselves (BS_Kraken2_Metagenomics, Penney_RYGB_Kraken2_Metagenomics, ...), so a
#strsplit(x,"_")[[1]][2] trick would silently grab the wrong substring.
studyKeys <- names(STUDIES)
rVals <- numeric(0); study1Out <- character(0); study2Out <- character(0)
for (X in studyKeys) {
  xKey <- paste0(X,"_Kraken2_Metagenomics")
  touching <- df[df$study1Key==xKey | df$study2Key==xKey, ]
  if (nrow(touching)==0) next
  other <- ifelse(touching$study1Key==xKey, touching$study2Key, touching$study1Key)
  other <- sub("_Kraken2_Metagenomics$", "", other)   # cosmetic: bare study name for the axis label
  rVals <- c(rVals, touching$r)
  study1Out <- c(study1Out, rep(X, nrow(touching)))
  study2Out <- c(study2Out, other)
}
df2<-data.frame(r=rVals,study1=study1Out,study2=study2Out)

theme_set(theme_gray(base_size = 14))
plot<-ggplot(data=df2,aes(x=factor(study2),y=r))+geom_boxplot()+geom_jitter(position = position_jitter(width=0.1))+facet_wrap(vars(study1))+
  labs(x="Studies",y="Spearman Coefficient")
pdf(paste0(outDir,taxaLevel,"_BoxPlotCorrelations.pdf"), width = 5,height = 5)
print(plot)
dev.off()

message("")
message("WGS cohorts actually compared: ", paste(present, collapse=", "))
message("All done!")
message("")

sessionInfo()
