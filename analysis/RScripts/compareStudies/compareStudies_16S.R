#Author: phase-2 extension (Rishi Gupta / Claude), rewritten from Farnaz Fouladi's compareStudies_16S.R
#Date: 2026
#Description: p-value versus p-value plots for 16S datasets
#
#The original hardcodes exactly 4 studies (BS/Assal/Ilhan/Afshar) as separate file-copy blocks,
#fixed studies/studyNames vectors, and a scatter-plot pagination loop assuming exactly 28 pairs
#(hardcoded page breaks at i=1,10,19). This is a generic rewrite over a STUDIES config so a
#phase-2 cohort is one new list entry. See CLAUDE.md's "Adding a new cohort" section and
#"Never assume" section — verify this against real phase-1 output (restricted via
#PHASE2_STUDIES) before trusting it for phase 2.
#
#compareStudies() itself (functions.R) already takes plain study/timepoint-suffix strings and
#needed no changes — only this script's per-study config and plotting scaffolding did.

library(ggplot2)
library(gridExtra)
library(ggsignif)
library(ggrepel)
library(stringr)

pipeRoot <- dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
outDir <- file.path(moduleDir, "output")
inputDir = file.path(moduleDir,"input")
dir.create(inputDir, showWarnings = FALSE)
message("Input tables will collected into an input folder: ", inputDir)

taxaLevel <- "Genus"
funcScript <- file.path(file.path(moduleDir, "resources"),"functions.R")
source(funcScript)

#study -> (DADA module name fragment, list of (column suffix, display label) pairs)
#Display labels for BS/Assal/Ilhan/Afshar are kept identical to the original script's
#studyNames text, since that text ends up in output plots/tables for phase 1.
STUDIES <- list(
  BS     = list(module="BSDADA",     timepoints=list(c("p1M","BS-1 month"), c("p6M","BS-6 months"))),
  Assal  = list(module="AssalDADA",  timepoints=list(c("p3M","Assal-3 months"), c("p1Y","Assal-1 year"), c("p2Y","Assal-2 years"))),
  Ilhan  = list(module="IlhanDADA",  timepoints=list(c("p6M","Ilhan-6 months"), c("p1Y","Ilhan-1 year"))),
  Afshar = list(module="AfsharDADA", timepoints=list(c("p6M","Afshar-Post surgery"))),
  Cholecystectomy    = list(module="CholecystectomyDADA",    timepoints=list(c("p6M","Cholecystectomy-6 months"), c("p12M","Cholecystectomy-12 months"))),
  IleocecalResection = list(module="IleocecalResectionDADA", timepoints=list(c("p1","IleocecalResection-1 month"), c("p3","IleocecalResection-3 months"), c("p6","IleocecalResection-6 months"))),
  Ileostomy          = list(module="IleostomyDADA",          timepoints=list(c("pPost","Ileostomy-Post surgery"))),
  SDT                = list(module="SDTDADA",                timepoints=list(c("pPost","SDT-Post surgery")))
)
active <- Sys.getenv("PHASE2_STUDIES", unset = "")
if (nzchar(active)) {
  keep <- trimws(strsplit(active, ",")[[1]])
  STUDIES <- STUDIES[keep]
}

#Gather each active study's MixedLinearModelResults file into inputDir (was 4 hand-written blocks)
for (study in names(STUDIES)) {
  cfg <- STUDIES[[study]]
  fileName <- paste0(taxaLevel,"_",study,"_MixedLinearModelResults.txt")
  moduleDirPath <- file.path(dir(pipeRoot, cfg$module, full.names = TRUE))
  srcFile <- file.path(file.path(file.path(moduleDirPath, "output"),"MixedLinearModels"),fileName)
  message("Grabbing input file: ", fileName, " from module: ", moduleDirPath)
  file.copy(from = srcFile, to = inputDir)
}

#Flatten STUDIES into the same studies/studyNames vectors the original built by hand
studies<-character(0); studyNames<-character(0)
for (study in names(STUDIES)) {
  for (tp in STUDIES[[study]]$timepoints) {
    studies<-c(studies, paste0(study,"-",tp[1]))
    studyNames<-c(studyNames, tp[2])
  }
}

r<-vector()
pval<-vector()
plotList<-list()
studyPairs<-vector()
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
      plot<-plotPairwiseStudies(df,xlab,ylab,r[count],pval[count])
      plotList[[count]]<-plot
      studyPairs[count]<-paste0(studies[s],"_",studies[s1])
      if(strsplit(studies[s],"-")[[1]][1] == strsplit(studies[s1],"-")[[1]][1])
        compariosn[[count]]<-"Same study"
      else
        compariosn[[count]]<-"Different study"

      if(strsplit(studies[s],"-")[[1]][2] == strsplit(studies[s1],"-")[[1]][2])
        compariosnTime[[count]]<-"Same timepoint"
      else
        compariosnTime[[count]]<-"Different timepoint"

      # add the specific images used in the main figure in a favs list (unchanged: these two
      # pairs are specific to the original paper's Figure 1c/1d and still exist regardless of
      # how many more studies are now being compared)
      if (xlab == paste0("BS-1 month"," vs. baseline") & ylab == paste0("Assal-3 months"," vs. baseline")){
        favs[1] <- count
      }
      if (xlab == paste0("Assal-3 months"," vs. baseline") & ylab == paste0("Assal-1 year"," vs. baseline")){
        favs[2] <- count
      }
    }
  }
}

df<-data.frame(studyPairs,compariosn,compariosnTime,pval,r)
pairOutFile <- file.path(outDir, paste0(taxaLevel,"_PairwiseComparison.txt"))
message("Saving table of pairwise comparisons: ", pairOutFile)
write.table(df, pairOutFile, sep="\t", row.names = FALSE)

#Wilcoxon test for comparing coefficients
p1<-wilcox.test(df$r[df$compariosn=="Different study" & df$compariosnTime=="Different timepoint"],
                df$r[df$compariosn=="Different study" & df$compariosnTime=="Same timepoint"])

p2<-wilcox.test(df$r[df$compariosn=="Different study" & df$compariosnTime=="Different timepoint"],
                df$r[df$compariosn=="Same study" & df$compariosnTime=="Different timepoint"])

p3<-wilcox.test(df$r[df$compariosn=="Different study" & df$compariosnTime=="Same timepoint"],
                df$r[df$compariosn=="Same study" & df$compariosnTime=="Different timepoint"])

p4<-wilcox.test(df$r[df$compariosn=="Different study"],
                df$r[df$compariosn=="Same study" & df$compariosnTime=="Different timepoint"])

outFileWilRes <- file.path(outDir, paste0(taxaLevel, "_wilcoxResults.txt"))
message("Creating output file: ", outFileWilRes)
adjusted<-p.adjust(c(p1$p.value,p2$p.value,p3$p.value,p4$p.value),method = "BH")
capture.output(p1, p2, p3, p4, adjusted, file = outFileWilRes)

plot1<-ggplot(data=df,aes(x=compariosn,y=r))+
  geom_boxplot(aes(color=compariosnTime),position = position_dodge(0.6), width = 0.5, size = 0.4,outlier.shape = NA)+
  geom_jitter(aes(color=compariosnTime),size=0.8,position = position_dodge(0.6))+labs(x="",y="Spearman Coefficient",color="")+
  geom_signif(y_position = c(0.75,0.8),xmin = c(0.9,1.1),xmax=c(2,2),annotations = c("*","*"),tip_length=0.03,vjust = 0.5,textsize =4)

#Generic pagination: 9 plots per page (3x3 grid), any remainder on its own final page.
#Replaces the original's hardcoded `for (i in c(1,10,19))` + separate trailing grid.arrange,
#which assumed exactly 28 pairs (3 full pages + 1 leftover). Produces identical pages for the
#phase-1-only case (28 pairs -> pages at i=1,10,19 plus a 1-plot page), verified below.
scatPlotOut <- file.path(outDir, paste0(taxaLevel,"_scatterPlots.pdf"))
message("Saving scatter plots to file: ", scatPlotOut)
pdf(scatPlotOut, width = 10, height = 10)
theme_set(theme_classic(base_size = 9))
n <- length(plotList)
for (i in seq(1, n, by = 9)) {
  pageEnd <- min(i+8, n)
  do.call(grid.arrange, c(plotList[i:pageEnd], ncol=3, nrow=3))
}
dev.off()

fig1Pics <- file.path(outDir, paste0(taxaLevel,"_mainFigure1cd.pdf"))
message("Saving specific plots to individual file to compare to figure 1, see file: ", fig1Pics)
pdf(fig1Pics, width = 5, height = 9)
theme_set(theme_classic(base_size = 9))
grid.arrange(plotList[[favs[1]]], plotList[[favs[2]]],ncol=1,nrow=2)
dev.off()

coefOut <- file.path(outDir, paste0(taxaLevel,"_coefficientsFromScatterPlots.pdf"))
message("Saving coefficients from scatter plots: ", coefOut)
pdf(coefOut, width = 5,height = 5)
theme_set(theme_classic(base_size = 14))
print(plot1)
dev.off()

#plot coefficients for each study
df$study1<-sapply(as.character(df$studyPairs),function(x){strsplit(x,"-")[[1]][1]})
df$study2<-sapply(sapply(as.character(df$studyPairs),function(x){strsplit(x,"-")[[1]][2]}),
                  function(x){strsplit(x,"_")[[1]][2]})

#Generic per-study coefficient extraction: every pair touching study X (as either study1 or
#study2), with the "other" study recorded for the x-axis facet. Replaces the original's
#hand-written r.BS/r.Assal/r.Ilhan/r.Afshar cumulative blocks — verified equivalent: for a
#study X, "all pairs where study1==X or study2==X" is exactly what those blocks compute (the
#pairwise loop always orders study1 before study2 in STUDIES iteration order, so a pair is
#never double-counted).
studyKeys <- names(STUDIES)
rVals <- numeric(0); study1Out <- character(0); study2Out <- character(0)
for (X in studyKeys) {
  touching <- df[df$study1==X | df$study2==X, ]
  if (nrow(touching)==0) next
  other <- ifelse(touching$study1==X, touching$study2, touching$study1)
  rVals <- c(rVals, touching$r)
  study1Out <- c(study1Out, rep(X, nrow(touching)))
  study2Out <- c(study2Out, other)
}
df2<-data.frame(r=rVals,study1=study1Out,study2=study2Out)

theme_set(theme_gray(base_size = 14))
plot<-ggplot(data=df2,aes(x=factor(study2),y=r))+geom_boxplot()+geom_jitter(position = position_jitter(width=0.1))+facet_wrap(vars(study1))+
  labs(x="Studies",y="Spearman Coefficient")
boxPlotFile <- file.path(outDir, paste0(taxaLevel,"_BoxPlotCorrelations.pdf"))
message("Saving box plots to file: ", boxPlotFile)
pdf(boxPlotFile, width = 5,height = 5)
print(plot)
dev.off()

message("")
message("All done!")
message("")

sessionInfo()
