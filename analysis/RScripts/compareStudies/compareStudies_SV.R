#Author: phase-2 extension (Rishi Gupta / Claude), rewritten from Farnaz Fouladi's compareStudies_SV.R
#Date: 2026
#Description: p-value versus p-value plots for 16S datasets (Sequence variants)
#
#IMPORTANT — do not add phase-2 cohorts to STUDIES below without verifying primer/region
#compatibility first. The original script ONLY ever compares BS and Assal at the exact-sequence
#(SV/ASV) level, deliberately excluding Ilhan and Afshar even though all 4 studies have SV-level
#MixedLinearModelResults files. Why: BS and Assal both use the 515F primer with truncLen=200,
#producing directly comparable 200bp ASVs over the same amplified region. Ilhan (truncLen=150,
#no primer trim) and Afshar (trimLeft=20, truncLen=200) produce ASVs of different length/region
#that cannot be meaningfully intersected at the exact-sequence level — genus-level comparison
#(compareStudies_16S.R) works across all 4 because genus classification is region-agnostic, but
#SV/ASV identity is not. As of this writing, all 4 phase-2 cohorts (03_dada2.sbatch) use
#PRIMER=NONE/TRIMLEFT=0 PLACEHOLDERS pending verification against real downloaded reads — none
#are confirmed to share BS/Assal's exact region, so none are included here by default. Add one
#only after confirming its actual primer/amplified-region matches (or accept that it changes
#what "same sequence" means and document that explicitly).

library(ggplot2)
library(gridExtra)
library(ggsignif)
library(ggrepel)
library(stringr)

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
dir.create(paste0(dirname(getwd()),"/input/"))
input = file.path(dirname(getwd()),"input/")
taxa<-"Seq"

#study -> (DADA module name fragment, list of (column suffix, display label) pairs).
#Confirmed region-compatible with each other: BS, Assal only (see note above).
STUDIES <- list(
  BS    = list(module="BSDADA",    timepoints=list(c("p1M","BS-1 month"), c("p6M","BS-6 months"))),
  Assal = list(module="AssalDADA", timepoints=list(c("p3M","Assal-3 months"), c("p1Y","Assal-1 year"), c("p2Y","Assal-2 years")))
)
#Same override mechanism as the other rewritten scripts, for when a phase-2 cohort's region IS
#confirmed compatible — e.g. PHASE2_STUDIES="BS,Assal,SomeConfirmedCohort" — but the default
#(unset) stays BS+Assal only, matching the original exactly.
active <- Sys.getenv("PHASE2_STUDIES", unset = "")
if (nzchar(active)) {
  keep <- trimws(strsplit(active, ",")[[1]])
  STUDIES <- STUDIES[keep]
}

for (study in names(STUDIES)) {
  cfg <- STUDIES[[study]]
  fileName <- paste0(taxa,"_",study,"_MixedLinearModelResults.txt")
  moduleDirPath <- file.path(dir(pipeRoot, cfg$module, full.names = TRUE))
  srcFile <- file.path(file.path(file.path(moduleDirPath, "output"),"MixedLinearModels"),fileName)
  file.copy(from = srcFile, to = input)
}

output = file.path(dirname(getwd()),"output/")
funcScript <- paste0(moduleDir,"/resources/functions.R")
source(funcScript)

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

for (s in 1:length(studies)){

  if (s!=length(studies)){

    otherStudies<-c((s+1):length(studies))

    for (s1 in otherStudies){

      path<-input
      df<-compareStudies(path,taxa,strsplit(studies[s],"-")[[1]][1],strsplit(studies[s1],"-")[[1]][1],strsplit(studies[s],"-")[[1]][2],strsplit(studies[s1],"-")[[1]][2])

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
      df<-compareStudies(path,taxa,strsplit(studies[s],"-")[[1]][1],strsplit(studies[s1],"-")[[1]][1],strsplit(studies[s],"-")[[1]][2],strsplit(studies[s1],"-")[[1]][2])
      xlab=paste0(studyNames[s]," vs. baseline")
      ylab=paste0(studyNames[s1]," vs. baseline")
      plot<-plotPairwiseStudiesMetagenomics(df,xlab,ylab,r[count],pval[count])
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
    }
  }

}

df<-data.frame(studyPairs,compariosn,compariosnTime,pval,r)
write.table(df, paste0(output,taxa,"_PairwiseComparison.txt"),sep="\t",row.names = FALSE)

p2<-wilcox.test(df$r[df$compariosn=="Different study" & df$compariosnTime=="Different timepoint"],
                df$r[df$compariosn=="Same study" & df$compariosnTime=="Different timepoint"])
capture.output(p2, file = paste0(output,taxa,"_wilcoxCoefficientResults.txt"))

plot1<-ggplot(data=df,aes(x=compariosn,y=r))+
  geom_boxplot(aes(color=compariosnTime),position = position_dodge(0.6), width = 0.5, size = 0.4,outlier.shape = NA)+
  geom_jitter(aes(color=compariosnTime),size=0.8,position = position_dodge(0.6))+labs(x="",y="Spearman Coefficient",color="")

#Generic pagination (see compareStudies_16S.R's note) — for BS+Assal's 10 pairs this still
#produces the original's exact page split (9 + 1), verified below.
pdf(paste0(output,taxa,"_scatterPlots.pdf"),width = 10,height = 10)
theme_set(theme_classic(base_size = 9))
n <- length(plotList)
for (i in seq(1, n, by = 9)) {
  pageEnd <- min(i+8, n)
  do.call(grid.arrange, c(plotList[i:pageEnd], ncol=3, nrow=3))
}
dev.off()

pdf(paste0(output,taxa,"_coefficientsFromScatterPlots.pdf"),width = 5,height = 5)
theme_set(theme_classic(base_size = 14))
print(plot1)
dev.off()

genusComparePath <- paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "CompareStudies16S"),"/output/")
df_g<-read.table(paste0(genusComparePath,"Genus_PairwiseComparison.txt"),sep="\t",header = TRUE)
df_g$taxanomy<-rep("Genus",nrow(df_g))
df$taxanomy<-rep("SV",nrow(df))
df_all<-rbind(df_g,df)

df_all$study1<-sapply(as.character(df_all$studyPairs),function(x){
  strsplit(x,"-")[[1]][1]})

df_all$study2<-sapply(sapply(as.character(df_all$studyPairs),function(x){
  strsplit(x,"_")[[1]][2]}),function(i){
    strsplit(i,"-")[[1]][1]
  })

#Restricts to the BS/Assal-only comparison regardless of how many total studies exist in
#Genus_PairwiseComparison.txt (which may include phase-2 cohorts) — this SV-vs-Genus comparison
#is only meaningful for the region-compatible pair, same restriction as STUDIES above.
df_sub<-df_all[(df_all$study1=="BS" | df_all$study1=="Assal") & (df_all$study2=="BS" | df_all$study2=="Assal"),  ]
plot2<-ggplot(data=df_sub,aes(x=compariosn,y=r))+
  geom_boxplot(aes(color=taxanomy),position = position_dodge(0.6), width = 0.5, size = 0.4,outlier.shape = NA)+
  geom_jitter(aes(color=taxanomy),size=0.8,position = position_dodge(0.6))+labs(x="",y="Spearman Coefficient",color="")

pdf(paste0(output,"coefficientsFromScatterPlots_compareSVandGenus.pdf"),width = 5,height = 5)
theme_set(theme_classic(base_size = 14))
print(plot2)
dev.off()

p1<-wilcox.test(df_sub$r[df_sub$compariosn=="Different study" & df_sub$taxanomy=="Genus"],
                df_sub$r[df_sub$compariosn=="Different study" & df_sub$taxanomy=="SV"])

p2<-wilcox.test(df_sub$r[df_sub$compariosn=="Different study" & df_sub$taxanomy=="Genus"],
                df_sub$r[df_sub$compariosn=="Same study" & df_sub$taxanomy=="Genus"])

p3<-wilcox.test(df_sub$r[df_sub$compariosn=="Different study" & df_sub$taxanomy=="Genus"],
                df_sub$r[df_sub$compariosn=="Same study" & df_sub$taxanomy=="SV"])

p4<-wilcox.test(df_sub$r[df_sub$compariosn=="Different study" & df_sub$taxanomy=="SV"],
                df_sub$r[df_sub$compariosn=="Same study" & df_sub$taxanomy=="Genus"])

p5<-wilcox.test(df_sub$r[df_sub$compariosn=="Different study" & df_sub$taxanomy=="SV"],
                df_sub$r[df_sub$compariosn=="Same study" & df_sub$taxanomy=="SV"])

p6<-wilcox.test(df_sub$r[df_sub$compariosn=="Same study" & df_sub$taxanomy=="Genus"],
                df_sub$r[df_sub$compariosn=="Same study" & df_sub$taxanomy=="SV"])

adjustedp<-p.adjust(c(p1$p.value,p2$p.value,p3$p.value,p4$p.value,p5$p.value,p6$p.value),method = "BH")

capture.output(p1, p2, p3, p4, p5, p6, adjustedp, file = paste0(output, taxa,"_wilcoxResults.txt"))
