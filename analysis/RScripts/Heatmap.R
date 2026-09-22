#Author: phase-2 extension (Rishi Gupta / Claude), rewritten from Farnaz Fouladi's Heatmap.R
#Date: 2026
#Description: Heatmap for 16S datasets
#
#The original hardcodes exactly 4 studies (BS/Assal/Ilhan/Afshar) as separate read blocks and a
#fixed-length matrix() call. This is a generic rewrite over a STUDIES config (same shape as
#compareStudies_16S.R's, using short suffixes instead of full display labels since Heatmap.R's
#row names were always the short "Study-Suffix" form). See CLAUDE.md's "Never assume" section —
#verified against real phase-1 output before trusting for phase 2.

library(ComplexHeatmap)
library(stringr)

pipeRoot <- dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
outputDir <- file.path(moduleDir,"output")
resourceDir <- file.path(moduleDir, "resources")
funcScript <- file.path(resourceDir, "functions.R")
source(funcScript)
taxa <- "Genus"

#study -> (DADA module name fragment, list of column suffixes as they appear in
#p<suffix>/s<suffix> columns of that study's MixedLinearModelResults file)
STUDIES <- list(
  BS                 = list(module="BSDADA",                 suffixes=c("1M","6M")),
  Assal              = list(module="AssalDADA",               suffixes=c("3M","1Y","2Y")),
  Ilhan              = list(module="IlhanDADA",               suffixes=c("6M","1Y")),
  Afshar             = list(module="AfsharDADA",              suffixes=c("6M")),
  Cholecystectomy    = list(module="CholecystectomyDADA",     suffixes=c("6M","12M")),
  IleocecalResection = list(module="IleocecalResectionDADA",  suffixes=c("1","3","6")),
  Ileostomy          = list(module="IleostomyDADA",           suffixes=c("Post")),
  SDT                = list(module="SDTDADA",                 suffixes=c("Post"))
)
#Afshar's display suffix in row names is "Post" (not "6M") even though its data columns are
#named p6M/s6M — matches the original script's studies vector exactly.
DISPLAY_SUFFIX <- list(Afshar = list("6M"="Post"))

#8-color palette: BS/Assal/Ilhan/Afshar keep their original colors exactly (phase-1 output is
#unaffected); phase-2 cohorts get new, distinct colors.
STUDY_COLORS <- c(
  BS="yellow", Assal="pink", Ilhan="blue", Afshar="grey",
  Cholecystectomy="darkgreen", IleocecalResection="orange", Ileostomy="purple", SDT="brown"
)

active <- Sys.getenv("PHASE2_STUDIES", unset = "")
if (nzchar(active)) {
  keep <- trimws(strsplit(active, ",")[[1]])
  STUDIES <- STUDIES[keep]
}

studies<-character(0)
tables<-list()
for (study in names(STUDIES)) {
  cfg <- STUDIES[[study]]
  inDir <- file.path(file.path(dir(pipeRoot, cfg$module, full.names = TRUE),"output"),"MixedLinearModels")
  inFile <- paste0(taxa, "_", study, "_MixedLinearModelResults.txt")
  message("Reading file: ", file.path(inDir, inFile))
  myT <- read.table(file.path(inDir, inFile), sep="\t",header = TRUE,row.names = 1)
  for (suf in cfg$suffixes) {
    myT[[paste0("logp_",suf)]] <- getlog10p(myT[[paste0("p",suf)]], myT[[paste0("s",suf)]])
  }
  tables[[study]] <- myT
  for (suf in cfg$suffixes) {
    dispSuf <- suf
    if (!is.null(DISPLAY_SUFFIX[[study]]) && !is.null(DISPLAY_SUFFIX[[study]][[suf]])) {
      dispSuf <- DISPLAY_SUFFIX[[study]][[suf]]
    }
    studies <- c(studies, paste0(study,"-",dispSuf))
  }
}

#Intersect bug names across every active study (strict intersection, same as the original's
#nested intersect() calls, generalized to any number of studies).
bugs <- Reduce(intersect, lapply(tables, rownames))

rows <- list()
for (study in names(STUDIES)) {
  cfg <- STUDIES[[study]]
  myT_c <- tables[[study]][bugs,]
  for (suf in cfg$suffixes) {
    rows[[length(rows)+1]] <- myT_c[[paste0("logp_",suf)]]
  }
}

df<-matrix(unlist(rows),ncol = length(bugs),nrow=length(studies),byrow = TRUE)
rownames(df)<-studies
colnames(df)<-bugs

#Heatmap
df<-t(df)
studyOf<-sapply(colnames(df), function(x){strsplit(x,"-")[[1]][1]})

col = list(Study = STUDY_COLORS[unique(studyOf)])

ha <- HeatmapAnnotation(
  Study = studyOf, col = col
)

pdf(file.path(outputDir, paste0(taxa,"_HeatmapAndCluster.pdf")),height = 8)
Heatmap(df,top_annotation = ha,row_names_gp = gpar(fontsize = 6),heatmap_height = unit(20, "cm"),
        name = "log10 p-value")

dev.off()
