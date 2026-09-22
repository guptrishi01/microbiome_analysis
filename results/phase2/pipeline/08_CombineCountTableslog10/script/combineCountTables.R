#Author: phase-2 extension (Rishi Gupta / Claude), rewritten from Farnaz Fouladi's combineCountTables.R
#Date: 2026
#Description: Join 16S table across studies. This script works only on the genus level.
#
#The original hardcodes exactly 4 studies (BS/Assal/Ilhan/Afshar) as separate variables with
#literal matrix row-index arithmetic — see CLAUDE.md's "Adding a new cohort" section. This is
#a generic rewrite over a STUDIES list so a phase-2 cohort is one new list entry, not a new
#hand-written block. Behavior for the original 4 studies is unchanged — verified by running
#this script with only those 4 entries against the real phase-1 pipeline output and diffing
#the result byte-for-byte against results/phase1/pipeline/*_CombineCountTables*/output/
#(see scripts/ for the verification harness, or tests/verify_combine_count_tables.sh).

library(stringr)

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
output = file.path(dirname(getwd()),"output/")

args <- commandArgs(trailingOnly = TRUE)
taxa<-"Genus"
normalized <- args[1]

funcScript <- paste0(moduleDir,"/resources/functions.R")
source(funcScript)

#study -> (TaxaClass module name fragment, timepoint column name, subject ID column name)
#Order matters only for matching legacy row ordering when comparing against phase-1 output.
STUDIES <- list(
  BS                 = list(module="BSTaxaClass",                 timepointCol="Timepoint", idCol="PatientID"),
  Assal              = list(module="AssalTaxaClass",               timepointCol="time",      idCol="ID"),
  Ilhan              = list(module="IlhanTaxaClass",               timepointCol="Group",     idCol="ID"),
  Afshar             = list(module="AfsharTaxaClass",              timepointCol="time",      idCol="ID"),
  Cholecystectomy    = list(module="CholecystectomyTaxaClass",     timepointCol="Group",     idCol="ID"),
  IleocecalResection = list(module="IleocecalResectionTaxaClass",  timepointCol="Timepoint", idCol="PatientID"),
  Ileostomy          = list(module="IleostomyTaxaClass",           timepointCol="time",      idCol="ID"),
  SDT                = list(module="SDTTaxaClass",                 timepointCol="time",      idCol="ID")
)
#PHASE2_STUDIES env var restricts to a subset (e.g. "BS,Assal,Ilhan,Afshar" for a phase-1-only
#verification run); unset/empty means all studies in STUDIES above.
active <- Sys.getenv("PHASE2_STUDIES", unset = "")
if (nzchar(active)) {
  keep <- trimws(strsplit(active, ",")[[1]])
  STUDIES <- STUDIES[keep]
}

tables <- list()
for (study in names(STUDIES)) {
  cfg <- STUDIES[[study]]
  path <- paste0(pipeRoot,"/",str_subset(dir(pipeRoot), cfg$module),"/output/")
  t1 <- read.table(paste0(path,taxa,"_norm_table_",study,".txt"),sep="\t",header=TRUE)
  if (normalized=="relab") {
    t1n <- getRelativeAbundance(t1)
  } else {
    t1n <- getNormalizedCountTable(t1)
  }
  tables[[study]] <- list(counts=t1n, meta=getMetaData(t1), cfg=cfg)
}

#Get a list of all taxa in all count tables
bugs <- list()
for (study in names(tables)) {
  bugs <- getListOfBug(tables[[study]]$counts, study, bugs)
}
namesOfTheList <- names(bugs)
namesOfTheBugs <- sapply(namesOfTheList, function(x) strsplit(x,"_Study")[[1]][1])
uniqueBugs <- unique(namesOfTheBugs)

totalRows <- sum(sapply(tables, function(x) nrow(x$counts)))
allRownames <- unlist(lapply(tables, function(x) rownames(x$counts)))
mat <- matrix(NA, nrow=totalRows, ncol=length(uniqueBugs), dimnames=list(allRownames, uniqueBugs))

rowOffset <- 0
for (study in names(tables)) {
  n <- nrow(tables[[study]]$counts)
  rowsHere <- (rowOffset+1):(rowOffset+n)
  for (bugName in uniqueBugs) {
    key <- paste0(bugName,"_Study",study)
    if (key %in% namesOfTheList) {
      mat[rowsHere, bugName] <- bugs[[key]]
    } else {
      mat[rowsHere, bugName] <- 0
    }
  }
  rowOffset <- rowOffset + n
}

df <- as.data.frame(mat)

#NOTE: ID is deliberately NOT as.character()'d here, matching the original script's own
#c(meta.BS$PatientID, meta.Assal$ID, meta.Ilhan$ID, meta.Afshar$ID) exactly (which lacks the
#as.character() wrap that `timepoint` gets). Combining multiple R factors with base c() collapses
#them to integer codes, not string labels — almost certainly an unintentional quirk in the
#original, but the published phase-1 metaData_merged.txt reflects it, so byte-identical
#reproduction requires keeping it for BS/Assal/Ilhan/Afshar. Verified via a direct diff against
#the real phase-1 output (see the module docstring at the top of this file).
meta_all <- data.frame(
  time = unlist(lapply(tables, function(x) x$meta$prepost)),
  timepoint = unlist(lapply(names(tables), function(s) as.character(tables[[s]]$meta[[tables[[s]]$cfg$timepointCol]]))),
  ID = do.call(c, lapply(names(tables), function(s) tables[[s]]$meta[[tables[[s]]$cfg$idCol]])),
  Study = rep(names(tables), sapply(tables, function(x) nrow(x$counts))),
  Sample_ID = allRownames
)
rownames(meta_all) <- meta_all$Sample_ID

write.table(df,paste0(output,taxa,"_countTable_merged_",normalized,".txt"),sep="\t",quote = FALSE)
write.table(meta_all,paste0(output,"metaData_merged.txt"),sep="\t",quote = FALSE)
