#Author: WGS-taxonomy extension (Rishi Gupta / Claude), modeled on Farnaz Fouladi's
#        BS_kraken2_Metagenomics.R and this repo's Lazaro_GB_kraken2_16S.R.
#Date: 2026
#Description: Gopalakrishnan et al. 2025 / COCKTAIL study dataset (PRJNA904651, GBP arm
#             only — the VLED diet arm is excluded per this project's scope decision):
#             compare the gut microbiome at each post-op timepoint vs. baseline using
#             mixed linear models, on Kraken2+Bracken species-level count tables.
#
#NOT RUNNABLE YET — see cohorts/COCKTAIL_GBP/NEEDS_METADATA.md: every sample's time, ID,
#and Arm are UNKNOWN placeholders (the BioProject's ENA aliases are anonymized hashes by
#design; per-subject mapping requires a request to the COCKTAIL Steering Committee via
#the corresponding author). The paper's own visit schedule, for when that mapping
#arrives, is: Week-3 (pre-run-in-diet baseline), Week0 (post-run-in, pre-surgery), Week2,
#Week6, Month4, Year1, Year2 — seven timepoints, not the simple Pre/Post pattern used
#elsewhere in this repo. This script is written against that eventual 7-level structure
#so it's ready to run once metaData.txt is filled in; it is NOT ready now (every bug's
#mean(bug>0)>0.1 filter will still evaluate, but the lme() calls below will fail or
#produce meaningless output against all-UNKNOWN time/ID columns).
#
#Decide the correct reference level once real data exists: Week-3 captures the pre-diet
#state (used for the within-group baseline in Gopalakrishnan 2025's own PCoA/Fig 1C), but
#Week0 (post-run-in, immediately pre-surgery) may be the more appropriate "Pre" reference
#for a pure surgery-effect comparison — this repo defaults to Week-3 below only because
#that's what the paper itself uses as its reference row; revisit if that's not the
#comparison actually wanted.

library(nlme)

pipeRoot <- dirname(dirname(getwd()))
input <- paste0(pipeRoot,"/input/RYGB_COCKTAIL_GBP_WGS/Kraken2/")
meta_input <- paste0(pipeRoot,"/input/RYGB_COCKTAIL_GBP_WGS/Metadata/")
output <- file.path(dirname(getwd()),"output/")
dir.create(paste0(output,"MixedLinearModels/"))
dir.create(paste0(output,"MixedLinearModels/Kraken2/"))
output <- paste0(output,"MixedLinearModels/Kraken2/")

TIME_LEVELS <- c("Week-3","Week0","Week2","Week6","Month4","Year1","Year2")

taxa<-c("Phylum","Class","Order","Family","Genus","Species")

for (t in taxa){

  kraken<-read.table(paste0(input,"RYGB_COCKTAIL_GBP_2026_metagenomics_taxaCount_norm_Log10_",t,".tsv"),
                     sep="\t",header = TRUE,row.names = 1,check.names = FALSE)

  meta<-read.table(paste0(meta_input,"metaData.txt"),sep="\t",header = TRUE)
  meta<-meta[meta$Run %in% rownames(kraken),]
  kraken<-kraken[match(meta$Run,rownames(kraken)),]

  # R's any() on a zero-length logical vector returns FALSE, not NA/error — so if the
  # join above produced 0 matching rows (e.g. a renamed/typo'd Run column, or a count
  # table built for a different sample set), every UNKNOWN/Arm check below would
  # silently evaluate to FALSE and this stop() would never fire, continuing into
  # mean(bug>0) on a 0-row column (NaN) and if(NA) deep in the per-taxon loop instead of
  # this clear message. Guard it explicitly first.
  if (nrow(meta) == 0) {
    stop("metaData.txt has 0 rows matching this rank's count table's sample IDs — ",
         "check cohorts/COCKTAIL_GBP/metaData.txt's Run column against the Kraken2 output.")
  }
  if (any(meta$time == "UNKNOWN") || any(meta$ID == "UNKNOWN") || any(meta$Arm == "UNKNOWN")) {
    stop("metaData.txt still has UNKNOWN time/ID/Arm values — see cohorts/COCKTAIL_GBP/",
         "NEEDS_METADATA.md. Resolve the subject/timepoint mapping before running this script.")
  }
  if (any(meta$Arm != "GBP")) {
    stop("cohorts/COCKTAIL_GBP/metaData.txt has a non-GBP Arm value — this cohort should ",
         "be GBP-only (the VLED diet arm is explicitly out of scope, see this file's own ",
         "header comment); check metaData.txt wasn't filled in with the unfiltered BioProject.")
  }
  meta$time<-factor(meta$time,levels = TIME_LEVELS)

  if(t=="Species"){
    myT<-cbind(kraken,meta)
    write.table(myT,paste0(output,"COCKTAIL_GBP_speciesMetadata.txt"),sep="\t",quote = FALSE)
  }

  pval<-vector()
  coefs<-list()
  for (lvl in TIME_LEVELS[-1]) coefs[[paste0("p",lvl)]]<-vector()
  for (lvl in TIME_LEVELS[-1]) coefs[[paste0("s",lvl)]]<-vector()
  bugName<-vector()
  index<-1

  for (i in seq_len(ncol(kraken))){

    bug<-kraken[,i]

    if (mean(bug>0)>0.1){

      df<-data.frame(bug,meta)

      fit<-anova(lme(bug~time,method="REML",random=~1|ID,data=df))
      sm<-summary(lme(bug~time,method="REML",random=~1|ID,data=df))

      pval[index]<-fit$`p-value`[2]
      for (j in seq_along(TIME_LEVELS[-1])) {
        lvl <- TIME_LEVELS[-1][j]
        coefs[[paste0("p",lvl)]][index] <- sm$tTable[j+1,5]
        coefs[[paste0("s",lvl)]][index] <- sm$tTable[j+1,1]
      }

      bugName[index]<-colnames(kraken)[i]
      index<-index+1
    }
  }

  df<-data.frame(bugName,pval)
  for (nm in names(coefs)) df[[nm]]<-coefs[[nm]]
  df<-df[order(df$pval),]
  df$Adjustedpval<-p.adjust(df$pval,method = "BH")
  for (lvl in TIME_LEVELS[-1]) {
    df[[paste0("Adjustedp",lvl)]]<-p.adjust(df[[paste0("p",lvl)]],method = "BH")
  }

  write.table(df,paste0(output,t,"_COCKTAIL_GBP_Kraken2_Metagenomics_MixedLinearModelResults.txt"),sep="\t",row.names = FALSE,quote = FALSE)
}
