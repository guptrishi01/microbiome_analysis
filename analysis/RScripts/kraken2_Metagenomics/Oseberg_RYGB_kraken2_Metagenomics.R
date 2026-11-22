#Author: WGS-taxonomy extension (Rishi Gupta / Claude), modeled on Farnaz Fouladi's
#        BS_kraken2_Metagenomics.R and this repo's Lazaro_GB_kraken2_16S.R.
#Date: 2026
#Description: Olsson et al. 2026 / Oseberg trial dataset (PRJEB59861), RYGB arm only
#             (n=39 patients, 78 runs): compare the gut microbiome at 12 months vs.
#             baseline using mixed linear models, on Kraken2+Bracken species-level
#             (and other ranks) count tables.
#
#RESOLVED 2026-11-14: this BioProject originally pooled RYGB and sleeve-gastrectomy
#patients from one randomized trial with no public per-patient procedure label. Arm was
#resolved via the paper's own Zenodo source-data deposit (10.5281/zenodo.19203868, linked
#from Olsson et al. 2026's Data Availability Statement) — a `mappP1.rds` data frame with
#per-patient `op2` (RYGB/SG) factor, verified against the paper's own reported arm sizes
#(39 RYGB + 38 SG) and diffed against every patient ID and Run accession already in
#cohorts/Oseberg/metaData.txt with zero mismatches. cohorts/Oseberg/ was then split into
#cohorts/Oseberg_RYGB/ and cohorts/Oseberg_SG/, and this script (the former unsplit
#Oseberg_kraken2_Metagenomics.R) into this RYGB-only copy and its SG-only sibling, per
#this repo's "one study, one procedure" rule (CLAUDE.md).

library(nlme)

pipeRoot <- dirname(dirname(getwd()))
input <- paste0(pipeRoot,"/input/RYGB_Oseberg_RYGB_WGS/Kraken2/")
meta_input <- paste0(pipeRoot,"/input/RYGB_Oseberg_RYGB_WGS/Metadata/")
output <- file.path(dirname(getwd()),"output/")
dir.create(paste0(output,"MixedLinearModels/"))
dir.create(paste0(output,"MixedLinearModels/Kraken2/"))
output <- paste0(output,"MixedLinearModels/Kraken2/")

taxa<-c("Phylum","Class","Order","Family","Genus","Species")

for (t in taxa){

  kraken<-read.table(paste0(input,"RYGB_Oseberg_RYGB_2026_metagenomics_taxaCount_norm_Log10_",t,".tsv"),
                     sep="\t",header = TRUE,row.names = 1,check.names = FALSE)

  meta<-read.table(paste0(meta_input,"metaData.txt"),sep="\t",header = TRUE)
  meta<-meta[meta$Run %in% rownames(kraken),]
  kraken<-kraken[match(meta$Run,rownames(kraken)),]

  # R's any() on a zero-length logical vector returns FALSE, not NA/error — so if the
  # join above produced 0 matching rows, the Arm check below would silently evaluate to
  # FALSE and never fire. Guard it explicitly first.
  if (nrow(meta) == 0) {
    stop("metaData.txt has 0 rows matching this rank's count table's sample IDs — ",
         "check cohorts/Oseberg_RYGB/metaData.txt's Run column against the Kraken2 output.")
  }
  if (any(meta$Arm != "RYGB")) {
    stop("cohorts/Oseberg_RYGB/metaData.txt has a non-RYGB Arm value — this cohort should ",
         "be RYGB-only post-split; check it wasn't accidentally pointed at unfiltered data.")
  }
  # Real timepoints decoded from ENA sample_alias (see NEEDS_METADATA.md); B (baseline) is reference.
  meta$time<-factor(meta$time,levels = c("B","12M"))

  if(t=="Species"){
    myT<-cbind(kraken,meta)
    write.table(myT,paste0(output,"Oseberg_RYGB_speciesMetadata.txt"),sep="\t",quote = FALSE)
  }

  pval<-vector()
  p12M<-vector()
  s12M<-vector()
  bugName<-vector()
  index<-1

  for (i in seq_len(ncol(kraken))){

    bug<-kraken[,i]

    if (mean(bug>0)>0.1){

      df<-data.frame(bug,meta)

      fit<-anova(lme(bug~time,method="REML",random=~1|ID,data=df))
      sm<-summary(lme(bug~time,method="REML",random=~1|ID,data=df))

      pval[index]<-fit$`p-value`[2]
      p12M[index]<-sm$tTable[2,5]
      s12M[index]<-sm$tTable[2,1]

      bugName[index]<-colnames(kraken)[i]
      index<-index+1
    }
  }

  df<-data.frame(bugName,pval,p12M,s12M)
  df<-df[order(df$pval),]
  df$Adjustedpval<-p.adjust(df$pval,method = "BH")
  df$Adjustedp12M<-p.adjust(df$p12M,method = "BH")

  write.table(df,paste0(output,t,"_Oseberg_RYGB_Kraken2_Metagenomics_MixedLinearModelResults.txt"),sep="\t",row.names = FALSE,quote = FALSE)
}
