#Author: WGS-taxonomy extension (Rishi Gupta / Claude), modeled on Farnaz Fouladi's
#        BS_kraken2_Metagenomics.R and this repo's Lazaro_GB_kraken2_16S.R.
#Date: 2026
#Description: Hayashi et al. 2026 dataset (PRJNA1329344, laparoscopic sleeve gastrectomy
#             arm only — the intragastric balloon placement arm is excluded per this
#             project's scope decision, since IGBP is endoscopic, not a surgical
#             procedure): compare the gut microbiome at 6 months vs. baseline using
#             mixed linear models, on Kraken2+Bracken species-level count tables.
#
#NOT RUNNABLE YET — see cohorts/Hayashi_LSG/NEEDS_METADATA.md: every sample's time, ID,
#and Arm are UNKNOWN placeholders (ENA aliases are generic SampleNNN strings with no
#encoded info). The paper describes a simple two-timepoint design (baseline, 6 months)
#for 13 LSG patients (26 of the BioProject's 38 total runs, once the 6 IGBP patients'
#12 runs are excluded) — this script is written against that expected Pre/Post structure
#so it's ready to run once metaData.txt is filled in; it is NOT ready now.

library(nlme)

pipeRoot <- dirname(dirname(getwd()))
input <- paste0(pipeRoot,"/input/RYGB_Hayashi_LSG_WGS/Kraken2/")
meta_input <- paste0(pipeRoot,"/input/RYGB_Hayashi_LSG_WGS/Metadata/")
output <- file.path(dirname(getwd()),"output/")
dir.create(paste0(output,"MixedLinearModels/"))
dir.create(paste0(output,"MixedLinearModels/Kraken2/"))
output <- paste0(output,"MixedLinearModels/Kraken2/")

taxa<-c("Phylum","Class","Order","Family","Genus","Species")

for (t in taxa){

  kraken<-read.table(paste0(input,"RYGB_Hayashi_LSG_2026_metagenomics_taxaCount_norm_Log10_",t,".tsv"),
                     sep="\t",header = TRUE,row.names = 1,check.names = FALSE)

  meta<-read.table(paste0(meta_input,"metaData.txt"),sep="\t",header = TRUE)
  meta<-meta[meta$Run %in% rownames(kraken),]
  kraken<-kraken[match(meta$Run,rownames(kraken)),]

  if (any(meta$time == "UNKNOWN") || any(meta$ID == "UNKNOWN") || any(meta$Arm == "UNKNOWN")) {
    stop("metaData.txt still has UNKNOWN time/ID/Arm values — see cohorts/Hayashi_LSG/",
         "NEEDS_METADATA.md. Resolve the subject/timepoint/arm mapping (and filter to ",
         "LSG-only) before running this script.")
  }
  meta$time<-factor(meta$time,levels = c("Pre","Post"))

  if(t=="Species"){
    myT<-cbind(kraken,meta)
    write.table(myT,paste0(output,"Hayashi_LSG_speciesMetadata.txt"),sep="\t",quote = FALSE)
  }

  pval<-vector()
  pPost<-vector()
  sPost<-vector()
  bugName<-vector()
  index<-1

  for (i in 1:ncol(kraken)){

    bug<-kraken[,i]

    if (mean(bug>0)>0.1){

      df<-data.frame(bug,meta)

      fit<-anova(lme(bug~time,method="REML",random=~1|ID,data=df))
      sm<-summary(lme(bug~time,method="REML",random=~1|ID,data=df))

      pval[index]<-fit$`p-value`[2]
      pPost[index]<-sm$tTable[2,5]
      sPost[index]<-sm$tTable[2,1]

      bugName[index]<-colnames(kraken)[i]
      index<-index+1
    }
  }

  df<-data.frame(bugName,pval,pPost,sPost)
  df<-df[order(df$pval),]
  df$Adjustedpval<-p.adjust(df$pval,method = "BH")
  df$AdjustedpPost<-p.adjust(df$pPost,method = "BH")

  write.table(df,paste0(output,t,"_Hayashi_LSG_Kraken2_Metagenomics_MixedLinearModelResults.txt"),sep="\t",row.names = FALSE,quote = FALSE)
}
