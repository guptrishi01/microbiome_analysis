#Author: WGS-taxonomy extension (Rishi Gupta / Claude), modeled on Farnaz Fouladi's
#        BS_kraken2_Metagenomics.R and this repo's Lazaro_GB_kraken2_16S.R.
#Date: 2026
#Description: Penney et al. 2022 dataset (PRJNA473348), sleeve-gastrectomy arm only
#             (confirmed subset: 14 patients, 31 runs — matches the paper's full reported
#             SG arm exactly; see cohorts/Penney_SG/NEEDS_METADATA.md): compare the gut
#             microbiome at each post-op timepoint vs. baseline using mixed linear models,
#             on Kraken2+Bracken species-level (and other ranks) count tables.
#
#RESOLVED 2026-11-14: this BioProject originally pooled RYGB and sleeve-gastrectomy
#patients with no public per-patient procedure label. A confirmed 14-patient SG subset
#was recovered from the paper's own Mendeley Data deposit (10.17632/t76nm3yfzh.3,
#Source Data/Figure 6/MixOmics_PairedVSG.xlsx's Sample_Name column; the paper calls this
#arm "VSG" — mapped to this repo's standard "SG" label for consistency with
#Oseberg_SG/Lazaro_SG/Chaudhari_SG) — verified against the ID space already decoded from
#ENA sample_alias (exact string match, no coincidence) and matches the paper's own
#reported n=14 SG in its "27 complete paired sets" figure exactly. ~54 other Penney
#patients with no Arm anywhere in the public deposit are not in this cohort — see
#NEEDS_METADATA.md. cohorts/Penney/ was split into cohorts/Penney_RYGB/ and
#cohorts/Penney_SG/, and this script (the former unsplit Penney_kraken2_Metagenomics.R)
#into this SG-only copy and its RYGB-only sibling, per this repo's "one study, one
#procedure" rule (CLAUDE.md).

library(nlme)

pipeRoot <- dirname(dirname(getwd()))
input <- paste0(pipeRoot,"/input/RYGB_Penney_SG_WGS/Kraken2/")
meta_input <- paste0(pipeRoot,"/input/RYGB_Penney_SG_WGS/Metadata/")
output <- file.path(dirname(getwd()),"output/")
dir.create(paste0(output,"MixedLinearModels/"))
dir.create(paste0(output,"MixedLinearModels/Kraken2/"))
output <- paste0(output,"MixedLinearModels/Kraken2/")

taxa<-c("Phylum","Class","Order","Family","Genus","Species")

for (t in taxa){

  kraken<-read.table(paste0(input,"RYGB_Penney_SG_2026_metagenomics_taxaCount_norm_Log10_",t,".tsv"),
                     sep="\t",header = TRUE,row.names = 1,check.names = FALSE)

  meta<-read.table(paste0(meta_input,"metaData.txt"),sep="\t",header = TRUE)
  meta<-meta[meta$Run %in% rownames(kraken),]
  kraken<-kraken[match(meta$Run,rownames(kraken)),]

  # R's any() on a zero-length logical vector returns FALSE, not NA/error — so if the
  # join above produced 0 matching rows, the Arm check below would silently evaluate to
  # FALSE and never fire. Guard it explicitly first.
  if (nrow(meta) == 0) {
    stop("metaData.txt has 0 rows matching this rank's count table's sample IDs — ",
         "check cohorts/Penney_SG/metaData.txt's Run column against the Kraken2 output.")
  }
  if (any(meta$Arm != "SG")) {
    stop("cohorts/Penney_SG/metaData.txt has a non-SG Arm value — this cohort should ",
         "be SG-only post-split; check it wasn't accidentally pointed at unfiltered data.")
  }
  # Real timepoints decoded from ENA sample_alias (see NEEDS_METADATA.md); Pre is baseline.
  meta$time<-factor(meta$time,levels = c("Pre","POST","1YrPOST"))

  if(t=="Species"){
    myT<-cbind(kraken,meta)
    write.table(myT,paste0(output,"Penney_SG_speciesMetadata.txt"),sep="\t",quote = FALSE)
  }

  pval<-vector()
  pPOST<-vector()
  p1YrPOST<-vector()
  sPOST<-vector()
  s1YrPOST<-vector()
  bugName<-vector()
  index<-1

  for (i in seq_len(ncol(kraken))){

    bug<-kraken[,i]

    if (mean(bug>0)>0.1){

      df<-data.frame(bug,meta)

      fit<-anova(lme(bug~time,method="REML",random=~1|ID,data=df))
      sm<-summary(lme(bug~time,method="REML",random=~1|ID,data=df))

      pval[index]<-fit$`p-value`[2]
      pPOST[index]<-sm$tTable[2,5]
      p1YrPOST[index]<-sm$tTable[3,5]
      sPOST[index]<-sm$tTable[2,1]
      s1YrPOST[index]<-sm$tTable[3,1]

      bugName[index]<-colnames(kraken)[i]
      index<-index+1
    }
  }

  df<-data.frame(bugName,pval,pPOST,p1YrPOST,sPOST,s1YrPOST)
  df<-df[order(df$pval),]
  df$Adjustedpval<-p.adjust(df$pval,method = "BH")
  df$AdjustedpPOST<-p.adjust(df$pPOST,method = "BH")
  df$Adjustedp1YrPOST<-p.adjust(df$p1YrPOST,method = "BH")

  write.table(df,paste0(output,t,"_Penney_SG_Kraken2_Metagenomics_MixedLinearModelResults.txt"),sep="\t",row.names = FALSE,quote = FALSE)
}
