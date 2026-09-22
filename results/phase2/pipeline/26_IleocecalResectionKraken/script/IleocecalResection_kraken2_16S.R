#Author: phase-2 extension (Rishi Gupta / Claude), modeled on Farnaz Fouladi's Assal_kraken2_16S.R
#Date: 2026
#Description: IleocecalResection dataset: compare the gut microbiome at each time point versus
#             baseline (month 0) using mixed linear models. Count tables classified by Kraken2.

library(nlme)

pipeRoot <- dirname(dirname(getwd()))
kraken_input <- paste0(pipeRoot,"/input/RYGB_IleocecalResection/Kraken2/")
meta_input <- paste0(pipeRoot,"/input/RYGB_IleocecalResection/Metadata/")
output <- file.path(dirname(getwd()),"output/")

dir.create(paste0(output,"MixedLinearModels/"))
dir.create(paste0(output,"MixedLinearModels/Kraken2/"))
output <- paste0(output,"MixedLinearModels/Kraken2/")

taxa<-c("Phylum","Class","Order","Family","Genus","Species")

for (t in taxa){

  myT<-read.table(paste0(kraken_input,"RYGB_IleocecalResection_2026_taxaCount_norm_Log10_",t,".tsv"),sep="\t",header = TRUE,row.names = 1,check.names = FALSE)

  meta<-read.table(paste0(meta_input,"metaData.txt"),sep="\t",header = TRUE)
  meta<-meta[meta$Run %in% rownames(myT),]
  myT<-myT[match(meta$Run,rownames(myT)),]
  meta$Timepoint<-factor(meta$Timepoint,levels = c(0,1,3,6))

  if(t=="Species"){
    myT2<-cbind(myT,meta)
    write.table(myT2,paste0(output,"IleocecalResection_speciesMetadata.txt"),sep="\t",quote = FALSE)
  }

  pval<-vector()
  p1<-vector()
  p3<-vector()
  p6<-vector()
  s1<-vector()
  s3<-vector()
  s6<-vector()
  bugName<-vector()
  index<-1

  for (i in 1:ncol(myT)){

    bug<-myT[,i]

    if (mean(bug>0)>0.1){

      df<-data.frame(bug,meta)

      fit<-tryCatch({anova(lme(bug~Timepoint,method="REML",random=~1|PatientID,data=df))},error=function(e){cat("LME ERROR:",conditionMessage(e),"\n");return(NA)})
      if(is.na(fit)[1]){ cat("  skipping non-converging column",i,"\n"); next }
      sm<-summary(lme(bug~Timepoint,method="REML",random=~1|PatientID,data=df))

      pval[index]<-fit$`p-value`[2]
      p1[index]<-sm$tTable[2,5]
      p3[index]<-sm$tTable[3,5]
      p6[index]<-sm$tTable[4,5]
      s1[index]<-sm$tTable[2,1]
      s3[index]<-sm$tTable[3,1]
      s6[index]<-sm$tTable[4,1]

      bugName[index]<-colnames(myT)[i]
      index<-index+1
    }
  }

  df<-data.frame(bugName,pval,p1,p3,p6,s1,s3,s6)
  df<-df[order(df$pval),]
  df$Adjustedpval<-p.adjust(df$pval,method = "BH")
  df$Adjustedp1<-p.adjust(df$p1,method = "BH")
  df$Adjustedp3<-p.adjust(df$p3,method = "BH")
  df$Adjustedp6<-p.adjust(df$p6,method = "BH")

  write.table(df,paste0(output,t,"_IleocecalResection_Kraken2_MixedLinearModelResults.txt"),sep="\t",row.names = FALSE,quote = FALSE)
}
