#Author: phase-3 extension (Rishi Gupta / Claude), modeled on the Cholecystectomy template
#Date: 2026
#Description: RYGB_SaoPaulo dataset: compare the gut microbiome at each time point versus
#             baseline using mixed linear models. Count tables classified by Kraken2.

library(nlme)

pipeRoot <- dirname(dirname(getwd()))
kraken_input <- paste0(pipeRoot,"/input/RYGB_SaoPaulo/Kraken2/")
meta_input <- paste0(pipeRoot,"/input/RYGB_SaoPaulo/Metadata/")
output <- file.path(dirname(getwd()),"output/")

dir.create(paste0(output,"MixedLinearModels/"))
dir.create(paste0(output,"MixedLinearModels/Kraken2/"))
output <- paste0(output,"MixedLinearModels/Kraken2/")

taxa<-c("Phylum","Class","Order","Family","Genus","Species")

for (t in taxa){

  myT<-read.table(paste0(kraken_input,"RYGB_SaoPaulo_2026_taxaCount_norm_Log10_",t,".tsv"),sep="\t",header = TRUE,row.names = 1,check.names = FALSE)

  meta<-read.table(paste0(meta_input,"metaData.txt"),sep="\t",header = TRUE)
  meta<-meta[meta$Run %in% rownames(myT),]
  myT<-myT[match(meta$Run,rownames(myT)),]
  meta$Group<-factor(meta$Group,levels = c("T1","T2","T3"))

  if(t=="Species"){
    myT2<-cbind(myT,meta)
    write.table(myT2,paste0(output,"RYGB_SaoPaulo_speciesMetadata.txt"),sep="\t",quote = FALSE)
  }

  pval<-vector()
  pT2<-vector()
  pT3<-vector()
  sT2<-vector()
  sT3<-vector()
  bugName<-vector()
  index<-1

  for (i in 1:ncol(myT)){

    bug<-myT[,i]

    if (mean(bug>0)>0.1){

      df<-data.frame(bug,meta)

      fit<-tryCatch({anova(lme(bug~Group,method="REML",random=~1|ID,data=df))},error=function(e){cat("LME ERROR:",conditionMessage(e),"\n");return(NA)})
      if(is.na(fit)[1]){ cat("  skipping non-converging column",i,"\n"); next }
      sm<-summary(lme(bug~Group,method="REML",random=~1|ID,data=df))

      pval[index]<-fit$`p-value`[2]
      pT2[index]<-sm$tTable[2,5]
      pT3[index]<-sm$tTable[3,5]
      sT2[index]<-sm$tTable[2,1]
      sT3[index]<-sm$tTable[3,1]

      bugName[index]<-colnames(myT)[i]
      index<-index+1
    }
  }

  df<-data.frame(bugName,pval,pT2,pT3,sT2,sT3)
  df<-df[order(df$pval),]
  df$Adjustedpval<-p.adjust(df$pval,method = "BH")
  df$AdjustedpT2<-p.adjust(df$pT2,method = "BH")
  df$AdjustedpT3<-p.adjust(df$pT3,method = "BH")

  write.table(df,paste0(output,t,"_RYGB_SaoPaulo_Kraken2_MixedLinearModelResults.txt"),sep="\t",row.names = FALSE,quote = FALSE)
}
