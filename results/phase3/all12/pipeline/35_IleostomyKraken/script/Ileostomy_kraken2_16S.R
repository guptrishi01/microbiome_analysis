#Author: phase-2 extension (Rishi Gupta / Claude), modeled on Farnaz Fouladi's Afshar_kraken2_16S.R
#Date: 2026
#Description: Ileostomy dataset (PRJNA1480144, IL arm): compare the gut microbiome
#             Pre vs. Post using mixed linear models. Count tables classified by Kraken2.

library(nlme)

pipeRoot <- dirname(dirname(getwd()))
kraken_input <- paste0(pipeRoot,"/input/RYGB_Ileostomy/Kraken2/")
meta_input <- paste0(pipeRoot,"/input/RYGB_Ileostomy/Metadata/")
output <- file.path(dirname(getwd()),"output/")

dir.create(paste0(output,"MixedLinearModels/"))
dir.create(paste0(output,"MixedLinearModels/Kraken2/"))
output <- paste0(output,"MixedLinearModels/Kraken2/")

taxa<-c("Phylum","Class","Order","Family","Genus","Species")

for (t in taxa){

  myT<-read.table(paste0(kraken_input,"RYGB_Ileostomy_2026_taxaCount_norm_Log10_",t,".tsv"),sep="\t",header = TRUE,row.names = 1,check.names = FALSE)

  meta<-read.table(paste0(meta_input,"metaData.txt"),sep="\t",header = TRUE)
  meta<-meta[meta$Run %in% rownames(myT),]
  myT<-myT[match(meta$Run,rownames(myT)),]
  meta$time<-factor(meta$time,levels = c("Pre","Post"))

  if(t=="Species"){
    myT2<-cbind(myT,meta)
    write.table(myT2,paste0(output,"Ileostomy_speciesMetadata.txt"),sep="\t",quote = FALSE)
  }

  pval<-vector()
  pPost<-vector()
  sPost<-vector()
  bugName<-vector()
  index<-1

  for (i in 1:ncol(myT)){

    bug<-myT[,i]

    if (mean(bug>0)>0.1){

      df<-data.frame(bug,meta)

      fit<-tryCatch({anova(lme(bug~time,method="REML",random=~1|ID,data=df))},error=function(e){cat("LME ERROR:",conditionMessage(e),"\n");return(NA)})
      if(is.na(fit)[1]){ cat("  skipping non-converging column",i,"\n"); next }
      sm<-summary(lme(bug~time,method="REML",random=~1|ID,data=df))

      pval[index]<-fit$`p-value`[2]
      pPost[index]<-sm$tTable[2,5]
      sPost[index]<-sm$tTable[2,1]

      bugName[index]<-colnames(myT)[i]
      index<-index+1
    }
  }

  df<-data.frame(bugName,pval,pPost,sPost)
  df<-df[order(df$pval),]
  df$Adjustedpval<-p.adjust(df$pval,method = "BH")
  df$AdjustedpPost<-p.adjust(df$pPost,method = "BH")

  write.table(df,paste0(output,t,"_Ileostomy_Kraken2_MixedLinearModelResults.txt"),sep="\t",row.names = FALSE,quote = FALSE)
}
