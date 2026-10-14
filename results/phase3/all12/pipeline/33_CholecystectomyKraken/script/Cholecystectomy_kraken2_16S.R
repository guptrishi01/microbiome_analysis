#Author: phase-2 extension (Rishi Gupta / Claude), modeled on Farnaz Fouladi's Ilhan_kraken2_16S.R
#Date: 2026
#Description: Cholecystectomy dataset: compare the gut microbiome at each time point versus
#             baseline using mixed linear models. Count tables classified by Kraken2.
#Note: reads metadata from Metadata/metaData.txt (not a duplicated Kraken2/ copy, unlike
#      Ilhan's original script) — see CLAUDE.md's phase-2 metadata-copy note.

library(nlme)

pipeRoot <- dirname(dirname(getwd()))
kraken_input <- paste0(pipeRoot,"/input/RYGB_Cholecystectomy/Kraken2/")
meta_input <- paste0(pipeRoot,"/input/RYGB_Cholecystectomy/Metadata/")
output <- file.path(dirname(getwd()),"output/")

dir.create(paste0(output,"MixedLinearModels/"))
dir.create(paste0(output,"MixedLinearModels/Kraken2/"))
output <- paste0(output,"MixedLinearModels/Kraken2/")

taxa<-c("Phylum","Class","Order","Family","Genus","Species")

for (t in taxa){

  myT<-read.table(paste0(kraken_input,"RYGB_Cholecystectomy_2026_taxaCount_norm_Log10_",t,".tsv"),sep="\t",header = TRUE,row.names = 1,check.names = FALSE)

  meta<-read.table(paste0(meta_input,"metaData.txt"),sep="\t",header = TRUE)
  meta<-meta[meta$Run %in% rownames(myT),]
  myT<-myT[match(meta$Run,rownames(myT)),]
  meta$time<-meta$Group

  meta1<-meta[meta$env_material=="fecal" & (meta$time=="Baseline" | meta$time=="6M" | meta$time=="12M"),]
  myT1<-myT[meta$env_material=="fecal" & (meta$time=="Baseline" | meta$time=="6M" | meta$time=="12M"),]

  meta1$time<-factor(meta1$time,levels = c("Baseline","6M","12M"))

  if(t=="Species"){
    myT2<-cbind(myT1,meta1)
    write.table(myT2,paste0(output,"Cholecystectomy_speciesMetadata.txt"),sep="\t",quote = FALSE)
  }

  pval<-vector()
  p6M<-vector()
  p12M<-vector()
  s6M<-vector()
  s12M<-vector()
  bugName<-vector()
  index<-1

  for (i in 1:ncol(myT1)){

    bug<-myT1[,i]

    if (mean(bug>0)>0.1){

      df<-data.frame(bug,meta1)

      fit<-tryCatch({anova(lme(bug~time,method="REML",random=~1|ID,data=df))},error=function(e){cat("LME ERROR:",conditionMessage(e),"\n");return(NA)})
      if(is.na(fit)[1]){ cat("  skipping non-converging column",i,"\n"); next }
      sm<-summary(lme(bug~time,method="REML",random=~1|ID,data=df))

      pval[index]<-fit$`p-value`[2]
      p6M[index]<-sm$tTable[2,5]
      p12M[index]<-sm$tTable[3,5]
      s6M[index]<-sm$tTable[2,1]
      s12M[index]<-sm$tTable[3,1]

      bugName[index]<-colnames(myT1)[i]
      index<-index+1
    }
  }

  df<-data.frame(bugName,pval,p6M,p12M,s6M,s12M)
  df<-df[order(df$pval),]
  df$Adjustedpval<-p.adjust(df$pval,method = "BH")
  df$Adjustedp6M<-p.adjust(df$p6M,method = "BH")
  df$Adjustedp12M<-p.adjust(df$p12M,method = "BH")

  write.table(df,paste0(output,t,"_Cholecystectomy_Kraken2_MixedLinearModelResults.txt"),sep="\t",row.names = FALSE,quote = FALSE)
}
