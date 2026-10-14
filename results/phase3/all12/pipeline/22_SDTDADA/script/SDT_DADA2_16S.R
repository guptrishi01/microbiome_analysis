#Author: phase-2 extension (Rishi Gupta / Claude), modeled on Farnaz Fouladi's Afshar_DADA2_16S.R
#Date: 2026
#Description: SDT dataset (PRJNA1480144, SD arm): compare the gut microbiome
#             Pre vs. Post using mixed linear models. Count tables classified by DADA2.

library(nlme)
library(stringr)

pipeRoot <- dirname(dirname(getwd()))
input <- paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "SDTTaxaClass"),"/output/")
output <- file.path(dirname(getwd()),"output/")

dir.create(paste0(output,"MixedLinearModels/"))
output <- paste0(output,"MixedLinearModels/")
taxa<-c("Phylum","Class","Order","Family","Genus","SV")

for (t in taxa){

  dada2<-read.table(paste0(input,t,"_norm_table_SDT.txt"),sep="\t",header = TRUE,row.names = 1,check.names = FALSE)

  finishAbundanceIndex<-which(colnames(dada2)=="BioProject")-1
  myT<-dada2[,1:finishAbundanceIndex]
  meta<-dada2[,(finishAbundanceIndex+1):ncol(dada2)]
  meta$time<-factor(meta$time,levels = c("Pre","Post"))

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

  write.table(df,paste0(output,t,"_SDT_MixedLinearModelResults.txt"),sep="\t",row.names = FALSE,quote = FALSE)
}
