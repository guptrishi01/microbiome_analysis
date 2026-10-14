#Author: phase-2 extension (Rishi Gupta / Claude), modeled on Farnaz Fouladi's Ilhan_DADA2_16S.R
#Date: 2026
#Description: Cholecystectomy dataset: compare the gut microbiome at each time point versus
#              baseline using mixed linear models. Count tables are classified by DADA2.

library(nlme)
library(stringr)

pipeRoot <- dirname(dirname(getwd()))
input <- paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "CholecystectomyTaxaClass"),"/output/")
output <- file.path(dirname(getwd()),"output/")

dir.create(paste0(output,"MixedLinearModels/"))
output <- paste0(output,"MixedLinearModels/")

taxa<-c("Phylum","Class","Order","Family","Genus","SV")

for (t in taxa){

  dada2<-read.table(paste0(input,t,"_norm_table_Cholecystectomy.txt"),sep="\t",header = TRUE,row.names = 1,check.names = FALSE)

  finishAbundanceIndex<-which(colnames(dada2)=="BioProject")-1
  myT<-dada2[,1:finishAbundanceIndex]
  meta<-dada2[,(finishAbundanceIndex+1):ncol(dada2)]
  meta$Group<-factor(meta$Group,levels = c("Baseline","6M","12M"))

  pval<-vector()
  p6M<-vector()
  p12M<-vector()
  s6M<-vector()
  s12M<-vector()
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
      p6M[index]<-sm$tTable[2,5]
      p12M[index]<-sm$tTable[3,5]
      s6M[index]<-sm$tTable[2,1]
      s12M[index]<-sm$tTable[3,1]
      bugName[index]<-colnames(myT)[i]
      index<-index+1
    }
  }

  df<-data.frame(bugName,pval,p6M,p12M,s6M,s12M)
  df<-df[order(df$pval),]
  df$Adjustedpval<-p.adjust(df$pval,method = "BH")
  df$Adjustedp6M<-p.adjust(df$p6M,method = "BH")
  df$Adjustedp12M<-p.adjust(df$p12M,method = "BH")

  write.table(df,paste0(output,t,"_Cholecystectomy_MixedLinearModelResults.txt"),sep="\t",row.names = FALSE,quote = FALSE)
}
