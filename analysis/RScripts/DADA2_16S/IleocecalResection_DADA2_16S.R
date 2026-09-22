#Author: phase-2 extension (Rishi Gupta / Claude), modeled on Farnaz Fouladi's Assal_DADA2_16S.R
#Date: 2026
#Description: IleocecalResection dataset: compare the gut microbiome at each time point versus
#              baseline (month 0) using mixed linear models. Count tables classified by DADA2.

library(nlme)
library(stringr)

pipeRoot <- dirname(dirname(getwd()))
input <- paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "IleocecalResectionTaxaClass"),"/output/")
output <- file.path(dirname(getwd()),"output/")

dir.create(paste0(output,"MixedLinearModels/"))
output <- paste0(output,"MixedLinearModels/")

taxa<-c("Phylum","Class","Order","Family","Genus","SV")

for (t in taxa){

  dada2<-read.table(paste0(input,t,"_norm_table_IleocecalResection.txt"),sep="\t",header = TRUE,row.names = 1,check.names = FALSE)

  finishAbundanceIndex<-which(colnames(dada2)=="BioProject")-1
  myT<-dada2[,1:finishAbundanceIndex]
  meta<-dada2[,(finishAbundanceIndex+1):ncol(dada2)]
  meta$Timepoint<-factor(meta$Timepoint,levels = c(0,1,3,6))

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

      fit<-anova(lme(bug~Timepoint,method="REML",random=~1|PatientID,data=df))
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

  write.table(df,paste0(output,t,"_IleocecalResection_MixedLinearModelResults.txt"),sep="\t",row.names = FALSE,quote = FALSE)
}
