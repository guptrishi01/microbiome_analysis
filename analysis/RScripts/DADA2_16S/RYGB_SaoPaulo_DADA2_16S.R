#Author: phase-3 extension (Rishi Gupta / Claude), modeled on the Cholecystectomy template
#Date: 2026
#Description: RYGB_SaoPaulo dataset: compare the gut microbiome at each time point versus
#              baseline using mixed linear models. Count tables are classified by DADA2.

library(nlme)
library(stringr)

pipeRoot <- dirname(dirname(getwd()))
input <- paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "RYGB_SaoPauloTaxaClass"),"/output/")
output <- file.path(dirname(getwd()),"output/")

dir.create(paste0(output,"MixedLinearModels/"))
output <- paste0(output,"MixedLinearModels/")

taxa<-c("Phylum","Class","Order","Family","Genus","SV")

for (t in taxa){

  dada2<-read.table(paste0(input,t,"_norm_table_RYGB_SaoPaulo.txt"),sep="\t",header = TRUE,row.names = 1,check.names = FALSE)

  finishAbundanceIndex<-which(colnames(dada2)=="BioProject")-1
  myT<-dada2[,1:finishAbundanceIndex]
  meta<-dada2[,(finishAbundanceIndex+1):ncol(dada2)]
  meta$Group<-factor(meta$Group,levels = c("T1","T2","T3"))

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

      fit<-anova(lme(bug~Group,method="REML",random=~1|ID,data=df))
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

  write.table(df,paste0(output,t,"_RYGB_SaoPaulo_MixedLinearModelResults.txt"),sep="\t",row.names = FALSE,quote = FALSE)
}
