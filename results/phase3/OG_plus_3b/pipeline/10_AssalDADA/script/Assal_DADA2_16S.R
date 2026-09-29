#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: Assal dataset:Compare the gut microbiome at each time point versus
#              baseline using mixed linear models.Count tables are classified by 
#             DADA2.

#BioLockJ configuration: Ali Sorgen
#Date: 09-30-2020

#Libraries
library(nlme)
library(stringr)

pipeRoot <- dirname(dirname(getwd()))
input <- paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "AssalTaxaClass"),"/output/")
output <- file.path(dirname(getwd()),"output/")

dir.create(paste0(output,"MixedLinearModels/"))
output <- paste0(output,"MixedLinearModels/")

taxa<-c("Phylum","Class","Order","Family","Genus","SV","Seq")

for (t in taxa){
  
  dada2<-read.table(paste0(input,t,"_norm_table_Assal.txt"),sep="\t",header = TRUE,row.names = 1,check.names = FALSE)
  
  finishAbundanceIndex<-which(colnames(dada2)=="Assay.Type")-1
  myT<-dada2[,1:finishAbundanceIndex]
  meta<-dada2[,(finishAbundanceIndex+1):ncol(dada2)]
  meta$time<-factor(meta$time,levels = c("Pre","3M","1Y","2Y"))
  
  pval<-vector()
  p3M<-vector()
  p1Y<-vector()
  p2Y<-vector()
  s3M<-vector()
  s1Y<-vector()
  s2Y<-vector()
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
      p3M[index]<-sm$tTable[2,5]
      p1Y[index]<-sm$tTable[3,5]
      p2Y[index]<-sm$tTable[4,5]
      s3M[index]<-sm$tTable[2,1]
      s1Y[index]<-sm$tTable[3,1]
      s2Y[index]<-sm$tTable[4,1]
      
      bugName[index]<-colnames(myT)[i]
      index<-index+1
    }
  }
  
  df<-data.frame(bugName,pval,p3M,p1Y,p2Y,s3M,s1Y,s2Y)
  df<-df[order(df$pval),]
  df$Adjustedpval<-p.adjust(df$pval,method = "BH")
  df$Adjustedp3M<-p.adjust(df$p3M,method = "BH")
  df$Adjustedp1Y<-p.adjust(df$p1Y,method = "BH")
  df$Adjustedp2Y<-p.adjust(df$p2Y,method = "BH")
  
  write.table(df,paste0(output,t,"_Assal_MixedLinearModelResults.txt"),sep="\t",row.names = FALSE,quote = FALSE)
}  
