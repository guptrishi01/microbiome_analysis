#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: BS dataset:Compare the gut microbiome at each time point versus
#              baseline using mixed linear models. Count tables are classified by 
#             DADA2.

#BioLockJ configuration: Ali Sorgen
#Date: 09-30-2020

#Libraries
library(nlme)
library(stringr)

pipeRoot <- dirname(dirname(getwd()))
input <- paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "BSTaxaClass"),"/output/")
output <- file.path(dirname(getwd()),"output/")

dir.create(paste0(output,"MixedLinearModels/"))
output <- paste0(output,"MixedLinearModels/")

taxa<-c("Phylum","Class","Order","Family","Genus","SV","Seq")

for (t in taxa){
  
  
  dada2<-read.table(paste0(input,t,"_norm_table_BS.txt"),sep="\t",header = TRUE,row.names = 1,check.names = FALSE)
   
  
  finishAbundanceIndex<-which(colnames(dada2)=="Barcode.ID" )-1
  myT<-dada2[,1:finishAbundanceIndex]
  meta<-dada2[,(finishAbundanceIndex+1):ncol(dada2)]

  pval<-vector()
  p1M<-vector()
  p6M<-vector()
  s1M<-vector()
  s6M<-vector()
  bugName<-vector()
  index<-1
  
  for (i in 1:ncol(myT)){
    
    bug<-myT[,i]
    
    if (mean(bug>0)>0.1){
      
      df<-data.frame(bug,meta)
      
      fit<-tryCatch({anova(lme(bug~factor(Timepoint),method="REML",random=~1|PatientID,data=df))},error=function(e){cat("LME ERROR:",conditionMessage(e),"\n");return(NA)})
      if(is.na(fit)[1]){ cat("  skipping non-converging column",i,"\n"); next }
      sm<-summary(lme(bug~factor(Timepoint),method="REML",random=~1|PatientID,data=df))
      
      pval[index]<-fit$`p-value`[2]
      p1M[index]<-sm$tTable[2,5]
      p6M[index]<-sm$tTable[3,5]
      s1M[index]<-sm$tTable[2,1]
      s6M[index]<-sm$tTable[3,1]
      bugName[index]<-colnames(myT)[i]
      index<-index+1
    }
  }
  
  df<-data.frame(bugName,pval,p1M,p6M,s1M,s6M)
  df<-df[order(df$pval),]
  df$Adjustedpval<-p.adjust(df$pval,method = "BH")
  df$Adjustedp1M<-p.adjust(df$p1M,method = "BH")
  df$Adjustedp6M<-p.adjust(df$p6M,method = "BH")
  
  write.table(df,paste0(output,t,"_BS_MixedLinearModelResults.txt"),sep="\t",row.names = FALSE,quote = FALSE)
}  
  
  