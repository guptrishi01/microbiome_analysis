#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: Afshar dataset:Compare the gut microbiome at each time point versus
#             baseline using mixed linear models. Count tables is classified by 
#             Kraken2.

#BioLockJ configuration: Ali Sorgen
#Date: 10-06-2020

#Libraries
library(nlme)

pipeRoot <- dirname(dirname(getwd()))
input <- paste0(pipeRoot,"/input/RYGB_Afshar2018/Kraken2/")
output <- file.path(dirname(getwd()),"output/")

dir.create(paste0(output,"MixedLinearModels/"))
dir.create(paste0(output,"MixedLinearModels/Kraken2/"))
output <- paste0(output,"MixedLinearModels/Kraken2/")

taxa<-c("Phylum","Class","Order","Family","Genus","Species")

for (t in taxa){
  
  myT<-read.table(paste0(input,"RYGB_Afshar_2020Apr23_taxaCount_norm_Log10_",t,".tsv"),
                  sep="\t",header = TRUE,row.names = 1,check.names = FALSE)
  
  meta<-read.table(paste0(input,"SraRunTable_Afshar.txt"),sep=",",header = TRUE)
  myT<-myT[match(meta$Run,rownames(myT)),] 
  
  meta$ID<-as.factor(sapply(as.character(meta$title),function(x){strsplit(x,"_")[[1]][1]}) )
  meta$time<-as.factor(sapply(as.character(meta$title),function(x){strsplit(x,"_")[[1]][2]}) )
  
  #Select longitudinal samples
  meta1<-meta[!is.na(meta$time) & meta$time!="BOCABS",]
  myT1<-myT[!is.na(meta$time) & meta$time!="BOCABS",]
  
  meta1$time<-as.factor(as.character(meta1$time))
  meta1$time<-relevel(meta1$time,ref = "Pre")
  
  if(t=="Species"){
    myT2<-cbind(myT1,meta1)
    write.table(myT2,paste0(output,"Afshar_speciesMetadata.txt"),sep="\t",quote = FALSE)
  }
  
  pval<-vector()
  p6M<-vector()
  s6M<-vector()
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
        s6M[index]<-sm$tTable[2,1]
      
      
      bugName[index]<-colnames(myT1)[i]
      index<-index+1
    }
  }
  
  df<-data.frame(bugName,pval,p6M,s6M)
  df<-df[order(df$pval),]
  df$Adjustedpval<-p.adjust(df$pval,method = "BH")
  df$Adjustedp6M<-p.adjust(df$p6M,method = "BH")
  
  write.table(df,paste0(output,t,"_Afshar_Kraken2_MixedLinearModelResults.txt"),sep="\t",row.names = FALSE,quote = FALSE)
}


