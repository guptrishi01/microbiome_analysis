#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: BS dataset:Compare the gut microbiome at each time point versus
#             baseline using mixed linear models. Count tables is classified by 
#             Kraken2.
#Note: Streptomycetales and Streptomycetaceae produce singular convergence error 
#in mixed linear model

#BioLockJ configuration: Ali Sorgen
#Date: 10-06-2020

#Libraries
library(nlme)

pipeRoot <- dirname(dirname(getwd()))
input <- paste0(pipeRoot,"/input/RYGB_BS/")
output <- file.path(dirname(getwd()),"output/")

dir.create(paste0(output,"MixedLinearModels/"))
dir.create(paste0(output,"MixedLinearModels/Kraken2/"))
output <- paste0(output,"MixedLinearModels/Kraken2/")

metaData<-paste0(input,"Metadata/16S_metadata.txt")
taxa<-c("Phylum","Class","Order","Family","Genus","Species")


for (t in taxa){
  
  myT<-read.table(paste0(input,"Kraken2/RYGB_BS_2020Oct04_taxaCount_norm_Log10_",t,".tsv"),
                  sep="\t",header = TRUE,row.names = 1,check.names = FALSE)
  
  
  meta<-read.table(metaData,sep="\t",header = TRUE,row.names = 1)
  
  myT<-myT[match(rownames(meta),rownames(myT)),]
  
  if(t=="Species"){

    myT2<-cbind(myT,meta)
    write.table(myT2,paste0(output,"BS_speciesMetadata.txt"),sep="\t",quote = FALSE)
  }
  
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
      
      fit<-tryCatch({anova(lme(bug~factor(Timepoint),method="REML",random=~1|PatientID,data=df))},
                    error=function(e){cat("ERROR :",conditionMessage(e), "\n")
                      return(NA)})
      if(is.na(fit)[1]){
        pval[index]<-NA
        p1M[index]<-NA
        p6M[index]<-NA
        s1M[index]<-NA
        s6M[index]<-NA
        bugName[index]<-colnames(myT)[i]
        index<-index+1
        print(paste(t,colnames(myT)[i]))
        
      } else{
        
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
  }
  
  df<-data.frame(bugName,pval,p1M,p6M,s1M,s6M)
  df<-df[order(df$pval),]
  df$Adjustedpval<-p.adjust(df$pval,method = "BH")
  df$Adjustedp1M<-p.adjust(df$p1M,method = "BH")
  df$Adjustedp6M<-p.adjust(df$p6M,method = "BH")
  
  write.table(df,paste0(output,t,"_BS_Kraken2_MixedLinearModelResults.txt"),sep="\t",row.names = FALSE,quote = FALSE)
}

