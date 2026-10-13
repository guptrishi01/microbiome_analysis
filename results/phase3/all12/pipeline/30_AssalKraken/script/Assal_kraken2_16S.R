#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: Assal dataset:Compare the gut microbiome at each time point versus
#             baseline using mixed linear models. Count tables is classified by 
#             Kraken2.

#BioLockJ configuration: Ali Sorgen
#Date: 10-06-2020

#Libraries
library(nlme)

pipeRoot <- dirname(dirname(getwd()))
input <- paste0(pipeRoot,"/input/RYGB_Assal2020/")
output <- file.path(dirname(getwd()),"output/")

dir.create(paste0(output,"MixedLinearModels/"))
dir.create(paste0(output,"MixedLinearModels/Kraken2/"))
output <- paste0(output,"MixedLinearModels/Kraken2/")

taxa<-c("Phylum","Class","Order","Family","Genus","Species")

for (t in taxa){
  
  myT<-read.table(paste0(input,"Kraken2/RYGB_Assal_2020Apr23_taxaCount_norm_Log10_",t,".tsv"),
                  sep="\t",header = TRUE,row.names = 1,check.names = FALSE)
  
  meta<-read.table(paste0(input,"Metadata/SraRunTable_Assal.txt"),sep=",",header = TRUE)
  myT<-myT[match(meta$Run,rownames(myT)),] 
  
  meta$ID<-as.factor(sapply(as.character(meta$Sample.Name),function(x){strsplit(x,"[.]")[[1]][1]}) )
  meta$time<-as.factor(sapply(as.character(meta$Sample.Name),function(x){strsplit(x,"[.]")[[1]][2]}) )
  meta$time<-factor(meta$time,levels = c("Pre","3M","1Y","2Y","3Y"))
  
  #Removing 3Y samples
  myT<-myT[meta$time!="3Y",]
  meta<-meta[meta$time!="3Y",]
  
  if (t=="Species"){
    
    #Removing low abundant Ehrlichia ruminantium(469) and 
    #Maricaulis maris (780) as they give convergence error
    
    myT<-myT[,-c(469,780)]
    myT2<-cbind(myT,meta)
    write.table(myT2,paste0(output,"Assal_speciesMetadata.txt"),sep="\t",quote = FALSE)
  }
  if(t=="Genus"){
    
    #Removing low abundant "Ehrlichia" (289), "Maricaulis" (465) 
    # as they give convergence error
    myT<-myT[,-c(289,465)]
  }
  
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
  
  write.table(df,paste0(output,t,"_Assal_Kraken2_MixedLinearModelResults.txt"),sep="\t",row.names = FALSE,quote = FALSE)
}

