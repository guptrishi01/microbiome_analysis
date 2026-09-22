#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: Ilhan dataset:Compare the gut microbiome at each time point versus
#             baseline using mixed linear models. Count tables is classified by 
#             Kraken2.

#BioLockJ configuration: Ali Sorgen
#Date: 10-06-2020

#Libraries
library(nlme)

pipeRoot <- dirname(dirname(getwd()))
input <- paste0(pipeRoot,"/input/RYGB_Ilhan2020/Kraken2/")
output <- file.path(dirname(getwd()),"output/")

dir.create(paste0(output,"MixedLinearModels/"))
dir.create(paste0(output,"MixedLinearModels/Kraken2/"))
output <- paste0(output,"MixedLinearModels/Kraken2/")

taxa<-c("Phylum","Class","Order","Family","Genus","Species")


for (t in taxa){
  
  myT<-read.table(paste0(input,"RYGB_Ilhan_2020Apr23_taxaCount_norm_Log10_",t,".tsv"),sep="\t",header = TRUE,row.names = 1,check.names = FALSE)
  
  meta<-read.table(paste0(input,"metaData.txt"),sep="\t",header = TRUE)
  meta<-meta[meta$Run %in% rownames(myT),]
  myT<-myT[match(meta$Run,rownames(myT)),] 
  meta$time<-meta$Group
  
  meta1<-meta[meta$env_material=="fecal" & (meta$time=="Baseline" | meta$time=="6M" | meta$time=="12M"),]
  myT1<-myT[meta$env_material=="fecal" & (meta$time=="Baseline" | meta$time=="6M" | meta$time=="12M"),]
  
  meta1$time<-factor(meta1$time,levels = c("Baseline","6M","12M"))
  
  if(t=="Species"){
    myT2<-cbind(myT1,meta1)
    write.table(myT2,paste0(output,"Ilhan_speciesMetadata.txt"),sep="\t",quote = FALSE)
  }
  
  pval<-vector()
  p6M<-vector()
  p1Y<-vector()
  s6M<-vector()
  s1Y<-vector()
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
      p1Y[index]<-sm$tTable[3,5]
      s6M[index]<-sm$tTable[2,1]
      s1Y[index]<-sm$tTable[3,1]
      
      bugName[index]<-colnames(myT1)[i]
      index<-index+1
      
    }
  }
  
  df<-data.frame(bugName,pval,p6M,p1Y,s6M,s1Y)
  df<-df[order(df$pval),]
  df$Adjustedpval<-p.adjust(df$pval,method = "BH")
  df$Adjustedp6M<-p.adjust(df$p6M,method = "BH")
  df$Adjustedp1Y<-p.adjust(df$p1Y,method = "BH")
  
  write.table(df,paste0(output,t,"_Ilhan_Kraken2_MixedLinearModelResults.txt"),sep="\t",row.names = FALSE,quote = FALSE)
}

