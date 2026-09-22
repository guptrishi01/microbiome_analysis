#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: Functions


#Generates taxonomic tables
getTaxaTable<-function(svTable,taxaTable,taxa){
  colnames(svTable)<-taxaTable[,taxa]
  svTable<-t(svTable)
  tab<-rowsum(svTable,group=rownames(svTable))
  tab<-t(tab)
  colnames(tab)[ncol(tab)]<-"others"
  return(tab)
}

#Normalizes taxonomic tables
norm<-function(table){
  table<-table[rowSums(table)>1000,]
  average<-sum(rowSums(table))/nrow(table)
  table<-sweep(table,1,rowSums(table),"/")
  table<-log10(table*average + 1)
  return(table)
}

#This function gives the relative abundance when the input is log10 normalized count
getRelativeAbundance<-function(table){
  finishAAbundanceIndex<-which(colnames(table)=="others")
  tab1<-10^table[,1:finishAAbundanceIndex]
  tab2<-tab1-1
  tab3<-sweep(tab2,1,rowSums(tab2),"/")
  return(tab3)
}

#This function gives the count table (log10 normalized count)
getNormalizedCountTable<-function(table){
  finishAAbundanceIndex<-which(colnames(table)=="others")
  tab1<-table[,1:finishAAbundanceIndex]
  return(tab1)
}

#This function gives the metadata that is attached to the count table
getMetaData<-function(table){
  finishAAbundanceIndex<-which(colnames(table)=="others")
  tab<-table[,(finishAAbundanceIndex+1):ncol(table)]
  return(tab)
}

#This function gives a list of taxa in the count table for a given study
getListOfBug<-function(table,study,list){
  for (i in 1:ncol(table)){
    list[[paste0(colnames(table)[i],"_Study",study)]]<-table[,i]
    }
  return(list)
}

#This function generates an ordination plot
getPCO<-function(countTable,metaData,variable,names,colors){
  
  pco<-capscale(countTable~1,distance = "bray")
  percentVariance<-pco$CA$eig/sum(pco$CA$eig)*100
  
  pcoa=ordiplot(pco,choices=c(1,2),type="none",cex.lab=1,display="sites",xlab=paste0("PCoA ",format(percentVariance[1],digits = 3),"%"),
                ylab=paste0("PCoA ",format(percentVariance[2],digits = 3),"%"))
  col1<-colors[factor(metaData[,variable])]
  points(pcoa,"sites",col=adjustcolor(col1,alpha.f = 0.5),pch=16)
  
  variableLevels = length(levels(factor(metaData[,variable])))
  for (i in 1:variableLevels){
    ordiellipse(pcoa, factor(metaData[,variable]), kind="se", conf=0.95, lwd=2, draw = "lines", col=colors[i],
                show.groups=levels(factor(metaData[,variable]))[i],label=T,font=2,cex=0.7)
  }
  legend("bottomright",legend = names,col=colors[1:variableLevels],pch=16,cex=0.7)
  
}

#This function corrects the sign of p-values based on coefficient sign

getlog10p<-function(pval,coefficient){
  
  log10p<-sapply(1:length(pval), function(x){
    
    if (coefficient[x]>0) return(-log10(pval[x]))
    else return(log10(pval[x]))
  })
}

#This function gives the log10 p-values from studies or time points that we want to compare
compareStudies<-function(path,taxa,study1,study2,study1Time,study2Time){
  
  table1 = file.path(path, paste0(taxa,"_",study1,"_MixedLinearModelResults.txt"))
  message("Read table: ", table1)
  myT1<-read.table(table1,
                   sep="\t",header = TRUE,check.names = FALSE,quote = "",comment.char = "")
  table2 = file.path(path, paste0(taxa,"_",study2,"_MixedLinearModelResults.txt"))
  message("Read table: ", table2)
  myT2<-read.table(table2,
                   sep="\t",header = TRUE,check.names = FALSE,quote = "",comment.char = "")
  
  common<-intersect(myT1[,"bugName"],myT2[,"bugName"]) 
  myT1c<-myT1[myT1[,"bugName"] %in% common,]
  myT2c<-myT2[myT2[,"bugName"] %in% common,]
  myT2c<-myT2c[match(myT1c[,"bugName"],myT2c[,"bugName"]),]
  
  naColumns<-c(which(is.na(myT1c[,study1Time])),which(is.na(myT2c[,study2Time])))
  if(length(naColumns)>0){
    myT1c<-myT1c[-naColumns,]
    myT2c<-myT2c[-naColumns,]
  }
  
  s1<-paste0("s",substr(study1Time,2,nchar(study1Time)))
  s2<-paste0("s",substr(study2Time,2,nchar(study2Time)))
  
  pval1<-getlog10p(myT1c[,study1Time],myT1c[,s1])
  pval2<-getlog10p(myT2c[,study2Time],myT2c[,s2])
  
  df<-data.frame(pval1,pval2,bugName=myT2c[,"bugName"])
  
  return(df) 
}

#This function generate p-value versus p-value plot
plotPairwiseStudies<-function(df,xlab,ylab,coeficient,p){
  
  df1<-df[(df[,"pval1"]<log10(0.05) & df[,"pval2"]<log10(0.05)) | (df[,"pval1"]>-log10(0.05) & df[,"pval2"]>-log10(0.05)),]
  
  if (p==0){
    p1 = "< 2.2e-16"
  } else {
    p1 = paste("=",format(p,digits = 3))
  }
  
  theme_set(theme_classic())
  
  plot<-ggplot(data=df,aes(x=pval1,y=pval2))+geom_point(size=0.8)+
    geom_hline(yintercept = 0,linetype="dashed", color = "red")+
    geom_vline(xintercept = 0,linetype="dashed", color = "red")+
    labs(x=xlab,y=ylab,
         title = paste0("Spearman coefficient = ",format(coeficient,digits =3),"\nAdjusted p ",p1))+
    geom_text_repel(data=df1,aes(x=pval1,y=pval2,label=bugName),segment.colour="red",size=2.5,min.segment.length = 0,
                    segment.color="grey",segment.size=0.2)
  
}

#This function generate p-value versus p-value plot for the comparisons for metagenomic data
plotPairwiseStudiesMetagenomics<-function(df,xlab,ylab,coeficient,p){
  
  if (p==0){
    p1 = "< 2.2e-16"
  } else {
    p1 = paste("=",format(p,digits = 3))
  }
  
  theme_set(theme_classic())
  
  plot<-ggplot(data=df,aes(x=pval1,y=pval2))+geom_point(size=0.5,color=alpha("black",0.5))+
    geom_hline(yintercept = 0,linetype="dashed", color = "red")+
    geom_vline(xintercept = 0,linetype="dashed", color = "red")+
    labs(x=xlab,y=ylab,
         title = paste0("Spearman coefficient = ",format(coeficient,digits =3),"\nAdjusted p ",p1))
}

#This function performs Spearman test between studies
correlationBetweenStudies<-function(df){
  
  myList<-list()
  test<-cor.test(df[,"pval1"],df[,"pval2"],method = "spearman")
  myList[[1]]<-test$estimate
  myList[[2]]<-test$p.value
  return(myList)
}









