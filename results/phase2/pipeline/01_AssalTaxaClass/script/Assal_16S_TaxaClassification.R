#Author: Farnaz Fouladi
#Date: 10-01-2020
#Description: This script generates taxonomic tables

#BioLockJ configuration: Ali Sorgen
#Date: 10-06-2020

pipeRoot <- dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
path <- paste0(pipeRoot,"/input/RYGB_Assal2020/")
output <- file.path(dirname(getwd()),"output/")

metaData<-paste0(path,"Metadata/SraRunTable_Assal.txt")
funcScript <- paste0(moduleDir,"/resources/functions.R")
source(funcScript)

dada<-read.table(paste0(path,"DADA2/ForwardReads.txt"),sep="\t",header=TRUE)
taxa<-read.table(paste0(path,"DADA2/taxForwardReads.txt"),sep="\t",header=TRUE)
meta<-read.table(metaData,sep=",",header=TRUE,row.names = 1)

#Modify Escherichia name
taxa$Genus = as.factor(taxa$Genus)
typoIndex = which(levels(taxa$Genus)=="Esherichica/Shigella")
levels(taxa$Genus)[typoIndex]="Escherichia/Shigella"

dada<-dada[match(rownames(meta),rownames(dada)),]

meta$time<-as.factor(sapply(as.character(meta$Sample.Name),function(x){strsplit(x,"[.]")[[1]][2]}))
meta$ID<-as.factor(sapply(as.character(meta$Sample.Name),function(x){strsplit(x,"[.]")[[1]][1]}))
#Remove time point 3 years due to small sample size
meta1<-meta[meta$time!="3Y",]
dada<-dada[meta$time!="3Y",]
meta1$time<-factor(meta1$time,levels = c("Pre","3M","1Y","2Y"))
meta1$prepost<-sapply(meta1$time,function(x){if (x=="Pre") return(0) else return(1)})

if(sum(rownames(dada)==rownames(meta1))!=nrow(dada)) stop("Error")


for(t in c("Phylum","Class","Order","Family","Genus")){
  
  t1<-getTaxaTable(dada,taxa,t)
  t1_norm<-norm(t1)
  t1_normMeta<-cbind(t1_norm,meta1)
  write.table(t1_normMeta,paste0(output,t,"_norm_table_Assal.txt"),sep = "\t",row.names = TRUE,quote = FALSE)
}

#Write normalized sequence variant tables 
dada1<-norm(dada)
dada1_meta<-cbind(dada1,meta1)
write.table(dada1_meta,paste0(output,"Seq_norm_table_Assal.txt"),sep="\t",row.names = TRUE,quote = FALSE)

#SV table
num<-c(1:nrow(taxa))
taxanomy<-apply(taxa,1,function(x){paste0(x[1],"_",x[2],"_",x[3],"_",x[4],"_",x[5],"_",x[6])})
taxanomy<-paste0(taxanomy,"_",num)
colnames(dada1)<-taxanomy
dada1_meta<-cbind(dada1,meta1)
write.table(dada1_meta,paste0(output,"SV_norm_table_Assal.txt"),sep="\t",row.names = TRUE,quote = FALSE)



