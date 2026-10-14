#Author: Farnaz Fouladi
#Date: 10-01-2020
#Description: This script generates taxonomic tables

#BioLockJ configuration: Ali Sorgen
#Date: 10-06-2020

pipeRoot <- dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
path <- paste0(pipeRoot,"/input/RYGB_Ilhan2020/")
output <- file.path(dirname(getwd()),"output/")

metaData<-paste0(path,"Metadata/metaData_Ilhan.txt")
funcScript <- paste0(moduleDir,"/resources/functions.R")
source(funcScript)

dada<-read.table(paste0(path,"DADA2/ForwardReads.txt"),sep="\t",header=TRUE)
taxa<-read.table(paste0(path,"DADA2/taxForwardReads.txt"),sep="\t",header=TRUE)
meta<-read.table(metaData,sep="\t",header=TRUE,row.names = 1)

#Modify Escherichia name
taxa$Genus = as.factor(taxa$Genus)
typoIndex = which(levels(taxa$Genus)=="Esherichica/Shigella")
levels(taxa$Genus)[typoIndex]="Escherichia/Shigella"

dada<-dada[match(rownames(meta),rownames(dada)),]

meta1<-meta[meta$env_material=="fecal" & (meta$Group=="Baseline" |meta$Group=="6M" | meta$Group=="12M"), ]
dada<-dada[meta$env_material=="fecal" & (meta$Group=="Baseline" |meta$Group=="6M" | meta$Group=="12M"), ]
meta1$Group<-factor(meta1$Group,levels = c("Baseline","6M","12M"))
meta1$prepost<-sapply(meta1$Group,function(x){if (x=="Baseline") return(0) else return(1)})

for(t in c("Phylum","Class","Order","Family","Genus")){
  
  t1<-getTaxaTable(dada,taxa,t)
  t1_norm<-norm(t1)
  t1_normMeta<-cbind(t1_norm,meta1)
  write.table(t1_normMeta,paste0(output,t,"_norm_table_Ilhan.txt"),sep = "\t",row.names = TRUE,quote = FALSE)
}

#SV table
dada1<-norm(dada)
num<-c(1:nrow(taxa))
taxanomy<-apply(taxa,1,function(x){paste0(x[1],"_",x[2],"_",x[3],"_",x[4],"_",x[5],"_",x[6])})
taxanomy<-paste0(taxanomy,"_",num)
colnames(dada1)<-taxanomy
dada1_meta<-cbind(dada1,meta1)
write.table(dada1_meta,paste0(output,"SV_norm_table_Ilhan.txt"),sep="\t",row.names = TRUE,quote = FALSE)

