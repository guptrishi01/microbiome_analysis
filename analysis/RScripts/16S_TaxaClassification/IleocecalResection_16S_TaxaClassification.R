#Author: phase-2 extension (Rishi Gupta / Claude), modeled on Farnaz Fouladi's BS_16S_TaxaClassification.R
#Date: 2026
#Description: This script generates taxonomic tables
#Cohort: PRJNA514452 (ileocecal resection), 0/1/3/6 months, 54 subjects.
#Fecal-only: cohorts/IleocecalResection/metaData.txt already excludes the 320 biopsy-specimen
#(Bp/iRs/niRs) runs at construction time, so no env_material filter is needed here (same
#simplicity as BS's own template, whose metadata is likewise pre-filtered to fecal).

pipeRoot <- dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
path <- paste0(pipeRoot,"/input/RYGB_IleocecalResection/")
output <- file.path(dirname(getwd()),"output/")

metaData<-paste0(path,"Metadata/metaData.txt")
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
meta$Timepoint<-factor(meta$Timepoint,levels = c(0,1,3,6))

for(t in c("Phylum","Class","Order","Family","Genus")){

  t1<-getTaxaTable(dada,taxa,t)
  t1_norm<-norm(t1)
  t1_normMeta<-cbind(t1_norm,meta)
  write.table(t1_normMeta,paste0(output,t,"_norm_table_IleocecalResection.txt"),sep = "\t",row.names = TRUE,quote = FALSE)
}

#SV table
dada1<-norm(dada)
num<-c(1:nrow(taxa))
taxanomy<-apply(taxa,1,function(x){paste0(x[1],"_",x[2],"_",x[3],"_",x[4],"_",x[5],"_",x[6])})
taxanomy<-paste0(taxanomy,"_",num)
colnames(dada1)<-taxanomy
dada1_meta<-cbind(dada1,meta)
write.table(dada1_meta,paste0(output,"SV_norm_table_IleocecalResection.txt"),sep="\t",row.names = TRUE,quote = FALSE)
