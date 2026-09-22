#Author: phase-2 extension (Rishi Gupta / Claude), modeled on Farnaz Fouladi's Afshar_16S_TaxaClassification.R
#Date: 2026
#Description: This script generates taxonomic tables
#Cohort: PRJNA1480144, IL (ileostomy) arm only. Pre/Post, 44 subjects (21 paired; the rest
#Pre-only — see CLAUDE.md / README for the attrition caveat). SD (SDT) arm is a separate
#cohort (SDT_16S_TaxaClassification.R), not merged here — one study, one procedure.
#Note: unlike Afshar's original script, our metaData.txt already has clean time/ID columns
#(no title-string parsing needed).

pipeRoot <- dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
path <- paste0(pipeRoot,"/input/RYGB_Ileostomy/")
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

meta1<-meta
meta1$time<-relevel(factor(meta1$time),ref="Pre")
meta1$prepost<-sapply(meta1$time,function(x){if (x=="Pre") return(0) else return(1)})

for(t in c("Phylum","Class","Order","Family","Genus")){

  t1<-getTaxaTable(dada,taxa,t)
  t1_norm<-norm(t1)
  t1_normMeta<-cbind(t1_norm,meta1)
  write.table(t1_normMeta,paste0(output,t,"_norm_table_Ileostomy.txt"),sep = "\t",row.names = TRUE,quote = FALSE)
}

#SV table
dada1<-norm(dada)
num<-c(1:nrow(taxa))
taxanomy<-apply(taxa,1,function(x){paste0(x[1],"_",x[2],"_",x[3],"_",x[4],"_",x[5],"_",x[6])})
taxanomy<-paste0(taxanomy,"_",num)
colnames(dada1)<-taxanomy
dada1_meta<-cbind(dada1,meta1)
write.table(dada1_meta,paste0(output,"SV_norm_table_Ileostomy.txt"),sep="\t",row.names = TRUE,quote = FALSE)
