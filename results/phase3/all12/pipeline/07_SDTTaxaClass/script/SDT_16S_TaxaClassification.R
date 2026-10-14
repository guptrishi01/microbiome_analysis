#Author: phase-2 extension (Rishi Gupta / Claude), modeled on Farnaz Fouladi's Afshar_16S_TaxaClassification.R
#Date: 2026
#Description: This script generates taxonomic tables
#Cohort: PRJNA1480144, SD (SDT) arm only. Pre/Post, 34 subjects (11 paired; the rest
#Pre-only — see CLAUDE.md / README for the attrition caveat). IL (ileostomy) arm is a
#separate cohort (Ileostomy_16S_TaxaClassification.R), not merged here.

pipeRoot <- dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
path <- paste0(pipeRoot,"/input/RYGB_SDT/")
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
  #norm() drops any sample with <=1000 total reads -- meta1 must be re-subset to match per
  #rank rather than assumed to always align.
  t1_normMeta<-cbind(t1_norm,meta1[rownames(t1_norm), , drop=FALSE])
  write.table(t1_normMeta,paste0(output,t,"_norm_table_SDT.txt"),sep = "\t",row.names = TRUE,quote = FALSE)
}

#SV table
dada1<-norm(dada)
num<-c(1:nrow(taxa))
taxanomy<-apply(taxa,1,function(x){paste0(x[1],"_",x[2],"_",x[3],"_",x[4],"_",x[5],"_",x[6])})
taxanomy<-paste0(taxanomy,"_",num)
colnames(dada1)<-taxanomy
dada1_meta<-cbind(dada1,meta1[rownames(dada1), , drop=FALSE])
write.table(dada1_meta,paste0(output,"SV_norm_table_SDT.txt"),sep="\t",row.names = TRUE,quote = FALSE)
