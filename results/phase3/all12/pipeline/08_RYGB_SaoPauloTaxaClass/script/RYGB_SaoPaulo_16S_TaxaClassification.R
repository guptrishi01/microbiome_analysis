#Author: phase-3 extension (Rishi Gupta / Claude), modeled on the Cholecystectomy/Ilhan template
#Date: 2026
#Description: This script generates taxonomic tables
#Cohort: PRJNA1364303 (Sao Paulo, RYGB +/- exercise), T1/T2/T3, 31 subjects with >=2
#timepoints (of 32 total RYGB-arm subjects; filtered to paired at metadata-construction
#time, same convention as IleocecalResection). Mixed with a mouse arm and single-timepoint
#LEAN/CTRL controls in the same BioProject -- excluded at metadata-construction time via
#sample_alias pattern matching (RYGB/RYGB+ET only), not here.

pipeRoot <- dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
path <- paste0(pipeRoot,"/input/RYGB_SaoPaulo/")
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
meta$Group<-factor(meta$Group,levels = c("T1","T2","T3"))
meta$prepost<-sapply(meta$Group,function(x){if (x=="T1") return(0) else return(1)})

for(t in c("Phylum","Class","Order","Family","Genus")){

  t1<-getTaxaTable(dada,taxa,t)
  t1_norm<-norm(t1)
  #norm() drops any sample with <=1000 total reads -- meta must be re-subset to match per
  #rank, not assumed to always align (see CLAUDE.md's row-count-mismatch note from phase 2).
  t1_normMeta<-cbind(t1_norm,meta[rownames(t1_norm), , drop=FALSE])
  write.table(t1_normMeta,paste0(output,t,"_norm_table_RYGB_SaoPaulo.txt"),sep = "\t",row.names = TRUE,quote = FALSE)
}

#SV table
dada1<-norm(dada)
num<-c(1:nrow(taxa))
taxanomy<-apply(taxa,1,function(x){paste0(x[1],"_",x[2],"_",x[3],"_",x[4],"_",x[5],"_",x[6])})
taxanomy<-paste0(taxanomy,"_",num)
colnames(dada1)<-taxanomy
dada1_meta<-cbind(dada1,meta[rownames(dada1), , drop=FALSE])
write.table(dada1_meta,paste0(output,"SV_norm_table_RYGB_SaoPaulo.txt"),sep="\t",row.names = TRUE,quote = FALSE)
