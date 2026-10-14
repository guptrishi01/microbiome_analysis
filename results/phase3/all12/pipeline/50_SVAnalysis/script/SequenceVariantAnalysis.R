#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: Compare BS and Assal 16S datasets at the sequence variant level.

#BioLockJ configuration: Ali Sorgen
#Date: 09-30-2020
# rm(list=ls())

#Libraries
library(ggplot2)
library(gridExtra)
library(ggsignif)
library(ggrepel)
library(stringr)

pipeRoot = dirname(dirname(getwd()))

# input<-"./input/"
input <- paste0(pipeRoot,"/input/")
output = file.path(dirname(getwd()),"output/")
# output<-"./output/"
taxa<-"Seq"
moduleDir <- dirname(getwd())
funcScript <- paste0(moduleDir,"/resources/functions.R")
source(funcScript)
# source("./Rcode/RYGB_IntegratedAnalysis/functions.R")

studies<-c("BS-p1M","BS-p6M","Assal-p3M","Assal-p1Y","Assal-p2Y")
studyNames<-c("BS-1 month","BS-6 months","Assal-3 months","Assal-1 year","Assal-2 years")

#Compare BS-1 month versus Assal-3 months
s=1
s1=3
path<-file.path(pipeRoot,dir(pipeRoot, pattern="CompareStudiesSV"),"input")
# path<-paste0(output,"MixedLinearModels/")
df<-compareStudies(path,taxa,strsplit(studies[s],"-")[[1]][1],strsplit(studies[s1],"-")[[1]][1],strsplit(studies[s],"-")[[1]][2],strsplit(studies[s1],"-")[[1]][2])

taxanomy<-read.table(paste0(input,"RYGB_BS/DADA2/taxForwardReads.txt"),sep="\t",header = TRUE)
tax<-taxanomy[df$bugName,]
df1<-cbind(df,tax)
df1$Genus<-factor(df1$Genus)

df2<-df1[(df1[,"pval1"]<log10(0.05) & df1[,"pval2"]<log10(0.05)) | (df1[,"pval1"]>-log10(0.05) & df1[,"pval2"]>-log10(0.05)),]

theme_set(theme_classic())
xlab="BS- 1 month versus baseline"
ylab="Assal- 3 months versus baseline"

plot<-ggplot(data=df1,aes(x=pval1,y=pval2))+geom_point(size=1,aes(color=Genus))+
  geom_hline(yintercept = 0,linetype="dashed", color = "red")+
  geom_vline(xintercept = 0,linetype="dashed", color = "red")+
  labs(x=xlab,y=ylab)+theme(legend.position = "none")+
  geom_text_repel(data=df2,aes(x=pval1,y=pval2,label=Genus),segment.colour="red",size=2.5,min.segment.length = 0,
                  segment.color="grey",segment.size=0.2)
pdf(paste0(output,"SV_BS_Assal.pdf"),width = 5,height = 5)
print(plot)
dev.off()

