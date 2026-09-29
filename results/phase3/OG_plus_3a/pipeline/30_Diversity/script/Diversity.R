#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: Compare divesity between 16S datasets.

#BioLockJ configuration: Ali Sorgen
#Date: 09-30-2020
rm(list=ls())

#Libraries
library(ggplot2)
library(nlme)
library(stringr)
library(vegan)

pipeRoot = dirname(dirname(getwd()))
output = file.path(dirname(getwd()),"output/")
# output<-"./output/"
input <- paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "CombineCountTableslog10"),"/output/")

taxa="Genus"

##Divesrity
myT1<-read.table(paste0(input,taxa,"_countTable_merged_log10.txt"),sep="\t",header = TRUE)
# myT1<-read.table(paste0(output,taxa,"_countTable_merged_log10.txt"),sep="\t",header = TRUE)
meta<-read.table(paste0(input,"metaData_merged.txt"),sep="\t",header = TRUE)

#Timepoints
meta$timepoint<-sapply(meta$timepoint,function(x){
  if (x=="0" | x=="Baseline" | x=="Pre") return("0")
  else if (x=="1") return("1")
  else if (x=="3M") return("3")
  else if (x=="6" | x=="6M" | x=="Post") return("6")
  else if (x=="12M" | x=="1Y") return("12")
  else if (x=="2Y") return("24")
  else if (x=='3Y') return("36")
})

meta$timepoint<-factor(meta$timepoint,levels = c("0","1","3","6","12","24"))

#Normalized count
myT2<-10^(myT1)-1
myT2$div<-diversity(myT2,index = "shannon")
myT2_c<-cbind(myT2,meta)

# pdf(paste0(output,"Diversity.pdf"),width = 5, height = 5, useDingbats = FALSE)
plot2<-ggplot(data=myT2_c,aes(x=timepoint,y=div))+geom_boxplot(outlier.shape = NA)+facet_wrap(~Study)+
  geom_jitter(position = position_jitter(width = 0.2),size=0.2)+labs(y="Shannon Divesrity Index",x="Time points (months)")
# print(plot2)
# dev.off()

ggsave(plot2, file=paste0(output,"Diversity.pdf"), width = 5, height = 5)

#ANOVA test: No significant result
s="Ilhan"
myT2_sub<-myT2_c[myT2_c$Study==s,]
anova1 <- anova(lm(myT2_sub$div~myT2_sub$timepoint))
anova2 <- anova(lme(div ~ timepoint, random = ~ 1 | ID, data=myT2_sub))

capture.output(anova1, anova2, file = paste0(output,taxa,"_anovaResults.txt"))
