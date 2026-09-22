#Reference:Wirbel J, Pyl PT, Kartal E, et al. Meta-analysis of fecal metagenomes
#reveals global microbial signatures that are specific for colorectal cancer. 
#Nat Med 2019;25:679-689.

#Description: Figures for Machine Learning 


#BioLockJ configuration: Ali Sorgen
#Date: 10-01-2020
# rm(list=ls())

# Libraries
library("tidyverse")
library("cowplot")
library("pROC")
library("yaml")
library("stringr")

pipeRoot = dirname(dirname(getwd()))
message(paste0("pipeRoot = ", pipeRoot))

output = file.path(dirname(getwd()),"output/")
message(paste0("output = ", output))

# output<-"./output/"
taxa<-"Genus"

# ml.method="randomForest"
# ml.method="lasso"
args <- commandArgs(trailingOnly = TRUE)
ml.method <- args[1]

metaPath <- paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "CombineCountTablesrelab"),"/output/")

predictPath <- paste0(pipeRoot,"/",str_subset(dir(pipeRoot), paste0("TrainModel",ml.method)),"/output/")

meta<-read.table(paste0(metaPath,"metaData_merged.txt"),sep="\t",header = TRUE)
# meta<-read.table(paste0(output,"metaData_merged.txt"),sep="\t",header = TRUE)
studies <- meta %>% pull(Study) %>% unique
# pred.matrix <- read.table(paste0(output,ml.method,"/predictions_",ml.method,".tsv"), sep='\t', check.names = FALSE)
pred.matrix <- read.table(paste0(predictPath,ml.method,"/predictions_",ml.method,".tsv"), sep='\t', check.names = FALSE)
pred.matrix$Sample_ID<-rownames(pred.matrix) 
pred.matrix <- as_tibble(pred.matrix)
df.all <- inner_join(meta, pred.matrix, by='Sample_ID')

# ##############################################################################
# Calculate AUROCs
auroc.all <- tibble()

for (study.train in studies){
  for (study.test in studies){
    predictor <- df.all %>%
      filter(Study == study.test) %>% 
      pull(study.train) 
     response <- df.all %>%
      filter(Study == study.test) %>% 
      pull(time)                  
    temp <- roc(predictor=predictor, response = response, ci=TRUE)
    
    auroc.all <- bind_rows(auroc.all, 
                           tibble(study.train=study.train, 
                                  study.test=study.test,
                                  AUC=c(temp$auc)))
    
  }
}

# ##############################################################################
# AUROC heatmap

col.scheme.heatmap <- c("black","slateblue4","slateblue2","skyblue1","darkslategray1")

plot.levels <- levels(factor(auroc.all$study.test))
studies<-levels(factor(auroc.all$study.test))
theme_set(theme_bw(base_size = 14))
g <- auroc.all %>% 
  mutate(study.test=factor(study.test, levels=plot.levels)) %>% 
  mutate(study.train=factor(study.train, levels=rev(plot.levels))) %>% 
  mutate(CV=study.train == study.test) %>%
  ggplot(aes(y=study.train, x=study.test, fill=AUC,size=CV ,color=CV)) +
  geom_tile() + theme_bw(base_size = 16) +
  # test in tiles
  geom_text(aes_string(label="format(AUC, digits=2)"), col='white', size=7)+
  # color scheme
  scale_fill_gradientn(colours = col.scheme.heatmap, limits=c(0.5, 1)) +
  # axis position/remove boxes/ticks/facet background/etc.
  scale_x_discrete(position='top') + 
  theme(axis.line=element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x.top = element_text(angle=45, hjust=.1), 
        panel.grid=element_blank(), 
        panel.border=element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_blank()) + 
  xlab('Test Set') + ylab('Training Set') + 
  scale_color_manual(values=c('grey', 'yellow'), guide=FALSE) + 
  scale_size_manual(values=c(0, 1), guide=FALSE)

# model average
g2 <- auroc.all %>% 
  filter(study.test != study.train) %>% 
  group_by(study.train) %>% 
  summarise(AUROC=mean(AUC)) %>% 
  mutate(study.train=factor(study.train, levels=rev(plot.levels))) %>% 
  ggplot(aes(y=study.train, x=1, fill=AUROC)) + 
  geom_tile() + theme_bw(base_size = 16) +
  geom_text(aes_string(label="format(AUROC, digits=2)"), col='white', size=7)+
  scale_fill_gradientn(colours = col.scheme.heatmap, limits=c(0.5, 1), 
                       guide=FALSE) + 
  scale_x_discrete(position='top') + 
  theme(axis.line=element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid=element_blank(), 
        panel.border=element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_blank()) + 
  xlab('Model Average') + ylab('')



pdf(paste0(output,'/performance_heatmap_', ml.method,'.pdf'), 
    width = 12, height = 7.5, useDingbats = FALSE)
plot_grid(g, g2, rel_widths = c(5/6, 2/6), align = 'h')
dev.off()

# ##############################################################################
# LOSO bars
df.plot.loso <- auroc.all %>% 
  filter(study.test != study.train) %>% 
  select(study.test, AUC) %>% 
  mutate(type='Single Study')

for (study in studies){
  predictor <- df.all %>%
    filter(Study == study) %>% 
    #filter(time=='post') %>% 
    pull(LOSO) 
  response <- df.all %>%
    filter(Study == study) %>% 
    #filter(time=='pre') %>% 
    pull(time)                  
  temp <- roc(predictor=predictor, response = response, ci=TRUE)
  
  df.plot.loso <- bind_rows(
    df.plot.loso, tibble(study.test=study, AUC=c(temp$auc), 
                         type='LOSO'))
  
}

temp <- df.plot.loso %>% 
  group_by(study.test, type) %>% 
  summarise(AUROC=mean(AUC), sd=sd(AUC)) %>% 
  ungroup() %>% 
  mutate(study.test = factor(study.test, levels = plot.levels)) %>% 
  mutate(type=factor(type,levels = c("Single Study","LOSO")))

g3 <- df.plot.loso %>% 
  mutate(study.test = factor(study.test, levels = plot.levels)) %>% 
  mutate(type=factor(type, levels=c('Single Study', 'LOSO'))) %>% 
  ggplot(aes(x=study.test, y=AUC, fill=type)) + 
  geom_bar(stat='summary', position=position_dodge(), 
           size=1, colour='black') + 
  geom_point(position = position_dodge(width = 0.9)) + 
  scale_y_continuous(breaks=seq(0, 1.1, by=0.1)) + 
  ylab('AUROC') + xlab('Study test') + 
  coord_cartesian(ylim=c(0, 1.05), expand=FALSE) + 
  theme(panel.grid.major.x=element_blank()) + 
  theme_classic() + 
  scale_fill_manual(values=c('white', 'darkgrey'), 
                    labels=c('Dataset average', 'LOSO'), name='')+
  geom_errorbar(data=temp, inherit.aes = FALSE, 
                aes(x=study.test, ymin=AUROC-sd, ymax=AUROC+sd,fill=type), 
                position = position_dodge(width = 0.9), width=0.2)




ggsave(g3, filename = paste0(output,'/loso_performance_', 
                             ml.method,'.pdf'),
       width = 5.5, height = 3)


# #######################
# End of script
# #######################












