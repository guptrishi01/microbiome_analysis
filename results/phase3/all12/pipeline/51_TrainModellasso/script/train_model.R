#Reference:Wirbel J, Pyl PT, Kartal E, et al. Meta-analysis of fecal metagenomes
#reveals global microbial signatures that are specific for colorectal cancer. 
#Nat Med 2019;25:679-689.

#Description: Machine Learning to Find Microbial Signature Across Cohorts

#BioLockJ configuration: Ali Sorgen
#Date: 10-01-2020
# rm(list=ls())

#libraries
library("tidyverse")
library("SIAMCAT")
library("yaml")
library("stringr")


start.time <- proc.time()[1]
set.seed(200)

pipeRoot = dirname(dirname(getwd()))
input <- paste0(pipeRoot,"/",str_subset(dir(pipeRoot), "CombineCountTablesrelab"),"/output/")
output = file.path(dirname(getwd()),"output/")
# output<-"./output/"
taxa<-"Genus"

#Loading tables
feat.all<-t(read.table(paste0(input,taxa,"_countTable_merged_relab.txt"),header = TRUE,quote=''))
# feat.all<-t(read.table(paste0(output,taxa,"_countTable_merged_relab.txt"),header = TRUE,quote=''))
meta<-read.table(paste0(input,"metaData_merged.txt"),header = TRUE,quote='')
# meta<-read.table(paste0(output,"metaData_merged.txt"),header = TRUE,quote='')
stopifnot(all(meta$Sample_ID %in% colnames(feat.all)))

#Set parameters
norm.method="log.std"
n.p<-list(log.n0=1e-05,sd.min.q=0.1,n.p=2,norm.margin=1)
num.folds=8
num.resample=10
#ml.method="randomForest"
# ml.method="lasso"
args <- commandArgs(trailingOnly = TRUE)
ml.method <- args[1]
dir.create(paste0(output,ml.method,"/"))
modsel.crit=list("pr")
min.nonzero.coeff=1
param.fs.ss=list(thres.fs=800,method.fs="AUC")
studies<-unique(meta$Study)

# Model Building
models <- list()
for (study in studies){
  # single study model
  meta.train <- meta %>%
    filter(Study == study)
  
  feat.train <- feat.all[,as.character(meta.train %>% pull(Sample_ID))]
  
  rownames(meta.train) <- meta.train$Sample_ID
  
  if(sum(colnames(feat.train)==rownames(meta.train))!=nrow(meta.train)) stop("error") 
  
  siamcat <- siamcat(feat=feat.train, meta=meta.train,
                     label = 'time', case=1)
  siamcat <- normalize.features(siamcat, norm.method = norm.method,
                                norm.param = n.p, feature.type = 'original',
                                verbose=3)
  siamcat <- create.data.split(siamcat, num.folds = num.folds,
                               num.resample = num.resample,stratify=TRUE)
  siamcat <- train.model(siamcat,
                         method = ml.method,
                         modsel.crit=modsel.crit,
                         min.nonzero.coeff = min.nonzero.coeff,
                         perform.fs = FALSE)
  siamcat <- make.predictions(siamcat)
  siamcat <- evaluate.predictions(siamcat)
  models[[study]] <- siamcat
  save(siamcat, file=paste0(output,ml.method,'/',study, '_', 
                            ml.method ,'_model','.RData'))
  
  write.table(siamcat@pred_matrix,paste0(output,ml.method,'/',study, '_', 
                            ml.method ,'_model','.txt'),sep="\t")
  
  cat("Successfully trained a single study model for study", study, '\n')
  
  # LOSO models
  meta.train <- meta %>%
    filter(Study != study)
  
  feat.train <- feat.all[,as.character(meta.train %>% pull(Sample_ID))]
  
  rownames(meta.train) <- meta.train$Sample_ID
  
  if(sum(colnames(feat.train)==rownames(meta.train))!=nrow(meta.train)) stop("error") 
  
  siamcat <- siamcat(feat=feat.train, meta=meta.train, 
                     label = 'time', case=1)
  siamcat <- normalize.features(siamcat, norm.method = norm.method, 
                                norm.param = n.p, feature.type = 'original', 
                                verbose=3)
  siamcat <- create.data.split(siamcat, num.folds = num.folds, 
                               num.resample = num.resample,stratify=FALSE)
  siamcat <- train.model(siamcat,
                         method = ml.method,
                         modsel.crit=modsel.crit,
                         min.nonzero.coeff = min.nonzero.coeff,
                         perform.fs = FALSE,
                         param.fs = param.fs.loso)
  siamcat <- make.predictions(siamcat)
  siamcat <- evaluate.predictions(siamcat)
  models[[paste0(study, '_LOSO')]] <- siamcat
  save(siamcat, file=paste0(output,ml.method,'/',study, '_loso_', 
                            ml.method, '_model','.RData'))
  
  write.table(siamcat@pred_matrix,paste0(output,ml.method,'/',study, '_loso_', 
                            ml.method, '_model','.txt'),sep="\t")
  
  cat("Successfully trained a LOSO model for study", study, '\n')
}

# ##############################################################################
# make Predictions
pred.matrix <- matrix(NA, nrow=nrow(meta), 
                      ncol=length(studies)+1, 
                      dimnames = list(meta$Sample_ID, 
                                      c(as.character(studies), 'LOSO')))

for (study in studies){
  
  # load model
  siamcat <- models[[study]]
  temp <- rowMeans(pred_matrix(siamcat))
  pred.matrix[names(temp), study] <- temp
  
  # predict other studies
  for (study_ext in setdiff(studies, study)){
    
    meta.test <- meta %>%
      filter(Study == study_ext)
    
    feat.test <- feat.all[,as.character(meta.test %>% pull(Sample_ID))]
    
    rownames(meta.test) <- meta.test$Sample_ID
    
    if(sum(colnames(feat.test)==rownames(meta.test))!=nrow(meta.test)) stop("error") 
    
    siamcat.test <- siamcat(feat=feat.test)
 
    siamcat.test <- make.predictions(siamcat, siamcat.holdout = siamcat.test)
    
    write.table(siamcat.test@pred_matrix,paste0(output,ml.method,'/',"trained_",study, '_PredictionFor_', study_ext,"_",
                                           ml.method,'.txt'),sep="\t")
    temp <- rowMeans(pred_matrix(siamcat.test))
    pred.matrix[names(temp), study] <- temp
    
  }
  
}

# ##############################################################################
# make LOSO Predictions
for (study in studies){
  
  # load model
  siamcat <- models[[paste0(study, '_LOSO')]]
  
  meta.test <- meta %>%
    filter(Study == study)
  
  feat.test <- feat.all[,as.character(meta.test %>% pull(Sample_ID))]
  
  rownames(meta.test) <- meta.test$Sample_ID
  
  if(sum(colnames(feat.test)==rownames(meta.test))!=nrow(meta.test)) stop("error") 
  
  siamcat.test <- siamcat(feat=feat.test)
  
  siamcat.test <- make.predictions(siamcat, siamcat.holdout = siamcat.test)
  
  write.table(siamcat.test@pred_matrix,paste0(output,ml.method,'/','LOSO _PredictionFor_', study,"_",
                                              ml.method,'.txt'),sep="\t")
  
  
  temp <- rowMeans(pred_matrix(siamcat.test))
  pred.matrix[names(temp), 'LOSO'] <- temp
}

# ##############################################################################
# save predictions
write.table(pred.matrix, file=paste0(output,ml.method,'/','predictions_', 
                                     ml.method,'.tsv'), 
            quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE)

cat('Successfully build models in',
    proc.time()[1]-start.time, 'second...\n')

# #######################
# End of script
# #######################


  
  
  
  
