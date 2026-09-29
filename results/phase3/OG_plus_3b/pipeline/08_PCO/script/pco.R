#Author: Farnaz Fouladi
#Date: 08-17-2020
#Description: Ordination plot.

#BioLockJ configuration: Ali Sorgen
#Date: 09-30-2020
# rm(list=ls())

#Libraries
library(vegan)
library(stringr)
set.seed(123)

pipeRoot <- dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
outputDir <- file.path(moduleDir,"output/")
inputModule <- dir(pipeRoot, "CombineCountTableslog10", full.names = TRUE)
input <- file.path(inputModule,"output")
message("Taking input files from folder: ", input)

taxaLevel <- "Genus"
colors <- c("red","blue","darkgreen","darkorange2","purple","hotpink","black","firebrick")
funcScript <- file.path(file.path(moduleDir, "resources"),"functions.R")
message("Using resource file: ", funcScript)
source(funcScript)

inFile <- file.path(input, paste0(taxaLevel, "_countTable_merged_log10.txt"))
message("Reading input table: ", inFile)
myT <- read.table(inFile,sep="\t",header = TRUE)

metaFile <- file.path(input,"metaData_merged.txt")
message("Reading metadata from file: ", metaFile)
meta <- read.table(metaFile,sep="\t",header = TRUE)

meta$TimeEachStudy <- paste0(meta$Study,"_",meta$time)

outFile <- file.path(outputDir, paste0(taxaLevel,"_PCO.pdf"))
message("Saving output to file: ", outFile)
pdf(outFile, width = 5,height = 10)
par(mfrow=c(2,1))

getPCO(myT,meta,"Study",names=levels(factor(meta$Study)),colors)
getPCO(myT,meta,"TimeEachStudy",names=levels(factor(meta$TimeEachStudy)),colors)

dev.off()

#PERMANOVA 
permanova <- adonis(myT~factor(meta$time)*factor(meta$Study))
permaOut <- file.path(outputDir, paste0(taxaLevel, "_permanovaResults.txt"))
message("Saving permanova summary to file: ", permaOut)
capture.output(permanova, file = permaOut)
# "Call:
# adonis(formula = myT ~ factor(meta$time) * factor(meta$Study)) 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
#                                       Df SumsOfSqs MeanSqs F.Model
# factor(meta$time)                      1    0.7379 0.73786  9.1916
# factor(meta$Study)                     3    6.7607 2.25357 28.0729
# factor(meta$time):factor(meta$Study)   3    0.2781 0.09269  1.1546
# Residuals                            180   14.4496 0.08028        
# Total                                187   22.2262                
#                                           R2 Pr(>F)    
# factor(meta$time)                    0.03320  0.001 ***
# factor(meta$Study)                   0.30418  0.001 ***
# factor(meta$time):factor(meta$Study) 0.01251  0.219    
# Residuals                            0.65012           
# Total                                1.00000           
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"

message("")
message("All done!")
message("")

sessionInfo()
