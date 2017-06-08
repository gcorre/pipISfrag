#Master script for R analysis

library(plyr,quietly = T,verbose=F,warn.conflicts=F)
library(dplyr,quietly = T,verbose=F,warn.conflicts=F)
library(reshape,quietly = T,verbose=F,warn.conflicts=F)
library(stringr,quietly = T,verbose=F,warn.conflicts=F)
library(colorRamps,quietly = T,verbose=F,warn.conflicts=F)
library(ggplot2,quietly = T,verbose=F,warn.conflicts=F)
library(gridExtra,quietly = T,verbose=F,warn.conflicts=F)
require(knitr)
require(scales) 
require(grid)
require(gridExtra)
require(tidyr)
require(wordcloud)
require(tm)
require(colorRamps)

# Retrieve arguments from command line ------------------------------------

# 

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Usage : Rscript script.r pwd", call.=FALSE)
}

print(args[1])

setwd(args[1])

# perfom the full annotation with intron/exon
source("/home/tempGC/Scripts/Full_annotation.r")

# perfom epigenetic annotation
source("/home/tempGC/Scripts/Epigenetic_ROC.R")

# Make a report
source("/home/tempGC/Scripts/Report.r")