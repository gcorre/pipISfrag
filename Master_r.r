#Master script for R analysis

for(package in c("plyr",
                  "dplyr",
                  "reshape",
                  "stringr",
                  "colorRamps",
                  "ggplot2",
                  "grid",
                  "gridExtra",
                  "knitr",
                  "scales",
                  "tidyr",
                  "wordcloud",
                  "tm",
                  "vegan",
                  "sonicLength",
                  "diverse",
                  "ineq",
                  "entropart",
                  "entropy")){
  if(!require(paste(package),character.only = TRUE,quietly = T,warn.conflicts=F)){
    install.packages(paste(package), dependencies = TRUE,repos = "https://cran.univ-paris1.fr/")
    library(paste(package))
  }
}




# Retrieve arguments from command line ------------------------------------

# 



args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Usage : Rscript script.r pwd", call.=FALSE)
}

print(args[1])

setwd(args[1])

info_file <- list.files(pattern = "_Info.txt")
  
Sample_name <- str_replace(grep(scan(info_file,sep = "\n",what = "character"),pattern = "Sample_",value = T),pattern = "Sample_.+:",replacement = "")

Annotation_path <- str_replace(grep(scan(info_file,sep = "\n",what = "character"),pattern = "Anno",value = T),pattern = "Annotation_path:",replacement = "")



cat("#########\nPerfom the full annotation with intron/exon\n########")
source("/home/tempGC/Scripts/Full_annotation.r")

cat("#########\nPerfom epigenetic annotation\n########")
source("/home/tempGC/Scripts/Epigenetic_ROC.R")

cat("#########\nMake a report\n########")
source("/home/tempGC/Scripts/Report.r")