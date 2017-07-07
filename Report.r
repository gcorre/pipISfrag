### Guillaume Corre
### Genethon
### Analysis of the VISA pipeline annotation
############################################
############################################
############################################


# Retrieve arguments from command line ------------------------------------

TAG <- str_extract(Sample_name,pattern = "TAG[0-9]+")
FILE <- list.files("10-Annotation/RefSeq/",pattern = "Full_annotation.txt",full.names = T)
full_annot <- read.delim(FILE,header=T)

FILE <- list.files("10-Annotation/RefSeq/",pattern = "_IntragenicIS.txt",full.names = T)
All_annot <- read.delim(FILE,header=F)
names(All_annot) <- c("Chr_IS","Start_IS","Stop_IS","Name_IS","Score_IS","Strand_IS","ReadCount_IS","Database", "Chr_Gene","Start_Gene","Stop_Gene","RefSeq","Score_Gene","Strand_Gene")

FILE <- list.files("10-Annotation/RefSeq/",pattern = "FragLength",full.names = T)
Corrected_FragLengthCount <- read.delim(FILE,header=T)


# Demultiplexing ----------------------------------------------------------
demux <- read.delim("./log/demultiplex-merged.log", sep ="\t",col.names = c("TAG","Count","File R1","File R2"),header=T)


# Prepare output table ----------------------------------------------------

sumtb <- data.frame("Sample"=Sample_name,
                    "UniqueIS" = full_annot %>% 
                      distinct(Chr_IS,Start_IS,Strand_IS)%>% 
                      nrow(),
                    "AnnotatedIS" = full_annot %>% nrow(),
                    "Targeted Genes" = full_annot %>% 
                      filter(str_detect(Regions,pattern = "Intragenic"))%>%
                      summarise("Targeted Genes"= n_distinct(GeneSymbol)),
                    "Reads" = full_annot %>% 
                      select(1:7)%>% 
                      distinct() %>% 
                      summarise("Reads count"=sum(ReadCount_IS)))

sumtb <- sumtb %>% 
  bind_cols(as.data.frame(t(as.matrix(prop.table(table(full_annot$Position))*100))))

sumtb <- format.data.frame(sumtb,big.mark = ',',big.interval=3L,dec=".",digits = 2, justify = "centre")


# Statistics fastq, bed, sam files --------------------------------------------------

## stat fastq
statFQ <- list.files(pattern = ".fqstat",full.names = T,recursive = T)

all_stat <- list()
for(file in statFQ){
  name <- str_replace(file,pattern = ".fqstat","")
  all_stat[[name]] <- read.delim(file,header = F,sep = "\t",row.names = NULL)
}

fastq_stat = ldply(all_stat,.id = "Sample") %>% 
  mutate(V1 = ifelse(str_detect(V1,"No reads in"),yes = "reads",as.character(V1))) %>% 
  cast(Sample ~ V1, value = "V2",fill = 0)

fqReadCount <- fastq_stat %>% select(Sample,reads)

## bed files length
bedfiles <- list.files(pattern = "*[0-9].+.bed$",recursive = T,full.names = T)
bedfiles <- grep(bedfiles,pattern = "10-random",invert = T,value = T)

all_BED <- list()
for(bedfile in bedfiles){
  name <- bedfile
  all_BED[[name]] <- data.frame(reads=length(scan(bedfile,what = "character",sep = "\n")))
}

bedISCount <- ldply(all_BED,.id = "Sample")


## sam files length

samfiles <- list.files(pattern = ".sam.samstat$",recursive = T,full.names = T)

all_sam <- list()
for(samfile in samfiles){
  name <- samfile
  temp <- scan(samfile,what = "character",sep = "\n")
  all_sam[[name]] <- data.frame(reads= as.numeric(str_extract(temp[grep(pattern = "total",x = temp)],pattern = "^[0-9]+"))- as.numeric(str_extract(temp[grep(pattern = "read1",x = temp)],pattern = "^[0-9]+")))
}
samISCount <- ldply(all_sam,.id = "Sample")



## Merge reads count table
all_counts <- bind_rows(fqReadCount,bedISCount,samISCount)


## annotate and compute percentages

correspondance <- read.delim("/home/tempGC/Ressources/correspondance_pipeline.txt",stringsAsFactors = F,allowEscapes = T) %>% select(-1)
correspondance$pattern <- str_replace(string=correspondance$pattern,pattern   = "TAG",replacement = TAG)


for(i in 1:nrow(correspondance)){
  print(i)
  x <- grep(all_counts$Sample,pattern = correspondance$pattern[i],value = T)
  if(length(x)>0){
  correspondance$"Sample"[i] <- x
  } else{
    correspondance$"Sample"[i] <- ""
  }
}

correspondance <- correspondance %>% select(Sample,description,parent,Step) %>% filter(Sample !="")

all_counts <- all_counts %>% left_join(correspondance, by = "Sample") %>% filter(!is.na(description)) %>% mutate(reads = as.numeric(reads))


for(i in 1: nrow(all_counts)){
  
  all_counts$pct_parent[i] <- round(all_counts$read[i] / sum(all_counts$reads[which(str_trim(all_counts$description)%in%unlist(str_split(all_counts$parent[i],pattern = ";")))]) * 100,digits = 1)
  all_counts$pct_total[i] <- round(all_counts$read[i] / all_counts$reads[which(all_counts$description=="R1_raw")] * 100,digits = 1)
}



##### draw pipeline numbers



patterns_reads <- paste("s/>",
                        all_counts$description,"</>",
                        paste(format(all_counts$reads,big.mark = ",",trim = T),paste(format(all_counts$pct_parent),"%",sep=""),paste(format(all_counts$pct_total),"%",sep=""),
                              sep="--"),
                        "</g",
                        sep="")


patterns <- c(patterns_reads, paste("s/TAG/",TAG,"/g",sep = ""),paste("s/VISASAMPLE/",Sample_name,"/g",sep = ""))
write.table(patterns,paste("QC_",Sample_name,"/pattern.txt",sep = ""),quote = F,row.names = F,col.names = F)



system(paste("sed -f ",paste("QC_",Sample_name,"/pattern.txt",sep = "")," /home/tempGC/ressources/empty_pipeline_reads.svg > pipeline_annotated.svg "))
       
system("rsvg-convert -f png -o pipeline_annotated.png pipeline_annotated.svg")



# Reporting ---------------------------------------------------------------

rmarkdown::render(input = "/home/tempGC/Scripts/squeleton_VISA_nrLMPCR.Rmd", output_format = "html_document" , output_dir = args[1], output_file = paste(Sample_name,"_report.html",sep = ""))

rm(list = grep(ls(),pattern = "sumtb|Corrected_FragLengthCount|Sample_name|full_annot|all_counts",invert = T,value = T))

save.image(file = paste("./10-Annotation/RefSeq/",Sample_name,".RData",sep = ""))