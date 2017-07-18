#### Script to merge all annotations from VISA pipeline
#### GC

# load all IS (qualitative) and performe the full annotation

# Load annotation tables -----------------------------------------------

FILES <- list.files(path = "./10-Annotation/RefSeq/",full.names = T)

# get annotation file with one line per IS and feature (1 line per transcript)
All_annot <- read.delim(grep(FILES,pattern = "_IntragenicIS.txt",value = T),header = F)
names(All_annot) <- c("Chr_IS","Start_IS","Stop_IS","Name_IS","Score_IS","Strand_IS","ReadCount_IS","Database", "Chr_Gene","Start_Gene","Stop_Gene","RefSeq","Score_Gene","Strand_Gene")

# get annotation table
refFlat <- read.delim(paste(Annotation_path,"RefSeq/refFlat.txt",sep = ""), header = F, col.names = c("GeneSymbol","RefSeq","Chr_Gene","Strand","TSS_Start","TSS_Stop","CDS_Start","CDS_Stop","ExonCount","Exon_Start","Exon_Stop"))


# get annotation in exons
exons <- read.delim(grep(FILES,pattern = "_IntronExonIS.txt",value = T), header = F)
names(exons) <- c("Chr_IS","Start_IS","Stop_IS","Name_IS","Score_IS","Strand_IS","ReadCount_IS", "Chr_Exon","Start_Exon","Stop_Exon","Name_Exon","GeneSymbol","Strand_Gene")



# Add gene symbol -------------------------------------------------

All_annot <- All_annot %>% 
  left_join(refFlat %>% select(1,2,3), by =c("RefSeq","Chr_Gene")) 



## collapse transcripts by gene symbol and classify IS as intra/inter/promoter

All_annot_collapsed <- All_annot %>%
  group_by(Chr_IS,Start_IS,Stop_IS,Name_IS,Score_IS,Strand_IS,ReadCount_IS, GeneSymbol) %>%
  summarise(Database = toString(basename(as.character(unique(Database)))), refseq = toString(RefSeq))%>%
  mutate(Regions=ifelse(str_detect(Database,"refFlat.bed"),yes = "Intragenic",no = ifelse(str_detect(Database,pattern = "Prom.bed"),yes = "Promoter",no = "Intergenic"))) %>% 
  ungroup()



## Add exon information to intragenic IS
Exons_collapsed <- exons %>% 
  group_by(Chr_IS,Start_IS,Stop_IS,Name_IS,Score_IS,Strand_IS,ReadCount_IS, GeneSymbol) %>% 
  summarise(Exons = toString(Name_Exon)) %>%
  mutate(Regions = "Intragenic")

All_annot_collapsed <- left_join(All_annot_collapsed,
                                 Exons_collapsed,
                                 by = c("Chr_IS","Start_IS","Stop_IS","Name_IS","Score_IS","Strand_IS","ReadCount_IS","GeneSymbol","Regions")) %>%
  as.data.frame()


All_annot_collapsed$Position <-ifelse((All_annot_collapsed$Regions=="Intragenic") & is.na(All_annot_collapsed$Exon),
                                      yes = "Intronic",
                                      no = ifelse((All_annot_collapsed$Regions=="Intragenic") & !is.na(All_annot_collapsed$Exon),
                                                  yes = "Exonic",
                                                  no = All_annot_collapsed$Regions))


rm(list = grep(ls(),pattern = "All_annot_collapsed|Sample_name|args",invert = T,value = T))

# Filter fragments with low count and INDELS -------------------------------
# load the quantitative output file with fragment length and count
FILE <- list.files("09-quantitativeIS",pattern = "_quantIS_cluster_sorted.bed", full.names = T)
intermediate_fragment_count <- read.delim(FILE,header=F)
names(intermediate_fragment_count) <- c("Chr_IS","Start_IS","Stop_IS","Name_IS","Score_IS","Strand_IS","Size","ReadCount","Cluster")



# collpase IS mapping at the same position and with the same fragment length
intermediate_fragment_count_sum <- intermediate_fragment_count %>% 
  group_by(Chr_IS,Start_IS,Stop_IS,Strand_IS,Size) %>% 
  summarise(ReadCount = sum(ReadCount)) %>% 
  ungroup() %>% 
  unite_("IS",c("Chr_IS","Start_IS","Strand_IS"))


 temp <- intermediate_fragment_count_sum %>% 
   group_by(IS) %>% 
   summarize(ReadCount_Raw_quant=sum(ReadCount),
             Fragments_Raw_quant = n())


## add a filter to remove low count fragments
# aggregate low count reads to fragment size with +/-1nt (probably an indel during sequencing or trimming)

intermediate_fragment_count_sum <- intermediate_fragment_count_sum %>% 
  group_by(IS) %>% 
  arrange(IS,Size) %>% 
  mutate(ratio = ReadCount/(lag(ReadCount))) %>% 
  mutate(ratio = ifelse(is.na(ratio),1,ratio))


intermediate_fragment_count_sum <- intermediate_fragment_count_sum %>% 
  mutate(Size = ifelse((ratio < 0.2 | ReadCount<=10 ),
                       yes = ifelse((lead(Size)-Size)<3 & !is.na(lead(Size)),
                                    yes = lead(Size),
                                    no= ifelse((Size-lag(Size)) <3 & !is.na(lag(Size)),
                                               yes = lag(Size),
                                               no = Size)),
                       no = Size)) 





Corrected_FragLengthCount <- intermediate_fragment_count_sum %>% 
  group_by(IS,Size) %>% 
  summarise(ReadCount = sum(ReadCount)) %>% ungroup()


temp2 <- Corrected_FragLengthCount %>% 
  group_by(IS) %>% 
  summarize(ReadCount_Corrected_quant=sum(ReadCount),
            Fragments_Corrected_quant = n(),
            ReadCount_Corrected_gt10reads_quant = sum(ReadCount[which(ReadCount>10)]),
            Fragments_Corrected_gt10reads_quant = length(which(ReadCount>10))) %>% 
  arrange(desc(Fragments_Corrected_quant))


temp <- left_join(temp,temp2, by = "IS")



# Apply the sonicLength correction ----------------------------------------
SonicAbundance <- estAbund(Corrected_FragLengthCount$IS,Corrected_FragLengthCount$Size)
temp <- temp %>% mutate(SonicAbundance_quanti=SonicAbundance$theta)


# Collapse all tables -----------------------------------------------------


All_annot_collapsed = left_join(All_annot_collapsed %>%
                                  mutate(Name_IS = paste(Chr_IS,Start_IS,Strand_IS,sep = "_")),
                                temp, by =c("Name_IS"="IS")) 

write.table(All_annot_collapsed,paste("./10-Annotation/RefSeq/",Sample_name,"_Full_annotation.txt",sep = ""),quote = F,row.names = F,sep="\t")

write.table(Corrected_FragLengthCount, paste("./10-Annotation/RefSeq/",Sample_name,"_FragLength_ReadCount_corrected.txt",sep = ""),quote = F,row.names = F,sep="\t")

rm(list = grep(ls(),pattern = "Sample_name|args",invert = T,value = T))

