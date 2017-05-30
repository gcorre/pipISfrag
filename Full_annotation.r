#### Script to merge all annotations from VISA pipeline
#### GC



require(plyr);require(dplyr);require(reshape) ;require(stringr)


# Retrieve arguments from command line ------------------------------------

# 

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Usage : Rscript script.r pwd", call.=FALSE)
}


print(args[1])


setwd(args[1])



# Process annotation tables -----------------------------------------------

FILES <- list.files(path = "./10-Annotation/RefSeq/",full.names = T)

Sample_name <- str_replace(basename(grep(FILES,pattern = "distTSS",value = T)),pattern = "_distTSS_IS.txt",replacement = "")

All_annot <- read.delim(grep(FILES,pattern = "_IntragenicIS.txt",value = T),header = F)
names(All_annot) <- c("Chr_IS","Start_IS","Stop_IS","Name_IS","Score_IS","Strand_IS","ReadCount_IS","Database", "Chr_Gene","Start_Gene","Stop_Gene","RefSeq","Score_Gene","Strand_Gene")

refFlat <- read.delim("../../dataset/hg19/2017-04-03/RefSeq/refFlat.txt", header = F, col.names = c("GeneSymbol","RefSeq","Chr_Gene","Strand","TSS_Start","TSS_Stop","CDS_Start","CDS_Stop","ExonCount","Exon_Start","Exon_Stop"))
# 
# HGNC <- read.delim("./dataset/hg19/2017-04-03/HGNC/HGNC.txt", header = T) %>% filter(Status == "Approved")

exons <- read.delim(grep(FILES,pattern = "_IntronExonIS.txt",value = T), header = F)
names(exons) <- c("Chr_IS","Start_IS","Stop_IS","Name_IS","Score_IS","Strand_IS","ReadCount_IS", "Chr_Exon","Start_Exon","Stop_Exon","Name_Exon","GeneSymbol","Strand_Gene")

quantitative <- read.delim(paste("./09-quantitativeIS/",Sample_name,"_quantIS.bed",sep = ""),header = F)
names(quantitative) <- c("Chr_IS","Start_IS","Stop_IS","Name_IS","Score_IS","Strand_IS","ReadCount_ISquant","Fragment_Count")

# Add gene symbol -------------------------------------------------

All_annot <- All_annot %>% 
  left_join(refFlat %>% select(1,2,3), by =c("RefSeq","Chr_Gene")) 

## collapse IS per gene symbol and classify as intra/inter/promoter
All_annot_collapsed <- All_annot %>%
  group_by(Chr_IS,Start_IS,Stop_IS,Name_IS,Score_IS,Strand_IS,ReadCount_IS, GeneSymbol) %>%
  summarise(Database = toString(basename(as.character(unique(Database)))), refseq = toString(RefSeq))%>%
  mutate(Regions=ifelse(str_detect(Database,"refFlat.bed"),yes = "Intragenic",no = ifelse(str_detect(Database,pattern = "Prom.bed"),yes = "Promoter",no = "Intergenic"))) %>% 
  ungroup()

## Add exon information to intragenic IS
Exons_collapsed <- exons %>% 
  group_by(Chr_IS,Start_IS,Stop_IS,Name_IS,Score_IS,Strand_IS,ReadCount_IS, GeneSymbol) %>% 
  summarise(Exons = toString(Name_Exon)) %>% mutate(Regions = "Intragenic")

All_annot_collapsed <- left_join(All_annot_collapsed,Exons_collapsed,by = c("Chr_IS","Start_IS","Stop_IS","Name_IS","Score_IS","Strand_IS","ReadCount_IS","GeneSymbol","Regions")) %>% as.data.frame()


All_annot_collapsed$Position <-ifelse((All_annot_collapsed$Regions=="Intragenic") & is.na(All_annot_collapsed$Exon),
                                      yes = "Intronic",
                                      no = ifelse((All_annot_collapsed$Regions=="Intragenic") & !is.na(All_annot_collapsed$Exon),
                                                  yes = "Exonic",
                                                  no = All_annot_collapsed$Regions))


All_annot_collapsed = full_join(All_annot_collapsed,quantitative %>% select(1,2,3,6,7,8),by = c("Chr_IS", "Start_IS", "Stop_IS", "Strand_IS")) %>% mutate(Name_IS = paste(Chr_IS,Start_IS,Strand_IS,sep = "_"))

write.table(All_annot_collapsed,paste("./10-Annotation/RefSeq/",Sample_name,"_Full_annotation.txt",sep = ""),quote = F,row.names = F,sep="\t")

