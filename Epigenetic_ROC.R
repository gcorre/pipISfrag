### Analyse epigenetic features

source("/home/tempGC/Scripts/AUC_curve.r")



AUC <- list()

MergedISrnd <- list.files("10-randomIS/",pattern = "slop|MergedQualIS_Rnd.bed",full.names = T)

MergedISrnd_df <- list()
for(file in MergedISrnd){
  
  MergedISrnd_df[[basename(file)]] <- read.delim(file, header=F)
}

MergedISrnd_df <- ldply(MergedISrnd_df,.id = "Window")

MergedISrnd_df <- MergedISrnd_df %>%
  mutate(Window = str_replace(Window,pattern = "MergedQualIS_Rnd.bed",replacement = "slop0kb.bed")) %>%
  mutate(Window = str_extract(str_extract(Window,pattern = '_[A-Za-z0-9]+.bed'),pattern = "[01]+[a-zA-Z]+"))


#################################
histones <- list.files("10-randomIS/",pattern = "b.tab",full.names = T)
histones_df <- list()
for(file in histones){
  
  histones_df[[basename(file)]] <- read.delim(file, header=T)
}

histones_df <- ldply(histones_df,.id = "Window")
names(histones_df) <- str_replace_all(names(histones_df),pattern = "^X[.]+|\\.$",replacement = "")

histones_df<- histones_df %>% 
  melt.data.frame(id.vars = c("Window","chr","start","end"),variable_name = "Feature") 

histones_df$Feature <- sapply(histones_df$Feature,function(x){
  unlist(str_split(x,pattern = "\\."))[3]
})

histones_df <- histones_df %>% mutate(Window = str_extract(str_extract(Window,pattern = '_[A-Za-z0-9]+.tab'),pattern = "[01]+[a-zA-Z]+"))



#################################


Epi_manual <- list.files("10-randomIS/",pattern = "_.+[01kM]+b.bed",full.names = T)
Epi_manual <- grep(Epi_manual,pattern = "(slop)",invert = T,value = T)
Epi_manual_df <- list()
for(file in Epi_manual){
  
  Epi_manual_df[[basename(file)]] <- read.delim(file, header=F,comment.char ="#" )
}

Epi_manual_df <- ldply(Epi_manual_df,.id = "Window") %>% 
  mutate(Feature = str_extract(str_extract(Window,pattern = '_[A-Za-z0-9]+.bed'),pattern = "[A-Za-z]+")) %>%
  mutate(Window = str_extract(str_extract(Window,pattern = '_[A-Za-z0-9]+.bed'),pattern = "[01]+[a-zA-Z]+")) %>% 
  select(c(1:4,17,8))

names(Epi_manual_df) <- names(histones_df)

#################################

FILE <- list.files("10-randomIS/",pattern = "distTSS.bed",full.names = T)
TSS <- read.delim(FILE,header=F)
TSS <- TSS %>% mutate(Feature = "DistTSS") %>% select(1,2,3,Feature,13)
TSS <- data.frame("Window"= "0kb",TSS)

names(TSS) <- names(histones_df)

#################################

Merged_Epi <- bind_rows(histones_df,Epi_manual_df, TSS) %>% 
  group_by(Window,chr,start,end,Feature) %>% 
  summarize(Sum = sum(value,na.rm = T))

Merged_Epi <- left_join(Merged_Epi,MergedISrnd_df, by = c("Window","chr" = "V1", "start"="V2", "end"="V3"))

Merged_Epi$V4 <- ifelse(str_detect(Merged_Epi$V4,pattern = "random"),yes = "random",no = "Experimental")


rm(list = grep(ls(),pattern = "Merged_Epi|Sample_name|Roc_curve|args",invert = T,value = T))

pdf("10-randomIS/AUC_curves.pdf")
AUC <- list()
for(slop in levels(factor(Merged_Epi$Window))){
  auc <- data.frame()
 for(feature in levels(factor(Merged_Epi$Feature))){
   
    sub = Merged_Epi %>% filter(Window == slop, Feature==feature, !is.na(Sum)) %>% as.data.frame()
    if(nrow(sub)==0 | length(levels(factor(sub$V4)))<2){next}
    sub$Sum <- ifelse(sub$Feature == "GCcontent",1-sub$Sum,sub$Sum) #reverse because AT% comes first in the bedtools nuc output
    Roc_curve(sub)
    rocres = pROC::roc(sub$V4 ~sub$Sum, direction="<",levels = c("random","Experimental"),auc = T,plot=F, algorithm = 3,percent = T,main = paste(slop,feature),xlim=c(100,0))
    auc <- rbind(auc, data.frame("Feature" = feature, AUC = as.numeric(rocres$auc)))
  }
  AUC[[slop]] <-auc 
  }

AUC <- ldply(AUC,.id = "Window")

a <- (matlab.like(length(levels(AUC$Feature))))

gg1 <- ggplot(AUC,aes(Sample_name, factor(paste(Feature,Window),levels = paste(sort(rep(unique(Feature),length(c("0kb","01kb","1kb","10kb","100kb","1Mb")))),c("0kb","01kb","1kb","10kb","100kb","1Mb"))),fill = AUC/100)) + 
  geom_tile() + 
  scale_fill_gradientn(name = 'AUC',colours = colorRampPalette(c("blue","white", "red"))(10),limits=c(0,1)) + 
  theme_bw() + 
  theme(axis.text.y = element_text(hjust = 0,color = a[sort(as.numeric(AUC$Feature))]))+
  xlab("") + ylab("")

print(gg1)

dev.off()

write.table(x = AUC,file = paste("10-randomIS/",Sample_name,"_AUC.txt",sep = ""),quote = F,row.names = F,sep="\t")

rm(list = grep(ls(),pattern = "Sample_name|Corrected_FragLengthCount|args",invert = T,value = T))
