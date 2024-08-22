# To compare DIA protein abundances with ProteomicsDB.

library(stringr)
library(tidyr)
library(ggplot2)
library(ggpubr)

setwd("/Users/ananth/Documents/DIANN/")

ProtDB  <- read.table("/Users/ananth/Documents/MaxQuant_Bechmarking/Human/ProteomicsDB/Proteomicsdb_API_response.tsv" , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
ProtDB$normalised.intensity <- as.numeric(ProtDB$normalised.intensity)
ProtDB$normalised.intensity <- format(round(ProtDB$normalised.intensity, 3), nsmall = 2)
ProtDB$normalised.intensity <- as.numeric(ProtDB$normalised.intensity)
ProtDB <- ProtDB[,-c(1)]
ProtDB <- ProtDB[order(ProtDB$tissue),]

ProtDB <- unique(ProtDB)

ProtDB$Organ <- ProtDB$tissue
ProtDB$Organ <- gsub("cardia", "heart", ProtDB$Organ, perl=TRUE)
ProtDB$Organ <- gsub("cerebral cortex|prefrontal cortex", "brain", ProtDB$Organ, perl=TRUE)
ProtDB$Organ <- gsub("colon muscle", "colon", ProtDB$Organ, perl=TRUE)

ProtDB_filtered <- subset(ProtDB , Organ == "adrenal gland" | Organ == "adrenal gland"| Organ == "brain" | 
                            Organ == "colon" | Organ == "esophagus" | Organ == "gall bladder" | 
                            Organ == "heart" | Organ == "kidney" | Organ == "liver" | Organ == "lung"|
                            Organ == "lymph node" | Organ == "ovary" |Organ == "pancreas" | Organ == "placenta"|
                            Organ == "prostate gland" | Organ == "rectum" | Organ == "salivary gland" | 
                            Organ == "spleen" | Organ == "stomach" | Organ == "testis" | Organ == "thyroid gland"|
                            Organ == "tonsil" | Organ == "urinary bladder")

ProtDB_filtered <- ProtDB_filtered[,-c(2)]
ProtDB_filtered$Organ <- as.character(ProtDB_filtered$Organ)
ProtDB_filtered$UniProt <- as.character(ProtDB_filtered$UniProt)
ProtDB_filtered$Organ <- str_replace(ProtDB_filtered$Organ, "^\\w{1}", toupper)
ProtDB_filtered$Organ <- gsub("$","_ProteomicsDB", ProtDB_filtered$Organ, perl=TRUE)
ProtDB_filtered_aggregate <- aggregate(ProtDB_filtered[,2], by=list(ProtDB_filtered$UniProt, ProtDB_filtered$Organ), median, na.rm =TRUE)

colnames(ProtDB_filtered_aggregate) <- c("UniProt","Organ","normalised.intensity")
#
ProtDB_filtered_wide <- spread(ProtDB_filtered_aggregate, Organ, normalised.intensity)

#read median intensities from DIA reanalysed datasets (ourstudy)
DIANN_abundances <- read.table(file = "/Users/ananth/Documents/DIANN/SupplementaryTable1_Canonical_protein_abundances_across_tissues_controlsamples.tsv", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

DIANN_abundances <- DIANN_abundances[-c(ncol(DIANN_abundances))]
colnames(DIANN_abundances) <- gsub("$", "_this study", colnames(DIANN_abundances))
colnames(DIANN_abundances)[1] <- "GeneID"

#read uniprot ENSG mappings
ID_maps <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Bin_intensity_analysis/SupplementaryFile_3_UniProt_IDs.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
ID_maps <- ID_maps[,c(1:3)]

DIANN_abundances <- merge(x=ID_maps, y=DIANN_abundances,
                            by.x=c("GeneName"), by.y=c("GeneID"), all.x=FALSE, all.y=FALSE)

Merged_data <- merge(x=DIANN_abundances, y=ProtDB_filtered_wide,
                     by.x=c("UniProt_ID"), by.y=c("UniProt"), all.x=FALSE, all.y=FALSE)


Brain <- cbind(Merged_data[, c(1,2,3,grep("Brain", colnames(Merged_data)))], Organ="Brain")
colnames(Brain) <- c("UniProtID","GeneID","GeneName","DIANN", "ProteomicsDB", "Organ")
Colon <- cbind(Merged_data[, c(1,2,3,grep("Colon", colnames(Merged_data)))], Organ="Colon")
colnames(Colon) <- c("UniProtID","GeneID","GeneName","DIANN", "ProteomicsDB", "Organ")
Esophagus <- cbind(Merged_data[, c(1,2,3,grep("Esophagus", colnames(Merged_data)))], Organ="Esophagus")
colnames(Esophagus) <- c("UniProtID","GeneID","GeneName","DIANN", "ProteomicsDB", "Organ")
Heart <- cbind(Merged_data[, c(1,2,3,grep("Heart", colnames(Merged_data)))], Organ="Heart")
colnames(Heart) <- c("UniProtID","GeneID","GeneName","DIANN", "ProteomicsDB", "Organ")
Liver <- cbind(Merged_data[, c(1,2,3,grep("Liver", colnames(Merged_data)))], Organ="Liver")
colnames(Liver) <- c("UniProtID","GeneID","GeneName","DIANN", "ProteomicsDB", "Organ")
Lung <- cbind(Merged_data[, c(1,2,3,grep("Lung", colnames(Merged_data)))], Organ="Lung")
colnames(Lung) <- c("UniProtID","GeneID","GeneName","DIANN", "ProteomicsDB", "Organ")
Pancreas <- cbind(Merged_data[, c(1,2,3,grep("Pancreas", colnames(Merged_data)))], Organ="Pancreas")
colnames(Pancreas) <- c("UniProtID","GeneID","GeneName","DIANN", "ProteomicsDB", "Organ")
Thyroid <- cbind(Merged_data[, c(1,2,3,grep("Thyroid", colnames(Merged_data)))], Organ="Thyroid")
colnames(Thyroid) <- c("UniProtID","GeneID","GeneName","DIANN", "ProteomicsDB", "Organ")


plotdata <- rbind(Brain, Colon, Esophagus, Heart, Liver, Lung, Pancreas,Thyroid)

# Supplementary Fig.3. Comparison of protein abundances between DIA and ProteomicsDB
ggplot(plotdata[complete.cases(plotdata),], aes(x=log2(DIANN), y=log2(ProteomicsDB))) + 
  geom_point(size=0.05, alpha=0.5) + 
  xlab("log2(iBAQ) this study")+
  ylab("log2(normalised intensity) ProteomicsDB")+
  theme_bw()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.5,
           label.y.npc = 0.3)+
  geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Protein abundance comparison (This study (DIA) vs ProteomicsDB)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))+
  facet_wrap(~Organ, scales = "free")




