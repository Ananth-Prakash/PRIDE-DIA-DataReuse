#Comparison with Jiang 2020 (TMT)
# Supplementary table 2
library(tidyr)
library(ggpubr)
library(ggplot2)

Jiang_2020 <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Normalised_protein_abundance_Jiang2020.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

colnames (Jiang_2020) <- gsub("\\.","", colnames(Jiang_2020), perl = TRUE)
#
Jiang_long <- gather(Jiang_2020, Organs, abundances, colnames(Jiang_2020)[3]:colnames(Jiang_2020)[ncol(Jiang_2020)], factor_key=TRUE)
Jiang_long$Organs <- gsub("\\d","", Jiang_long$Organs, perl = TRUE)
Jiang_long$Organs <- gsub("EsophagusGastroesophagealJunction|EsophagusMucosa|EsophagusMuscularis", "Esophagus", Jiang_long$Organs)
Jiang_long$Organs <- gsub("ArteryAorta|HeartAtrialAppendage|HeartLeftVentricle|ArteryCoronary", "Heart", Jiang_long$Organs)
Jiang_long$Organs <- gsub("ColonSigmoid|ColonTransverse", "Colon", Jiang_long$Organs)
Jiang_long$Organs <- gsub("BrainCerebellum|BrainCortex", "Brain", Jiang_long$Organs)
Jiang_long$Organs <- gsub("MinorSalivaryGland", "SalivaryGland", Jiang_long$Organs)
Jiang_long$Organs <- gsub("SmallIntestineTerminalIleum", "SmallIntestine", Jiang_long$Organs)

Jiang_long_aggregate <- aggregate(Jiang_long[,1:4], by=list(Jiang_long$geneidfull, Jiang_long$geneid, Jiang_long$Organs), median, na.rm =TRUE)

Jiang_long_aggregate <- Jiang_long_aggregate[,-c(4,5,6)]
colnames(Jiang_long_aggregate) <- c("GeneIDfull","GeneID","Organ","Abundance")
Jiang_long_aggregate <- Jiang_long_aggregate[!is.na(Jiang_long_aggregate$Abundance),]

Jiang_long_aggregate$Median_log2_Abundance <- log2(Jiang_long_aggregate$Abundance)

Jiang_long_aggregate <- Jiang_long_aggregate[,-c(4)]
Jiang_long_aggregate_wide <- spread(Jiang_long_aggregate, Organ, Median_log2_Abundance)
colnames(Jiang_long_aggregate_wide) <- gsub("$", "_Jiang_Etal_2020", colnames(Jiang_long_aggregate_wide), perl=TRUE)

#read median intensities from reanalysed datasets (ourstudy)
DIANN_abundances <- read.table(file = "/Users/ananth/Documents/DIANN/SupplementaryTable1_Canonical_protein_abundances_across_tissues_controlsamples.tsv", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

ID_maps <- read.table(file = "/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Bin_intensity_analysis/SupplementaryFile_3_UniProt_IDs.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
ID_maps <- ID_maps[,c(1:3)]

DIANN_abundances <- merge(x=ID_maps, y=DIANN_abundances,
                          by.x=c("GeneName"), by.y=c("Genes"), all.x=FALSE, all.y=FALSE)


Merged_data <- merge(x=DIANN_abundances, y=Jiang_long_aggregate_wide[-c(1)],
                     by.x=c("GeneID"), by.y=c("GeneID_Jiang_Etal_2020"), all.x=FALSE, all.y=FALSE)

Brain <- cbind(Merged_data[, c(1,2,grep("Brain", colnames(Merged_data)))], Organ="Brain")
colnames(Brain ) <- c("GeneID","GeneName","DIANN_abundances", "Jiang_etal_2020", "Organ")
Colon <- cbind(Merged_data[, c(1,2,grep("Colon", colnames(Merged_data)))], Organ="Colon")
colnames(Colon ) <- c("GeneID","GeneName","DIANN_abundances", "Jiang_etal_2020", "Organ")
Esophagus <- cbind(Merged_data[, c(1,2,grep("Esophagus", colnames(Merged_data)))], Organ="Esophagus")
colnames(Esophagus ) <- c("GeneID","GeneName","DIANN_abundances", "Jiang_etal_2020", "Organ")
Heart <- cbind(Merged_data[, c(1,2,grep("Heart", colnames(Merged_data)))], Organ="Heart")
colnames(Heart ) <- c("GeneID","GeneName","DIANN_abundances", "Jiang_etal_2020", "Organ")
Liver <- cbind(Merged_data[, c(1,2,grep("Liver", colnames(Merged_data)))], Organ="Liver")
colnames(Liver ) <- c("GeneID","GeneName","DIANN_abundances", "Jiang_etal_2020", "Organ")
Lung <- cbind(Merged_data[, c(1,2,grep("Lung", colnames(Merged_data)))], Organ="Lung")
colnames(Lung ) <- c("GeneID","GeneName","DIANN_abundances", "Jiang_etal_2020", "Organ")
Pancreas <- cbind(Merged_data[, c(1,2,grep("Pancreas", colnames(Merged_data)))], Organ="Pancreas")
colnames(Pancreas ) <- c("GeneID","GeneName","DIANN_abundances", "Jiang_etal_2020", "Organ")
SkeletalMuscle <- cbind(Merged_data[, c(1,2,grep("Skeletal", colnames(Merged_data)))], Organ="SkeletalMuscle")
colnames(SkeletalMuscle ) <- c("GeneID","GeneName","DIANN_abundances", "Jiang_etal_2020", "Organ")
Thyroid <- cbind(Merged_data[, c(1,2,grep("Thyroid", colnames(Merged_data)))], Organ="Thyroid")
colnames(Thyroid ) <- c("GeneID","GeneName","DIANN_abundances", "Jiang_etal_2020", "Organ")

plotdata <- rbind(Brain, Colon, Esophagus, Heart, Liver, Lung, Pancreas, SkeletalMuscle,Thyroid)


#Supplementary Fig. 4. Correrlation of protein abundances between DIA and Jiang et. al. 2020

ggplot(plotdata[complete.cases(plotdata),], aes(x=log2(DIANN_abundances), y=log2(Jiang_etal_2020))) + 
  geom_point(size=0.05, alpha=0.5) + 
  xlab("log2(iBAQ) this study")+
  ylab("log2(normalised abundance) Jiang etal 2020")+
  theme_bw()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.5,
           label.y.npc = 0.3)+
  geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Protein abundance comparison (This study (DIA) vs Jiang et al 2020 (TMT))")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15))+
  facet_wrap(~Organ, scales="free")


