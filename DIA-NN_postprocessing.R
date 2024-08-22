## To compare protien abundances in baseline tissues
## between DDA (MaxQuant) and DIA (DIA-NN) experiments
## Computing from reports.tsv file output of DIA-NN
## 
######################
### USE THIS CODE ####
######################

library(data.table)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggvenn)
library(viridis)
library(heatmap3)

setwd('/Users/ananth/Documents/DIANN/')

#Read input file
DIA_datasets <- data.frame(Dataset = c(
  "PXD012254",
  "PXD001506",
  "PXD001764",
  "PXD025705",
  "PXD004684",
  "PXD032076",
  "PXD025431",
  "PXD002732",
  "PXD019594",
  "PXD022872",
  "PXD034908",
  "PXD018830",
  "PXD039665",
  "PXD033060",
  "PXD031419"
))


###################################################################
## Read and combine all datasets
###################################################################

Rootdata <- "SwissProt"

combined_reports <- list()

for(i in 1: nrow(DIA_datasets)){
  
  datasetID <- DIA_datasets[i,]
  print(datasetID)
  tmp  <- read.table(file=paste(datasetID, Rootdata, "report_AllSamples.tsv", sep="/"), quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
  tmp$File.Name <- gsub(".*PXD","PXD",tmp$File.Name,perl=TRUE)
  tmp$File.Name <- gsub("/submitted.*","",tmp$File.Name,perl=TRUE)

  tmp <- tmp[,c("File.Name","Run","Protein.Ids","Genes","Genes.Normalised","Stripped.Sequence")]
  
  # Don't consider peptides which map to more than one gene/protein identifier (gene groups)
  # Peptide mappings should be unique to a gene/protein identifier
  tmp <- tmp[!grepl(";",tmp$Genes),]
  tmp <- tmp[!grepl(";",tmp$Protein.Ids),]
  # Remove contaminants
  tmp <- tmp[!grepl("SWISS-PROT|\\(Bos|S.avidinii",tmp$Genes),]
  tmp <- tmp[!grepl("Streptavidin",tmp$Protein.Ids),]
  
  tmp$Genes.Normalised <- as.numeric(tmp$Genes.Normalised)
  
  tmp <- tmp[!grepl("Library|POOL|REPLICATE",tmp$Run, ignore.case=TRUE),]
  
  combined_reports[[i]] <- tmp
}

reports_submatrix <- bind_rows(combined_reports)
reports_submatrix <- reports_submatrix[reports_submatrix$Genes!="" & reports_submatrix$Protein.Ids!="" & reports_submatrix$File.Name!="" ,]

rm(tmp, combined_reports)

### Remove (-2016.* from datasets PXD031419 & PXD022872) from sample name
reports_submatrix$Run <- gsub("-2016.*","",reports_submatrix$Run, perl=TRUE)
reports_submatrix$Run <- gsub("-\\D+.*","",reports_submatrix$Run, perl=TRUE)

###################################################################
## Summarise number of PSMs mapped to each gene for every run(sample) in a dataset
###################################################################
PSMs_per_gene <- reports_submatrix  %>%
  group_by(File.Name,Run,Genes) %>%  #Protein.Ids can be removed to give accurate count
  summarise(Total_PSMs = n(),
            Distinct_StrippedSeq_PSMs = n_distinct(Stripped.Sequence))


Merged_PSMs_per_gene <- merge(x=reports_submatrix, y=PSMs_per_gene,
                              by.x=c("File.Name","Run","Genes"),
                              by.y=c("File.Name","Run","Genes"))

rm(reports_submatrix)

###################################################################
## Add tissue information of each dataset
###################################################################
DIA_tissues_dataset_annot <- read.table(file="DIA_tissues_datasets.tsv", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)

Merged_PSMs_per_gene <- merge(x=DIA_tissues_dataset_annot, y=Merged_PSMs_per_gene,
                                by.x=c("Dataset"), by.y=c("File.Name"),
                                all.x=FALSE, all.y=FALSE)

Merged_PSMs_per_gene$Tissues <- gsub("Esophageal epithelium","Esophagus", Merged_PSMs_per_gene$Tissues)

###################################################################
## Read sdrf and annotate sample metadata to each run(sample)
###################################################################

combined_sdrf <- list()

for(i in 1: nrow(DIA_datasets)){
  
  datasetID <- DIA_datasets[i,]
  
  sdrf_file <- list.files(path=datasetID, pattern = "\\.sdrf.txt")
  
  sdrf  <- read.table(file=paste(datasetID,sdrf_file, sep="/") , quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
  
  sdrf$Dataset <- datasetID
  sdrf <- sdrf[,grep("Dataset|Source.Name|Characteristics.disease|Characteristics.sampling.site", colnames(sdrf))]
  sdrf$Source.Name <- gsub("iBAQ.","",sdrf$Source.Name)
  colnames(sdrf) <- gsub("Characteristics.","",colnames(sdrf))
  colnames(sdrf) <- gsub("\\.","",colnames(sdrf))
  
  
  combined_sdrf[[i]] <- sdrf
}
all_sdrf <- bind_rows(combined_sdrf)
all_sdrf <- subset(all_sdrf, select=-c(diseasestaging))

# Extract only normal (baseline) samples
control_sdrf <- all_sdrf[all_sdrf$disease=="normal"|all_sdrf$disease=="insulin sensitive"|all_sdrf$samplingsite=="lung tumor-adjacent tissues",]
control_sdrf <- control_sdrf[grepl("^NA",row.names(control_sdrf))==FALSE,]

Merged_PSMs_per_gene <- merge(x=Merged_PSMs_per_gene, y=control_sdrf,
                              by.x=c("Run","Dataset"), by.y=c("SourceName","Dataset"),
                              all.x=FALSE, all.y=FALSE)

###################################################################
## Compute iBAQ values from DIA-NN LFQ : divide by number of theoretical tryptic peptides
###################################################################
library("Biostrings")

if(Rootdata == "SwissProt"){
  fasta_file <- readAAStringSet("Human_OneProteinPerGeneSet_May2023_UP000005640_9606_PLUS_Contaminants.fasta")} else {
    fasta_file <- readAAStringSet("Human_UniProt_ReferenceProteome_withIsoforms_March2023_PLUS_Contaminants.fasta")
  }

fasta_seq <- as.data.frame(fasta_file)
fasta_seq <- tibble::rownames_to_column(fasta_seq, "name")
colnames(fasta_seq) <- c("name","seq")
fasta_seq$Gene <- gsub(".* GN=|PE=.*","", fasta_seq$name, perl=TRUE)
fasta_seq$Gene <- gsub("^sp.*|tr.*","", fasta_seq$Gene, perl=TRUE)
fasta_seq$Gene <- gsub(" .*","", fasta_seq$Gene, perl=TRUE)
fasta_seq$PE <- gsub(".*PE=","", fasta_seq$name, perl=TRUE)
fasta_seq$PE <- gsub("^sp.*|tr.*","", fasta_seq$PE, perl=TRUE)
fasta_seq$PE <- gsub(" .*","", fasta_seq$PE, perl=TRUE)
fasta_seq$name <- gsub(" .*","", fasta_seq$name, perl=TRUE)
fasta_seq$number_of_tryptic_peptides <- "NA"

fragments <- list()

# Compute number of tryptic peptides in each sequence
# Look for Trypsin cleavage sites ----
for(i in 1:nrow(fasta_seq)){
  
  seq_id <-fasta_seq[i,1]
  seq <- fasta_seq[i,2]
  seq_length <- nchar(fasta_seq[i,2])
  
  print(seq_id)
  #get coordinates of Trypsin cleavage sites (Lysine or Arginine, but not followed by Proline)
  repl_seq <- seq
  repl_seq <- gsub("KP|RP","xx", repl_seq, perl=TRUE)
  seq_frag_coords <- data.frame(gregexpr(pattern ='K|R',repl_seq))
  fasta_seq[i, "number_of_tryptic_peptides"] <- nrow(seq_frag_coords)+1
}

fasta_seq$protein_id <-  gsub("sp\\||tr\\|", "", fasta_seq$name, perl=TRUE) 
fasta_seq$protein_id <-  gsub("\\|.*", "", fasta_seq$protein_id, perl=TRUE) 
fasta_seq$number_of_tryptic_peptides <- as.numeric(fasta_seq$number_of_tryptic_peptides)

write.table(fasta_seq, file = "SwissProt_TrypticPeptides.txt", sep = "\t", row.names = FALSE, quote = FALSE )
  

Filtered_PSMs_per_gene <- Merged_PSMs_per_gene[Merged_PSMs_per_gene$Distinct_StrippedSeq_PSMs > 1,]

Normalised_Filtered_PSMs_per_gene <- merge(x=Filtered_PSMs_per_gene, y=fasta_seq,
                                         by.x=c("Genes","Protein.Ids"), by.y=c("Gene","protein_id"),
                                         all.x=FALSE, all.y=FALSE)


# Normalise LFQ by tryptic peptides to get iBAQ
Normalised_Filtered_PSMs_per_gene$Genes.Normalised.Tryptic.iBAQ <- Normalised_Filtered_PSMs_per_gene$Genes.Normalised/Normalised_Filtered_PSMs_per_gene$number_of_tryptic_peptides

Normalised_Filtered_PSMs_per_gene_subset <- Normalised_Filtered_PSMs_per_gene[,c("Genes","Protein.Ids","Run","Dataset","Tissues","Genes.Normalised.Tryptic.iBAQ")]

Normalised_Filtered_PSMs_per_gene_subset_unique <- unique(Normalised_Filtered_PSMs_per_gene_subset)

Normalised_Filtered_PSMs_per_gene_subset_wide <- spread(Normalised_Filtered_PSMs_per_gene_subset_unique, Tissues, Genes.Normalised.Tryptic.iBAQ)

write.table(Normalised_Filtered_PSMs_per_gene_subset_wide, file = "report_AllDatasets_ControlSamplesOnly.tsv", sep = "\t", row.names = FALSE, quote = FALSE )


#Aggregate DIA abundances over Runs and Datasets
Normalised_Filtered_PSMs_per_dataset_aggregate <- Normalised_Filtered_PSMs_per_gene %>%
                                                  group_by(Tissues, Genes, Dataset) %>%
                                                  summarise(Median_Genes.Normalised.Tryptic.iBAQ = median(Genes.Normalised.Tryptic.iBAQ, na.rm=TRUE)) 


Normalised_Filtered_PSMs_per_gene_aggregate <- Normalised_Filtered_PSMs_per_gene %>%
                                                group_by(Tissues, Genes) %>%
                                                summarise(Median_Genes.Normalised.Tryptic.iBAQ = median(Genes.Normalised.Tryptic.iBAQ, na.rm=TRUE)) 

Normalised_Filtered_PSMs_per_gene_protein_aggregate <- Normalised_Filtered_PSMs_per_gene %>%
  group_by(Tissues, Genes, Protein.Ids) %>%
  summarise(Median_Genes.Normalised.Tryptic.iBAQ = median(Genes.Normalised.Tryptic.iBAQ, na.rm=TRUE)) 

Normalised_Filtered_PSMs_per_gene_aggregate_wide <- spread(Normalised_Filtered_PSMs_per_gene_aggregate, Tissues, Median_Genes.Normalised.Tryptic.iBAQ)
Normalised_Filtered_PSMs_per_gene_aggregate_wide <- Normalised_Filtered_PSMs_per_gene_aggregate_wide %>% mutate(Present_in_number_of_tissues = rowSums(!is.na(select(., -Genes))))

Normalised_Filtered_PSMs_per_gene_protein_aggregate_wide <- spread(Normalised_Filtered_PSMs_per_gene_protein_aggregate, Tissues, Median_Genes.Normalised.Tryptic.iBAQ)
Normalised_Filtered_PSMs_per_gene_protein_aggregate_wide <- Normalised_Filtered_PSMs_per_gene_protein_aggregate_wide %>% 
                                                            ungroup %>%
                                                            mutate(Present_in_number_of_tissues = rowSums(!is.na(select(., -Genes,-Protein.Ids))))

write.table(Normalised_Filtered_PSMs_per_gene_protein_aggregate_wide, file = "SupplementaryTable1_Canonical_protein_abundances_across_tissues_controlsamples.tsv", sep = "\t", row.names = FALSE, quote = FALSE )


################################
### Plots per datasets & Tissues
################################

# Fig.2A. Identified canonical proteins across organs 
ggplot(Normalised_Filtered_PSMs_per_gene_aggregate, aes(x=Tissues)) + 
  geom_bar(stat="count") + 
  xlab("")+
  ylab("Identified canonical proteins")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Distribution of canonical proteins across organs")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15))

# Fig.2B. Abundance distribution of abundances per organ
ggplot(Normalised_Filtered_PSMs_per_gene_aggregate, aes(x=Tissues, y=Median_Genes.Normalised.Tryptic.iBAQ)) + 
  geom_boxplot() + 
  xlab("")+
  ylab("protein abundance (iBAQ)")+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Distribution of canonical proteins across organs")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15))

# Fig.2C. Abundance distribution of abundances per dataset
ggplot(Normalised_Filtered_PSMs_per_dataset_aggregate, aes(x=Dataset)) + 
  geom_bar(stat="count") + 
  xlab("")+
  ylab("Identified canonical proteins")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Distribution of canonical proteins across organs")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15))

# Fig.2D. Abundance distribution of abundances per dataset
ggplot(Normalised_Filtered_PSMs_per_dataset_aggregate, aes(x=Dataset, y=Median_Genes.Normalised.Tryptic.iBAQ)) + 
  geom_boxplot() + 
  xlab("")+
  ylab("protein abundance (iBAQ)")+
  scale_y_log10()+
  theme_bw()+
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 5))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Distribution of canonical proteins across organs")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15))


#Fig.3. Heatmap  Pearson's correlation betweeen samples

binned_data <- read.table(file="report_AllDatasets_ControlSamplesOnly_binned_eachsample.txt", sep="\t", quote="\"", header=TRUE, stringsAsFactors = FALSE, comment.char = "#")

cordata_long <- gather(binned_data, Tissues, bin.iBAQ, Brain:Heart, factor_key=TRUE)
cordata_long$samples <- paste(cordata_long$Dataset,cordata_long$Run,cordata_long$Tissues, sep=".")
cordata_long <- cordata_long[,c("Genes","bin.iBAQ","samples")]
cordata_long <- cordata_long[!is.na(cordata_long$bin.iBAQ),]
cordata <- spread(cordata_long,samples,bin.iBAQ)

par(mar=c(5,5,2,2),xaxs = "i",yaxs = "i",cex.axis=0.4,cex.lab=10)

cor_results1 <- rcorr(as.matrix(cordata[,-c(1)]), type = c("pearson"))
correl_matrix <- cor_results1$r
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
correl_summary <- flattenCorrMatrix(cor_results1$r, cor_results1$P) 

Pairwise_only_Brain <- with(correl_summary, correl_summary[ grepl("Brain", row) & grepl("Brain", column),])
Pairwise_only_Brain$R2 <- Pairwise_only_Brain$cor*Pairwise_only_Brain$cor

Organs <- data.frame(ID= factor(gsub(".*\\.", "", colnames(cordata)[-1], perl=TRUE)))
annotation <- data.frame(Organs = gsub(".*\\.", "", colnames(correl_matrix), perl=TRUE))
annotation$Organs <- gsub("epithelium", "Breast epithelium", annotation$Organs)
annotation$Organs <- gsub("muscle", "Skeletal muscle", annotation$Organs)
rownames(annotation) <- colnames(correl_matrix)
annotation$Datasets <- gsub("\\..*", "", rownames(annotation), perl=TRUE)

pheatmap(correl_matrix, show_rownames = F, show_colnames = F,
         annotation_col = annotation,
         border_color = NA,
         fontsize = 6)


#############################
### Read DDA iBAQ values  ###
#############################
DDA_all <- read.table(file="/Users/ananth/Documents/MaxQuant_Bechmarking/Human/DDA_manuscript/JPR/SupplementaryFiles/SupplementaryTable_2.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
Protein_geneid_map <- DDA_all[,c("Majority.protein.IDs","Gene.Symbol","EnsemblID")]

DDA_all <- DDA_all[,c("Gene.Symbol","EnsemblID","Colon","Duodenum","Esophagus","Heart","Kidney","Liver","Lung","Pancreas","Thyroid","UterineEndometrium","Brain")]
colnames(DDA_all) <- gsub("UterineEndometrium","Uterus", colnames(DDA_all))

DDA_all <- DDA_all[!grepl(";", DDA_all$EnsemblID),]
DDA_all <- DDA_all[!grepl(";", DDA_all$Gene.Symbol),]

DDA_all_aggregate <- aggregate(DDA_all[ ,3:ncol(DDA_all)], list("Gene.ID" = DDA_all$EnsemblID, "Gene.Symbol" = DDA_all$Gene.Symbol), median, na.rm =TRUE)

DDA_all_long <- gather(DDA_all_aggregate, Samples, ppb.iBAQ, colnames(DDA_all_aggregate)[3]:colnames(DDA_all_aggregate)[ncol(DDA_all_aggregate)])

DDA_all_long$Tissues <- DDA_all_long$Samples

colnames(DDA_all_long)[4] <- "Median_ppb.iBAQ"

##########################################################
### Merge DDA and DIA abundances on tissue and gene IDs ###
##########################################################

Merged_DDA_DIA <- merge(x=DDA_all_long, y=Normalised_Filtered_PSMs_per_gene_aggregate,
                        by.x=c("Gene.Symbol","Tissues"), by.y=c("Genes","Tissues"),
                        all.x=FALSE, all.y=FALSE)

length(unique(Merged_DDA_DIA$Gene.Symbol))

## Fig.4. Correlations between DDA iBAQ and DIA iBAQ protein abundances
ggplot(Merged_DDA_DIA, aes(x=Median_Genes.Normalised.Tryptic.iBAQ, y=Median_ppb.iBAQ)) + 
  geom_point(size=0.5, alpha=0.2) + 
  scale_x_log10()+
  scale_y_log10()+
  xlab("Median iBAQ (DIA)")+
  ylab("Median ppb.iBAQ (DDA)")+
  theme_bw()+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.x.npc = 0.5,
           label.y.npc = 0.3)+
  geom_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x, size=0.4) +
  ggtitle("Abundance comparison (canonical geneID) \nppb.iBAQ(DDA-RefProt) v/s iBAQ.tryptic_peptides(DIA-SwissProt 20KSet)")+
  theme(axis.text=element_text(size=9),
        axis.title=element_text(size=12))+
  facet_wrap(~Tissues, scales = "free")

