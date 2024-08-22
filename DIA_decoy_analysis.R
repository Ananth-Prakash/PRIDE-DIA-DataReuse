library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)

setwd("/Users/ananth/Documents/DIANN/")

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

Rootdata <- "SwissProt_with_ArabidopsisDecoys"

combined_reports <- list()
PSMs_dataset <- list()
Human_PSMs <- list()

for(i in 1: nrow(DIA_datasets)){
  
  datasetID <- DIA_datasets[i,]
  print(datasetID)
  tmp  <- read.table(file=paste(datasetID, Rootdata, "report.pg_matrix_with_ArabidopsisDecoys.tsv", sep="/"), quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
  colnames(tmp) <- gsub(".*mzML_files.","",colnames(tmp),perl=TRUE)

  #Remove contaminants
  tmp <- tmp[!grepl("SWISS-PROT",tmp$Genes),]
  tmp <- tmp[tmp$Genes!= "" & tmp$Genes!="(Bos" & tmp$Genes!="(S.avidinii)",]
  
  # Fetch protein groups which contains only Arabidopsis decoys (but not together with human proteins)
  decoys <- tmp[grepl("_ARATH",tmp$Protein.Names),]
  decoys <- decoys[!grepl("_HUMAN",decoys$Protein.Names),]
 
  decoys_long <- gather(decoys, samples, LFQ, colnames(decoys[6]):colnames(decoys[ncol(decoys)]), factor_key=TRUE)
  decoys_long$dataset <- rep(datasetID, nrow(decoys_long))
  decoys_long <- decoys_long[,c("Protein.Names","Genes","LFQ","dataset")]

  # Fetch human protein groups only
  human <- tmp[grepl("_HUMAN",tmp$Protein.Names),]
  human_long <- gather(human, samples, LFQ, colnames(human[6]):colnames(human[ncol(human)]), factor_key=TRUE)
  human_long$dataset <- rep(datasetID, nrow(human_long))
  human_long <- human_long[,c("Protein.Names","Genes","LFQ","dataset")]
  
  combined_reports[[i]] <- decoys_long
  Human_PSMs[[i]] <- human_long
  PSMs_dataset[[i]] <- data.frame(Dataset=datasetID, totalPSMs = nrow(tmp))
}
#total PSMs in each dataset
PSMs_dataset_counts <- bind_rows(PSMs_dataset)

#human proteins
human_submatrix <- bind_rows(Human_PSMs)

human_submatrix_aggregate <- aggregate(human_submatrix[ ,"LFQ"], list("ProteinNames"=human_submatrix$Protein.Names, "Genes"=human_submatrix$Genes, "Dataset"=human_submatrix$dataset), median, na.rm =TRUE)
human_submatrix_aggregate$Genes <- sapply(strsplit(human_submatrix_aggregate$Genes, ';'), function(x) toString(paste(gtools::mixedsort(x),collapse=";")))
human_submatrix_aggregate$ProteinNames <- sapply(strsplit(human_submatrix_aggregate$ProteinNames, ';'), function(x) toString(paste(gtools::mixedsort(x),collapse=";")))

human_submatrix_aggregate_wide <- spread(human_submatrix_aggregate, Dataset, x)

human_submatrix_aggregate_wide$present_in_number_of_samples <- "NA"
human_submatrix_aggregate_wide$present_in_number_of_samples <-  apply(human_submatrix_aggregate_wide[3:17], 1, function(x) length(which(!is.na(x))) )

human_pgroups_in_each_dataset <- data.frame(human_pcount=sapply(human_submatrix_aggregate_wide[3:17], function(x) length(which(!is.na(x))) ))
human_pgroups_in_each_dataset$Dataset <- row.names(human_pgroups_in_each_dataset)

# decoys
decoys_submatrix <- bind_rows(combined_reports)

decoys_submatrix_aggregate <- aggregate(decoys_submatrix[ ,"LFQ"], list("ProteinNames"=decoys_submatrix$Protein.Names, "Genes"=decoys_submatrix$Genes, "Dataset"=decoys_submatrix$dataset), median, na.rm =TRUE)
decoys_submatrix_aggregate$Genes <- sapply(strsplit(decoys_submatrix_aggregate$Genes, ';'), function(x) toString(paste(gtools::mixedsort(x),collapse=";")))
decoys_submatrix_aggregate$ProteinNames <- sapply(strsplit(decoys_submatrix_aggregate$ProteinNames, ';'), function(x) toString(paste(gtools::mixedsort(x),collapse=";")))

decoys_submatrix_aggregate_wide <- spread(decoys_submatrix_aggregate, Dataset, x)

decoys_submatrix_aggregate_wide$present_in_number_of_samples <- "NA"
decoys_submatrix_aggregate_wide$present_in_number_of_samples <-  apply(decoys_submatrix_aggregate_wide[3:17], 1, function(x) length(which(!is.na(x))) )

decoys_in_each_dataset <- data.frame(decoy_count=sapply(decoys_submatrix_aggregate_wide[3:17], function(x) length(which(!is.na(x))) ))
decoys_in_each_dataset$Dataset <- row.names(decoys_in_each_dataset)

decoys_in_each_dataset <- merge(x=decoys_in_each_dataset, y=PSMs_dataset_counts,
                                by.x=c("Dataset"), by.y=c("Dataset"))

All_counts <- merge(x=decoys_in_each_dataset,y=human_pgroups_in_each_dataset,
                    by.x=c("Dataset"), by.y=c("Dataset"))

decoys_in_each_dataset$normalised_percentage_decoy_count <- (decoys_in_each_dataset$decoy_count/decoys_in_each_dataset$totalPSMs)*100

# plot number of Arabidopsis protein group decoys present across datasets
ggplot(decoys_submatrix_aggregate_wide, aes(x=present_in_number_of_samples)) + 
  geom_bar(aes(fill="red") )+ 
  xlab("Present in number of datasets")+
  ylab("Arabidopsis decoys")+
  theme_bw()+
  theme(legend.position = "none")+
  stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-0.7) +
  ggtitle("Distribution of decoys (n=1367) collected from all datasets (n=15)")


# plot percentage of Arabidopsis protein group decoys in each dataset
# normalised by dividing by total PSMs in each dataset
ggplot(decoys_in_each_dataset, aes(x=Dataset,y=normalised_percentage_decoy_count)) + 
  geom_col(aes(fill="red") )+ 
  xlab("Dataset")+
  ylab("Percentage of Arabidopsis decoys in all PSMs")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(label=format(round(normalised_percentage_decoy_count,digits=1),nsmall=1)), vjust=-0.7) +
  ggtitle("Distribution of decoys (n=1367) in each dataset")

human_pgroups_count <- human_submatrix_aggregate_wide[,c("present_in_number_of_samples"), drop=FALSE]
human_pgroups_count$Group <- "Human proteins"
decoy_pgroups_count <- decoys_submatrix_aggregate_wide[,c("present_in_number_of_samples"), drop=FALSE]
decoy_pgroups_count$Group <- "Arabidopsis decoys"

pgroups_distribution <- rbind(human_pgroups_count, decoy_pgroups_count)

# Supplementary Fig.1. Arabidopsis decoys and Human protein group distribution present in datasets
ggplot(pgroups_distribution, aes(x=present_in_number_of_samples, y=..count.., fill=factor(Group))) + 
  geom_histogram(position = "dodge")+ 
  xlab("Present in number of datasets")+
  ylab("Protein group counts")+
  theme_bw()+
  theme(legend.position = c(0.8,0.7))+
  labs(fill = "Group")+
  stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-0.3, position = position_dodge(width = .9)) +
  ggtitle("Distribution of protein groups common in all datasets (n=15)")
