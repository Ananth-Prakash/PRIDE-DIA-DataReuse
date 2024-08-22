library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)

setwd("/Users/ananth/Documents/DIANN/")


Fraction_of_missing_values <- function(abundances){
  abundances[abundances == 0] <- NA
  sum_of_NAs <- sum(is.na(abundances))
  total_number_of_observations <- nrow(abundances)*ncol(abundances)
  fraction_of_missingness <- (sum_of_NAs/(total_number_of_observations))*100
  return(fraction_of_missingness)
}

Fraction_of_completeness <- function(abundances){
  num_of_samples <- ncol(abundances)
  df <- data.frame(completeness=apply(abundances, 1, function(x) length(which(!is.na(x)))) )
  df$fraction_of_completeness <- (df$completeness/num_of_samples)*100
  return(df)
}

complete_threshold=90

##############################################################
DDA_missing_reports <- list()
DDA_allsamples_ppb <- read.table(file="/Users/ananth/Documents/MaxQuant_Bechmarking/Human/Bin_intensity_analysis/Ppb_intensities_all_samples.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
colnames(DDA_allsamples_ppb) <- gsub("\\.","_", colnames(DDA_allsamples_ppb), perl=TRUE)
colnames(DDA_allsamples_ppb) <- gsub("Amygdala|AnteriorPitutaryGland|AnteriorTemporalLobe|CaudateNucleus|CerebellarHemisphericCortex|
                                      |Cerebellum|CorpusCallosum|DorsoLateralPreFrontalCortex|EntorhinalCortex|FrontalGyrus|FrontalPole|
                                      |InferiorParietalLobule|MiddleTemporalLobe|NeoCortex|OccipitalCortex|PinealGland|PitutaryHypophysis|
                                      |Precuneus|PreFrontalCortex|SubstantiaNigra|SuperiorTemporalGyrus|TemporalCortex|Thalamus|VisualCortex", "Brain", colnames(DDA_allsamples_ppb), perl=TRUE)
                                     
colnames(DDA_allsamples_ppb) <- gsub("Aorta|AorticValve|AtrialSeptum|InferiorVenaCava|
                                     |LeftAtrium|LeftVentricle|MitralValve|PulmonaryArtery|
                                     |PulmonaryValve|PulmonaryVein|RightAtrium|
                                     |RightVentricle|TricuspidValve|VentricularSeptum", "Heart", colnames(DDA_allsamples_ppb), perl=TRUE)

DDA_datasets <- data.frame(Dataset=colnames(DDA_allsamples_ppb[4:ncol(DDA_allsamples_ppb)]))
DDA_datasets$Tissue <- gsub(".*_","",DDA_datasets$Dataset,perl=TRUE)
DDA_datasets$Tissue <- gsub("\\..*","",DDA_datasets$Tissue,perl=TRUE)
DDA_datasets$Dataset <- gsub("_.*","",DDA_datasets$Dataset,perl=TRUE)
DDA_datasets <- unique(DDA_datasets)

for(i in 1:nrow(DDA_datasets)){
  Dataset <- DDA_datasets[i,c("Dataset")]
  Tissue <-  DDA_datasets[i,c("Tissue")]
  DDA_dataset_subset <- DDA_allsamples_ppb[,grep(Dataset,colnames(DDA_allsamples_ppb)), drop=FALSE]
  DDA_tissue_subset <- DDA_dataset_subset[,grep(Tissue,colnames(DDA_dataset_subset)),drop=FALSE]
  
  sample_abundances <- DDA_tissue_subset[rowSums(is.na(DDA_tissue_subset)) != ncol(DDA_tissue_subset), ,drop=FALSE]
  number_of_genes <- nrow(sample_abundances)
  number_of_samples <- ncol(sample_abundances)
  
  missingness <- Fraction_of_missing_values(sample_abundances)
  completeness <- Fraction_of_completeness(sample_abundances)
  completeness <- completeness[completeness$fraction_of_completeness >= complete_threshold,]
  norm_completeness_above_threshold <- (nrow(completeness)/(number_of_genes))*100
    
  report <- data.frame(Dataset=Dataset, Tissue=Tissue, 
                       Fraction_of_NA=missingness,
                       Normalised_completeness=norm_completeness_above_threshold)
  
  DDA_missing_reports[[i]] <- report
}

DDA_missingness_report <- bind_rows(DDA_missing_reports)
DDA_missingness_report$Type <- rep("DDA",nrow(DDA_missingness_report))
DDA_missingness_report[DDA_missingness_report$Fraction_of_NA == 0,]$Fraction_of_NA <- NA
DDA_missingness_report[DDA_missingness_report$Normalised_completeness == 100,]$Normalised_completeness <- NA

DDA_missingness <- DDA_missingness_report %>% group_by(Dataset) %>%
                   mutate(Fraction_of_NA_avg = mean(Fraction_of_NA, na.rm = TRUE),
                          Fraction_of_completeness_avg = mean(Normalised_completeness, na.rm = TRUE))

DDA_missingness <- unique(DDA_missingness[,c("Dataset","Fraction_of_NA_avg","Fraction_of_completeness_avg","Type")])

####### count number of samples and genes present in each dataset
DDA_per_dataset_report <- list()
DDA_datasets_only <- unique(DDA_datasets[,c("Dataset"), drop=FALSE])
for(i in 1:nrow(DDA_datasets_only)){
  Dataset <- DDA_datasets_only[i,c("Dataset")]
  DDA_dataset_subset <- DDA_allsamples_ppb[,grep(Dataset,colnames(DDA_allsamples_ppb)), drop=FALSE]
  sample_abundances <- DDA_dataset_subset[rowSums(is.na(DDA_dataset_subset)) != ncol(DDA_dataset_subset), ,drop=FALSE]
  number_of_genes <- nrow(sample_abundances)
  number_of_samples <- ncol(sample_abundances)
  report <- data.frame(Dataset=Dataset, Number_of_genes_present=number_of_genes,
                       Number_of_samples=number_of_samples)
  DDA_per_dataset_report[[i]] <- report
}
DDA_sample_genes_report <- bind_rows(DDA_per_dataset_report)
DDA_sample_genes_report$label <- paste0(DDA_sample_genes_report$Dataset," (",DDA_sample_genes_report$Number_of_genes_present,", ",DDA_sample_genes_report$Number_of_samples,")")

DDA_missingness <- merge(x=DDA_missingness, y=DDA_sample_genes_report,
                         by.x=c("Dataset"), by.y=c("Dataset"))


######################################################
####### Take DIA control sample abundances
combined_sdrf <- list()

DIA_dataset_list <- read.table(file="DIA_datasets_manuscript.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE)
for(i in 1: nrow(DIA_dataset_list)){
  
  datasetID <- DIA_dataset_list[i,]
  
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
control_sdrf$SourceName <- gsub("^","iBAQ.",control_sdrf$SourceName,perl=TRUE)


DIA_missing_reports <- list()
for(i in 1:nrow(DIA_dataset_list)){
  dataset <- DIA_dataset_list[i,]
  all_samples <-  read.table(file=paste("/Users/ananth/Documents/DIANN/",dataset,"/SwissProt/MappedToGeneID.txt",sep=""), quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE, check.names = FALSE)
  
  sdrf_dataset <- control_sdrf[control_sdrf$Dataset == dataset,]
  # Extract exact samples using word boundary
  control_samples <- paste0("\\b",sdrf_dataset[,c("SourceName")],"\\b")
  control_samples <- paste(control_samples, collapse="|")
  
  control_samples_abundance <- all_samples[ ,grep(control_samples, colnames(all_samples))]
  #Since only control samples are extracted, there could be NA's in all control samples for a gene.
  #therefore remove columns with all NAs
  control_samples_abundance <- control_samples_abundance[rowSums(is.na(control_samples_abundance)) != ncol(control_samples_abundance), ]
  number_of_genes <- nrow(control_samples_abundance)
  number_of_samples <- ncol(control_samples_abundance)
  
  missingness <- Fraction_of_missing_values(control_samples_abundance)
  
  completeness <- Fraction_of_completeness(control_samples_abundance)
  completeness <- completeness[completeness$fraction_of_completeness >= complete_threshold,]
  norm_completeness_above_threshold <- (nrow(completeness)/(number_of_genes))*100
  
  report <- data.frame(Dataset=dataset,Fraction_of_NA_avg=missingness,
                       Fraction_of_completeness_avg=norm_completeness_above_threshold,
                       Number_of_genes_present=number_of_genes,
                       Number_of_samples=number_of_samples)
  
  DIA_missing_reports[[i]] <- report
}

DIA_missingness <- bind_rows(DIA_missing_reports)
DIA_missingness$Type <- rep("DIA",nrow(DIA_missingness))
DIA_missingness$label <- paste0(DIA_missingness$Dataset," (",DIA_missingness$Number_of_genes_present,", ",DIA_missingness$Number_of_samples,")")
DIA_missingness <- DIA_missingness[,c("Dataset","Fraction_of_NA_avg","Fraction_of_completeness_avg","Type","Number_of_genes_present","Number_of_samples","label")]

Missingness_all <- rbind(DDA_missingness, DIA_missingness)

#Fractionation
Fractionation <- read.table(file="Fractionated_datasets.txt", quote = "\"", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", fill=TRUE, check.names = FALSE)
Missingness_all <- merge(x=Missingness_all, y=Fractionation,
                         by.x=c("Dataset","Type"), by.y=c("Dataset","Type"))

Missingness_all$label <- paste0(Missingness_all$Fractionation, " ", Missingness_all$label)
  
#Supplementary Fig.2A. Missing protein abundances
ggplot(Missingness_all, aes(x=reorder(label,Fraction_of_NA_avg), y=Fraction_of_NA_avg)) + 
  geom_col(aes(fill=Type))+ 
  xlab("Dataset")+
  ylab("\nPercentage of missing values")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Comparison of missing values (NA) across samples (control) in each dataset\n(Total NAs/Total observations) across all samples in all canon. prot. identified in a dataset")+
  facet_wrap(~Type, scale="free_x")

#Supplementary Fig.2B. Protein abundances present in 90% or more samples
ggplot(Missingness_all, aes(x=reorder(label,Fraction_of_NA_avg), y=Fraction_of_completeness_avg)) + 
  geom_col(aes(fill=Type))+ 
  xlab("Dataset")+
  ylab(paste0("Normalised percentage of canonical proteins with\n abundances observed in ",complete_threshold,"% or more samples"))+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Comparison of abundances across samples (control) in each dataset\n(Num. of canon. prot. where abund. observed in 90% or more samples/Total canon. prot.)")+
  facet_wrap(~Type, scale="free_x")