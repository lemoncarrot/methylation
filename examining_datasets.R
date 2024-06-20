#' @title Formatting datasets prior to analysis
#' 
#' @author Kevin Chen

######################################################

#' @title Format datasets for matrix reduction pipeline
#' @description 
#' transform beta matrices and metadata into a format suitable for processing
#'
#' @param dataset_name GEO Accession number (e.g. GSE124366)
#' beta matrix (already processed) : first column CpG site identifiers
#'                                   following columns samples
#' metadata (already processed) : first column "sample_id"
#'                                second column "age"
#' 
#' @return formatted dataset is returned : last column CpG site identifiers
#'                                         first row ages
#'          
#' @export

# Load required libraries
library(data.table)

# Define input paths
metadata_path <- "/Users/kchen/OneDrive/Documents/methylation/metadata18/"
beta_path <- "/Users/kchen/OneDrive/Documents/methylation/betas18/"
# Define output path
combined_path <- "/Users/kchen/OneDrive/Documents/methylation/combined18/"

formatDataset <- function(dataset_name) {
  # Read beta matrix
  beta <- read.csv(paste0(beta_path, dataset_name, "_beta.csv"))
  # Set row names as CpG site identifiers
  rownames(beta) <- beta[, 1]
  beta <- beta[, -1]
  # Read metadata
  metadata <- read.csv(paste0(metadata_path, dataset_name, "_metadata_cleaned.csv"))
  # Cut metadata to sample ID and age
  metadata <- metadata[, c(1, 2)]
  t_metadata <- as.data.frame(t(metadata))
  colnames(t_metadata) <- t_metadata[1,]
  t_metadata <- t_metadata[-1, ]
  if (all(colnames(t_metadata)==colnames(beta))) {
    # Combine beta matrix and metadata
    combined <- rbind(t_metadata, beta)
  }
  else {
    print("Error")
  }
  #Set CpG identifiers to last column "cpg" (for data.table support)
  combined$cpg <- rownames(combined)
  return(combined)
}
# Run code
datasets <- c("GSE124076") #list of datasets to be formatted

for (dataset in datasets) {
  temp <- formatDataset(dataset)
  write.csv(temp, paste0(combined_path, dataset, "_df.csv"), row.names=FALSE)
  rm(temp)
  gc()
  print(paste0("DONE WITH: ", dataset))
}

#DO NOT RUN // TESTING CODE
#GSE124366 <- formatDataset("GSE124366")
#write.csv(GSE124366, paste0(combined_path, "GSE124366_df.csv"), row.names=FALSE)
#rm(GSE124366)
#gc()


#' @title Examine CpG sites present in datasets
#' @description 
#' examine which CpG sites are present in datasets, with a cutoff
#'
#' @param dataset_name GEO Accession number (e.g. GSE124366)
#' 
#' @return list of CpG sites passing a cutoff (present in at least 7/9 datasets)
#'          
#' @export

# Load required libraries
library(ggplot2)
library(data.table)

beta_path <- "/Users/kchen/OneDrive/Documents/methylation/betas18/"

overlappingSites <- function(beta_path, datasets) {
  all_sites <- character()
  #total_samples <- 0
  
  #LOop through list of datasets to be examined
  for (element in datasets) {
    filepath <- paste0(beta_path, element, "_beta.csv")
    #sample_names <- fread(filepath, nrows = 1, header = FALSE)[[1]]
    site_names <- fread(filepath, select = 1, skip = 1)[[1]]
    #total_samples <- total_samples + length(sample_names)
    #all sites
    all_sites <- c(all_sites, site_names)
    
    print(paste("Processing: ", element))
  }
  
  #create table of all sites and their frequency across datasets
  site_counts <- table(all_sites)
  cutoff_number <- 11
  overlap <- names(site_counts)[site_counts >= cutoff_number]
  
  output <- data.frame(cpg = overlap, stringsAsFactors=FALSE)
  return(output)
}

# Run code
datasets <- c("GSE50660", "GSE90124", "GSE40279", "GSE67705", "GSE61256", "GSE73103", "GSE124366",
              "GSE114134", "GSE89253", "GSE53740", "GSE85568", "GSE106648", "GSE51057", "GSE124076")
overlaps <- overlappingSites(beta_path, datasets)
#359280 sites in at least 6 out of 7 training datasets
#write.csv(overlaps, "5_18_cpgfilter.csv")

#check to ensure beta values are present
path <- "/Users/kchen/OneDrive/Documents/methylation/combined18/"
temp <- readr::read_csv(paste0(path, "GSE124366", "_df.csv"))
temp <- temp[-1, ]
all_values_in_range <- sapply(temp, function(x) {
  if (is.numeric(x)) {
    return(all(x >= 0 & x <= 1.1, na.rm = TRUE))
  } else {
    return(NA)  # or return(FALSE) if you want to treat non-numeric data as failing the check
  }
})
print(all_values_in_range)
#GSE85568 could be problematic (values slightly above 1.000)
#GSE89253 has M-values




#GSE72308_metadata <- read.csv("/Users/kchen/OneDrive/Documents/methylation/metadata18/GSE72308_metadata_cleaned.csv")
# DELETED FUNCTION convertAgesToYears -> dataset was removed
#GSE72308_metadata <- convertAgesToYears(GSE72308_metadata)
GSE114134_metadata <- read.csv("/Users/kchen/OneDrive/Documents/methylation/metadata18/GSE114134_metadata_cleaned.csv")
GSE50660_metadata <- read.csv("/Users/kchen/OneDrive/Documents/methylation/metadata18/GSE50660_metadata_cleaned.csv")
GSE89253_metadata <- read.csv("/Users/kchen/OneDrive/Documents/methylation/metadata18/GSE89253_metadata_cleaned.csv")
GSE53740_metadata <- read.csv("/Users/kchen/OneDrive/Documents/methylation/metadata18/GSE53740_metadata_cleaned.csv")
GSE90124_metadata <- read.csv("/Users/kchen/OneDrive/Documents/methylation/metadata18/GSE90124_metadata_cleaned.csv")
GSE40279_metadata <- read.csv("/Users/kchen/OneDrive/Documents/methylation/metadata18/GSE40279_metadata_cleaned.csv") #hannum
GSE67705_metadata <- read.csv("/Users/kchen/OneDrive/Documents/methylation/metadata18/GSE67705_metadata_cleaned.csv")
GSE85568_metadata <- read.csv("/Users/kchen/OneDrive/Documents/methylation/metadata18/GSE85568_metadata_cleaned.csv")
GSE106648_metadata <- read.csv("/Users/kchen/OneDrive/Documents/methylation/metadata18/GSE106648_metadata_cleaned.csv")
GSE61256_metadata <- read.csv("/Users/kchen/OneDrive/Documents/methylation/metadata18/GSE61256_metadata_cleaned.csv")
GSE51057_metadata <- read.csv("/Users/kchen/OneDrive/Documents/methylation/metadata18/GSE51057_metadata_cleaned.csv")
GSE73103_metadata <- read.csv("/Users/kchen/OneDrive/Documents/methylation/metadata18/GSE73103_metadata_cleaned.csv")
#GSE70977_metadata <- read.csv("/Users/kchen/OneDrive/Documents/methylation/metadata18/GSE70977_metadata_cleaned.csv")
#GSE53051_metadata <- read.csv("/Users/kchen/OneDrive/Documents/methylation/metadata18/GSE53051_metadata_cleaned.csv")
GSE124076_metadata <- read.csv("/Users/kchen/OneDrive/Documents/methylation/metadata18/GSE124076_metadata_cleaned.csv")
GSE124366_metadata <- read.csv("/Users/kchen/OneDrive/Documents/methylation/metadata18/GSE124366_metadata_cleaned.csv")
#GSE101764_metadata <- read.csv("/Users/kchen/OneDrive/Documents/methylation/metadata18/GSE101764_metadata_cleaned.csv")

#plot age distribution
ggplot(data=df_train, aes(x=age)) + 
  geom_histogram(binwidth=1, boundary=0, color="black", fill="blue") +
  scale_x_continuous(limits=c(0, 100), breaks=seq(0, 100, by=5)) +
  theme_minimal() +
  labs(x="Age", y="Samples", title="Training Age Distribution")

ggplot(data=df_test, aes(x=age)) + 
  geom_histogram(binwidth=1, boundary=0, color="black", fill="blue") +
  scale_x_continuous(limits=c(0, 100), breaks=seq(0, 100, by=5)) +
  theme_minimal() +
  labs(x="Age", y="Samples", title="Testing Age Distribution")


df_train <- rbind(GSE124366_metadata[2], GSE40279_metadata[2],
                  GSE61256_metadata[2], GSE67705_metadata[2], 
                  GSE90124_metadata[2], GSE50660_metadata[2], 
                  GSE73103_metadata[2])

df_test <- rbind(GSE114134_metadata[2], GSE89253_metadata[2], 
                 GSE53740_metadata[2], GSE85568_metadata[2],
                 GSE106648_metadata[2], GSE51057_metadata[2], 
                 GSE124076_metadata[2])

#making correctly formatted GSE89253 off GEO directly
library(data.table)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(Rhdf5lib)

#download RAW.tar file from GEO
#extract tar manually
baseDir <- "/Users/kchen/OneDrive/Documents/methylation/GSE89253_RAW/"

# List all .idat.gz files in the directory
idat_files <- list.files(baseDir, pattern = "\\.idat\\.gz$", full.names = TRUE)

# Reformatting idat.gz to idat
sapply(idat_files, function(f) {
  out_file <- sub("\\.gz$", "", f)
  if (!file.exists(out_file)) {
    message("Decompressing ", f)
    R.utils::gunzip(f, destname = out_file, remove = FALSE)
  }
})

idat_files <- list.files(baseDir, pattern = "\\.idat$", full.names = TRUE)
rgSet <- read.metharray.exp(base = baseDir)
mSet <- preprocessQuantile(rgSet)
betaValues <- getBeta(mSet)
beta_df <- as.data.frame(betaValues)
# Restructure sample IDs from GSE1293487_1234_2134124 to GSE1293487
names(beta_df) <- gsub("_.*", "", names(beta_df))
#write.csv(beta_df, "GSE89253_beta.csv", row.names=TRUE)

#plotting distribution of beta values
beta_path <- "/Users/kchen/OneDrive/Documents/methylation/betas18/"
temp <- readr::read_csv(paste0(beta_path, "GSE124366_beta.csv"))
temp <- temp[, -1]
df_long <- melt(temp)
colnames(df_long) <- c("Sample", "BetaValue")
ggplot(df_long, aes(x = BetaValue)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  labs(title = "Beta Distribution GSE106648", x = "Beta Value", y = "Count") +
  theme_minimal()


temp <- readr::read_csv("13059_2013_3156_MOESM3_ESM.csv")
positive_rows <- which(temp[, 23]=="positive")
positive <- temp[positive_rows, ]
negative_rows <- which(temp[, 23]=="negative")
negative <- temp[negative_rows, ]



plotSingleSiteModified <- function(cpg_name, split_lists, training=TRUE, positive=TRUE) {
  if (positive) {
    xx <- "positive"
  } else {
    xx <- "negative"
  }
  subset_folder <- testing_subset_folder
  if (training) {
    subset_folder <- training_subset_folder
  }
  split <- checkWhichSplit(cpg_name, split_lists)
  temp <- readr::read_csv(paste0(subset_folder, "partition_" , split, ".csv"))
  temp <- as.data.frame(temp)
  if (training) {
    rownames(temp) <- temp$cpg...466 #466 for training, 207 for testing
    
  }
  if (!training) {
    rownames(temp) <- temp$cpg...207 #466 for training, 207 for testing
    
  }
  temp <- temp[, grep("GSM", names(temp), value=TRUE)]
  ages <- as.numeric(temp[1,])
  cols <- which(ages<=100 & ages>=0)
  betas <- as.numeric(temp[cpg_name, cols])
  ages <- ages[cols]
  data_to_plot <- data.frame(Ages=ages, Betas=betas)
  
  # Filter data for the specified age ranges
  data_25_45 <- data_to_plot[data_to_plot$Ages >= 25 & data_to_plot$Ages <= 45, ]
  data_30_50 <- data_to_plot[data_to_plot$Ages >= 30 & data_to_plot$Ages <= 50, ]
  data_35_55 <- data_to_plot[data_to_plot$Ages >= 35 & data_to_plot$Ages <= 55, ]
  
  
  # Calculate correlation coefficients
  corr_25_45 <- round(cor(data_25_45$Ages, data_25_45$Betas, use="complete.obs", method="pearson"), 2)
  print(corr_25_45)
  corr_30_50 <- round(cor(data_30_50$Ages, data_30_50$Betas, use="complete.obs", method="pearson"), 2)
  print(corr_30_50)
  corr_35_55 <- round(cor(data_35_55$Ages, data_35_55$Betas, use="complete.obs", method="pearson"), 2)
  print(corr_35_55)
  
  ggplot(data_to_plot, aes(x=Ages, y=Betas)) +
    geom_point() + 
    geom_smooth(data = data_25_45, method = "lm", se = FALSE, color = "blue") +
    geom_smooth(data = data_30_50, method = "lm", se = FALSE, color = "red") +
    geom_smooth(data = data_35_55, method = "lm", se=FALSE, color = 'green') +
    labs(x="Age", y="Beta value", title=paste0(cpg_name, " (", xx, ")")) + 
    theme_minimal() +
    ylim(0, 1.2) +
    annotate("text", x = 10, y = 1.1, label = paste("25-45 Corr:", corr_25_45), color = "blue", hjust = 0) +
    annotate("text", x = 40, y = 1.1, label = paste("30-50 Corr:", corr_30_50), color = "red", hjust = 0) +
    annotate("text", x = 70, y = 1.1, label = paste("35-55 Corr:", corr_35_55), color='green', hjust=0)
  ggsave(
    paste0("6_20_", cpg_name, "_", xx, ".png"),
    plot = last_plot(),
    bg="white",
    scale = 1,
    dpi = 300
  )
}

for (i in 1:3) {
  x <- sample.int(nrow(positive), 1)
  site_name <- positive[x, 1]
  site <- as.character(site_name)
  plotSingleSiteModified(site, split_lists, training=TRUE, positive=TRUE)
}

for (i in 1:3) {
  x <- sample.int(nrow(negative), 1)
  site_name <- negative[x, 1]
  site <- as.character(site_name)
  plotSingleSiteModified(site, split_lists, training=TRUE, positive=FALSE)
}
