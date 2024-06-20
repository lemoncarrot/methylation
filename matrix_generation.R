#' @title Matrix Reduction
#' 
#' @author Kevin Chen

######################################################

#' @title Reduction into matrix format for input into hierarchical clustering
#' @description
#'
#' 

# Load required packages
library(data.table)
library(dplyr)


#START PRE INITIALIZING CHUNK
# all datasets
training_datasets <- c("GSE50660", "GSE90124", "GSE40279", "GSE67705", "GSE61256", "GSE73103", "GSE124366")

testing_datasets <- c("GSE114134", "GSE89253", "GSE53740", "GSE85568", "GSE106648", "GSE51057", "GSE124076")
# Define input path
combined_path <- "/Users/kchen/OneDrive/Documents/methylation/combined18/"
# Read desired CpG sites from overlappingSites()
initial_cpg_list <- read.csv("6_15_cpgfilter1.csv")[, 2] #could change in future based on cpg used filtering process
# Extract number of sites
NUMCPGSITE <- length(initial_cpg_list)
# Define chunk length
sites_per_list <- 30000
# Define group number
num_groups <- ceiling(NUMCPGSITE/sites_per_list)
#11 partitions as of 6/15
# Split CpG sites into chunks
split_lists <- vector("list", num_groups)
for (i in 1:num_groups) {
  start_index <- (i-1)*sites_per_list + 1
  end_index <- min(i*sites_per_list, NUMCPGSITE)
  split_lists[[i]] <- initial_cpg_list[start_index:end_index]
}

training_subset_folder <- "/Users/kchen/OneDrive/Documents/methylation/training_split/"
testing_subset_folder <- "/Users/kchen/OneDrive/Documents/methylation/testing_split/"
#END PRE INITIALIZING CHUNK

standardizeDataset <- function(df_name, desired_cpgs) {
  df <- readr::read_csv(paste0(combined_path, df_name, "_df.csv"))
  
  # Extract ages
  ages <- df[1, ]
  df <- df[-1, ]
  template <- data.table(cpg=desired_cpgs)
  # Reorder and filter by reference CpG sites
  df <- merge(template, df, by.x="cpg", by.y="cpg", all.x=TRUE, sort=FALSE)
  # Reattach ages
  df <- rbind(ages, df)
  
  return(df)
}

makeChunks <- function(datasets, split_site_list, output_folder) {
  count <- 1
  for (site_list in split_site_list) {
    temp <- data.frame()
    for (dataset in datasets) {
      df <- standardizeDataset(dataset, site_list)
      if (ncol(temp)==0) {
        temp <- df
      }
      else {
        temp <- cbind(temp, df)
      }
    }
    write.csv(temp, file=paste0(output_folder, "partition_", count, ".csv"), row.names=TRUE)
    count <- count + 1
  }
}

# Run code
training_subset_folder <- "/Users/kchen/OneDrive/Documents/methylation/training_split/"
testing_subset_folder <- "/Users/kchen/OneDrive/Documents/methylation/testing_split/"
makeChunks(training_datasets, split_lists, training_subset_folder)
makeChunks(testing_datasets, split_lists, testing_subset_folder)


checkWhichSplit <- function(cpg_name, split_lists) {
  for (split_index in seq_along(split_lists)) {
    if (cpg_name %in% split_lists[[split_index]]) {
      return(split_index)
    }
  }
  return(NA)
}

plotSingleSite <- function(cpg_name, split_lists, training=TRUE) {
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
    labs(x="Age", y="Beta value", title=paste(cpg_name)) + 
    theme_minimal() +
    ylim(0, 1.2) +
    annotate("text", x = 10, y = 1.1, label = paste("25-45 Corr:", corr_25_45), color = "blue", hjust = 0) +
    annotate("text", x = 40, y = 1.1, label = paste("30-50 Corr:", corr_30_50), color = "red", hjust = 0) +
    annotate("text", x = 70, y = 1.1, label = paste("35-55 Corr:", corr_35_55), color='green', hjust=0)
  ggsave(
    paste0("6_20_", cpg_name, ".png"),
    plot = last_plot(),
    bg="white",
    scale = 1,
    dpi = 300
  )
}


x <- sample.int(NUMCPGSITE, 1)
rand_site <- initial_cpg_list[x]

site <- "cg01584473"

plotSingleSite(site, split_lists, training=TRUE)


filterCorrelationsAndChange <- function(num_splits, training=TRUE) {
  subset_folder <- testing_subset_folder
  if (training) {
    subset_folder <- training_subset_folder
  }
  #initialize output list 
  results_by_window <- list()
  for (start_age in seq(5, 70, by=5)) {
    end_age <- start_age + 20
    if (end_age <= 100) {
      results_by_window[[paste(start_age, "to", end_age)]] <- vector("list")
    }
  }
  
  for (j in 1:num_splits) {
    temp <- readr::read_csv(paste0(subset_folder, "partition_" , j, ".csv")) 
    temp <- as.data.frame(temp)
    if (!training) {
      rownames(temp) <- temp$cpg...207 #466 for training, 207 for testing
    }
    if (training) {
      rownames(temp) <- temp$cpg...466 #466 for training, 207 for testing
      
    }
    temp <- temp[, grep("GSM", names(temp), value=TRUE)]
    ages <- as.numeric(temp[1, ])
    betas <- as.matrix(temp[-1, ])


    for (start_age in seq(5, 70, by=5)) {
      end_age <- start_age + 20 
      print(paste0("split: ", j, " startage: ", start_age, " endage: ", end_age))
      age_window <- which(ages>=start_age & ages<=end_age)
      
      #check if window is within bounds (sample present)
      if (end_age > 100 || length(age_window) == 0) {
        next
      }
      
      betas_sub <- betas[, age_window, drop = FALSE]
      ages_sub <- ages[age_window]

      if (ncol(betas_sub) < 2 || length(ages_sub) < 2) {
        next  # Skip if there aren't enough data points to compute correlation
      }
      
      cor_results <- apply(betas_sub, 1, function(x) {
        if (sum(!is.na(x)) < 2) {
          return(NA)
        } else {
          return(cor(x, ages_sub, use="complete.obs", method="pearson"))
        }
      })
      
      beta_changes <- apply(betas_sub, 1, function(x) {
        max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
      })
      
      sig_indices <- which(abs(cor_results) > 0.4 & beta_changes>=0.2)
      if (length(sig_indices) > 0) {
        print("Filling results")
        selected_cpg_sites <- rownames(betas)[sig_indices]
        key <- paste(start_age, "to", end_age)
        results_by_window[[key]] <- unique(c(results_by_window[[key]], selected_cpg_sites))
      }
    }
  }
  return(results_by_window)
}

results_by_window <- filterCorrelationsAndChange(num_groups, training=TRUE)

#print length of each list
for (i in 1:length(results_by_window)) {
  print(paste0("Start age: ", (i-1)*5+5, " End age: ", (i-1)*5+25, " Num CpG sites r>0.4 & Abs Change>=0.15: ", length(results_by_window[[i]])))
}

#print total unique sites

all_cpg_sites <- unlist(results_by_window)
unique_cpg_sites <- unique(all_cpg_sites)
print(length(unique_cpg_sites))

sorted_unique_cpg_sites <- sort(unique_cpg_sites)
write.csv(sorted_unique_cpg_sites, "6_19_cpgfilter_changes.csv")

sites_to_investigate <- results_by_window["25 to 45"]
sites_to_investigate <- unlist(sites_to_investigate)
#length(sites_to_investigate)

sites_to_investigate <- results_by_window[c(1, 2)]
both <- intersect(sites_to_investigate[1], sites_to_investigate[2])
sites_to_investigate <- unlist(sites_to_investigate)
freq <- count(sites_to_investigate)
#random selection of sites

sites_to_investigate <- unlist(results_by_window["25 to 45"])
x <- sample.int(length(sites_to_investigate), 1)
rand_site <- sites_to_investigate[x]

plotSingleSite(rand_site, split_lists)

site1 <- results_by_window[["25 to 45"]]
site2 <- results_by_window[["30 to 50"]]

overlap <- intersect(site1, site2)



#old arbitrary variance filter
varianceFilter <- function(num_splits, sites, threshold=0.01, training=TRUE) {
  subset_folder <- testing_subset_folder
  if (training) {
    subset_folder <- training_subset_folder
  }
  
  all_sites <- list()
  
  for (j in 1:num_splits) {
    temp <- readr::read_csv(paste0(subset_folder, "partition_" , j, ".csv")) 
    temp <- as.data.frame(temp)
    des_rows <- which(temp$cpg...207 %in% sites)
    if (training) {
      des_rows <- which(temp$cpg...466 %in% sites)
    }
    rownames(temp) <- temp$cpg...466 #466 for training, 207 for testing
    temp <- temp[des_rows, ]
    temp <- temp[, grep("GSM", names(temp), value=TRUE)]

    variances <- apply(temp, 1, function(x) var(x, na.rm=TRUE))
    selected_variances <- variances[variances>threshold]
    selected_cpg_sites <- names(selected_variances)
    all_sites <- c(all_sites, selected_cpg_sites)
    
  }
  return(all_sites)
}

sites <- read.csv("6_19_cpgfilter_changes.csv")
sites <- sites$x
sites_filtered_2 <- varianceFilter(num_groups, sites, training=TRUE)
write.csv(sites_filtered_2, "6_19_cpgfilter2.csv")
    
#5_18_filtered_sites.csv is filtered with correlations
#5_28 is filtered with correlations AND arbitrary 0.015 variance
#comparison matrix is training_matrix_v3.csv, needs to be manipulated a bit
generateMatrix <- function(file_with_pre_filter, num_splits=num_groups, training=TRUE, average) {
  sites <- read.csv(file_with_pre_filter)
  sites <- sites[, -1]
  sites <- as.character(sites)
  #length: 288962
  NUMCPGSITE2 <- length(sites)
  #initialize output matrixs
  B <- matrix(NA, nrow=NUMCPGSITE2, ncol=101)
  age_range <- 0:100
  rownames(B) <- sites
  column_names <- character(101)
  for (i in age_range) {
    column_names[i+1] <- paste0("age_", i)
  }
  colnames(B) <- column_names
  subset_folder <- testing_subset_folder
  if (training) {
    subset_folder <- training_subset_folder
  }
  
  for (i in 1:num_splits) {
    temp <- readr::read_csv(paste0(subset_folder, "partition_", i, ".csv"))
    # Remove extra first column (automatic indices)
    temp <- temp[, -1]
    #filter for sites from correlation function
    #466 for training, 207 for testing
    if (training) {
      des_rows <- which(temp$cpg...466 %in% sites)
    }
    if (!training) {
      des_rows <- which(temp$cpg...207 %in% sites)
      
    }
    temp <- temp[c(1, des_rows), ] #keep ages in the first row
    # Extract CpG site list
    identifiers <- temp[-1, ncol(temp)]
    for (j in 0:100) {
      # DEBUG // print(paste0("Split: ", i, " Age: ", j))
      des_cols <- which(temp[1, ]==j)
      # DEBUG // print(paste0("Columns at Age ", j, " : ", length(des_cols)))
      # Check if there are samples at age j
      if (length(des_cols)>0) {
        # Filter for age j
        age_data <- temp[, des_cols, drop = FALSE]
        # Drop age row
        age_data <- age_data[-1, ]
        
        #indice of identifiers corresponds to indice of age_data
        #age_data has samples as column names, no other identifiers
        # DEBUG // print(paste0("Dim of age_data: ", dim(age_data)))
        #this step takes the longest
        setDT(age_data)
        # Vectorized averaging
        if (average) {
          aggregated <- age_data[, .(final = rowMeans(.SD, na.rm = TRUE)), by = .I]
        }
        if (!average) {
          aggregated <- age_data[, .(final = apply(.SD, 1, sd, na.rm = TRUE)), by = .I]
        }        
        # DEBUG // print(paste0("Averaging done"))
        aggregated <- as.data.frame(aggregated)
        if (training) {
          rownames(aggregated) <- as.character(identifiers$cpg...2441) #2441 for training, 1512 for testing
        }
        else {
          rownames(aggregated) <- as.character(identifiers$cpg...1512) #2441 for training, 1512 for testing
        }
        # Reassignment to matrix B
        B[rownames(aggregated), paste0("age_", j)] <- aggregated$final
        # DEBUG //print(paste0("Split: ", i, ", Age: ", j, " updated"))
      }
    }
  }
  return(B)
}

V <- generateMatrix("6_19_cpgfilter2.csv", num_splits=num_groups, training=TRUE, average=FALSE)
write.csv(V, "6_20_variance_matrix_v1.csv", row.names=TRUE)
B <- generateMatrix("6_19_cpgfilter2.csv", num_splits=num_groups, training=TRUE, average=TRUE)
write.csv(B, "6_20_training_matrix_v1.csv", row.names=TRUE)

transformAVG <- function(x) {
  if (!is.na(x) && x >= 0 && x <= 1) {
    return(2 * sqrt(x * (1 - x)))
  } else {
    return(NA)  # Handle non-numeric or out-of-bound values appropriately
  }
}

filterBySDComparison <- function(path_to_variance_matrix, path_to_avg_matrix, transformAVG=transformAVG) {
  final_site_list <- list()
  #assume already formatted correctly here
  var <- readr::read_csv(path_to_variance_matrix)
  #var <- var[, colSums(is.na(var)) != nrow(var)]
  cpg_site_names <- var[[1]]
  var <- var[, -1, with = FALSE]
  var <- as.data.frame(var)
  rownames(var) <- cpg_site_names
  avg <- readr::read_csv(path_to_avg_matrix)
  #avg <- avg[, colSums(is.na(avg)) != nrow(avg)]
  avg <- avg[, -1, with = FALSE]
  avg <- as.data.frame(avg)
  rownames(avg) <- cpg_site_names
  #t <- apply(avg, c(1, 2), transformAVG)
  #rownames(t) <- rownames(avg)
  #colnames(t) <- colnames(avg)
  for (i in 1:nrow(avg)) {
    #DEBUGGING
    print(paste("ROW", i))
    count <- 0
    for (j in 1:ncol(avg)) {
      avg_value <- avg[i, j]
      adj <- 0.5 * sqrt(avg_value * (1 - avg_value)) # TWEAK VALUE FUNCTION
      variance <- var[i, j]
      if (!is.na(avg_value) && !is.na(adj) && !is.na(variance) && (variance > adj)) {
        count <- count + 1
      }
    }
    if (count<2) {
      final_site_list <- append(final_site_list, rownames(avg)[i])
      print(paste("Added Site", rownames(avg)[i]))
    } else {
      print("Too High Variance")
    }
  }
  return(final_site_list)
}

second_filter <- filterBySDComparison("6_20_variance_matrix_v1.csv", "6_20_training_matrix_v1.csv") 
write.csv(second_filter, "6_20_cpgfilter3.csv")

Z <- generateMatrix("6_20_cpgfilter3.csv", num_splits=num_groups, training=TRUE, average=TRUE)
write.csv(Z, "6_20_training_matrix_v2.csv", row.names=TRUE)

P <- generateMatrix("6_20_cpgfilter3.csv", num_splits=num_groups, training=FALSE, average=TRUE)
write.csv(P, "6_20_testing_matrix_v2.csv", row.names=TRUE)




generateVarianceMatrix <- function(num_splits, sites, training=TRUE) {
  #start with list from correlation filter
  #need to generate and fill a matrix of variances at that age
  # (read all samples of that age and site, calculate SD, fill in entry in matrix)
  # (generate "ideal" variance matrix using 2sqrt(avg_beta * (1-avg_beta)))
  # (only add to final site list if site has SD more than ideal at 1 or no time points)
  
  
}


results_by_window <- filterCorrelations(num_groups)
all_cpg_sites <- unlist(results_by_window)
unique_cpg_sites <- unique(all_cpg_sites)
print(length(unique_cpg_sites))
sorted_unique_cpg_sites <- sort(unique_cpg_sites)
write.csv(sorted_unique_cpg_sites, "6_14_filtered_sites.csv")

V <- generateMatrix("6_14_filtered_sites.csv", num_splits=18, training=TRUE, average=FALSE)
write.csv(V, "variance_matrix_v2.csv", row.names=TRUE)
B <- generateMatrix("6_14_filtered_sites.csv", num_splits=18, training=TRUE, average=TRUE)
write.csv(B, "training_matrix_v5.csv", row.names=TRUE)

second_filter <- filterBySDComparison("variance_matrix_v2.csv", "training_matrix_v5.csv") 
write.csv(second_filter, "6_14_second_filter.csv")

Z <- generateMatrix("6_14_second_filter.csv", num_splits=18, training=TRUE, average=TRUE)
write.csv(Z, "6_14_filtered_matrix.csv", row.names=TRUE)

C <- generateMatrix("6_14_second_filter.csv", num_splits=18, training=FALSE, average=TRUE)
write.csv(C, "6_14_filtered_matrix_test.csv", row.names=TRUE)







