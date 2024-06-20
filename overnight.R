datasets <- c("GSE50660", "GSE90124", "GSE40279", "GSE67705", "GSE61256", "GSE73103", "GSE124366",
              "GSE114134", "GSE89253", "GSE53740", "GSE85568", "GSE106648", "GSE51057", "GSE124076")
overlaps <- overlappingSites(beta_path, datasets)
#321565 sites present in 11/14 datasets
write.csv(overlaps, "6_15_cpgfilter1.csv")

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
makeChunks(training_datasets, split_lists, training_subset_folder)
makeChunks(testing_datasets, split_lists, testing_subset_folder)

results_by_window <- filterCorrelations(num_groups, training=TRUE)

#print length of each list
for (i in 1:length(results_by_window)) {
  print(paste0("Start age: ", (i-1)*5 + 5, " End age: ", (i-1)*5+25, " Num CpG sites r>0.4: ", length(results_by_window[[i]])))
}

#print total unique sites

all_cpg_sites <- unlist(results_by_window)
unique_cpg_sites <- unique(all_cpg_sites)
print(length(unique_cpg_sites))

sorted_unique_cpg_sites <- sort(unique_cpg_sites)
write.csv(sorted_unique_cpg_sites, "6_14_cpgfilter2.csv")
#filter down to 224995

V <- generateMatrix("6_16_cpgfilter2v1.csv", num_splits=num_groups, training=TRUE, average=FALSE)
write.csv(V, "6_15_variance_matrix_v3.csv", row.names=TRUE)
B <- generateMatrix("6_16_cpgfilter2v1.csv", num_splits=num_groups, training=TRUE, average=TRUE)
write.csv(B, "6_15_training_matrix_v7.csv", row.names=TRUE)

second_filter <- filterBySDComparison("6_15_variance_matrix_v3.csv", "6_15_training_matrix_v7.csv") 
write.csv(second_filter, "6_16_cpgfilter3v1.csv")

D <- generateMatrix("6_16_cpgfilter3v1.csv", num_splits=num_groups, training=TRUE, average=TRUE)
write.csv(D, "6_16_filtered_matrix_train.csv", row.names=TRUE)

E <- generateMatrix("6_16_cpgfilter3v1.csv", num_splits=num_groups, training=FALSE, average=TRUE)
write.csv(E, "6_16_filtered_matrix_test.csv", row.names=TRUE)

#6/17/2024
library(dplyr)
library(data.table)
library(factoextra)
library(cluster)
library(corrplot)
library(tidyr)
library(tibble)
library(splines)
library(ggplot2)
library(dbscan)
library(kernlab)
library(apcluster)

# Load data
B <- readr::read_csv("6_16_filtered_matrix_train.csv")
G <- readr::read_csv("6_16_filtered_matrix_test.csv")

# Preprocessing
B <- B[, colSums(is.na(B)) != nrow(B)]
G <- G[, colSums(is.na(G)) != nrow(G)]
B <- na.omit(B)
G <- na.omit(G)
cpg_site_names <- B[[1]]
B <- B[, -1, with = FALSE]
B <- as.data.frame(B)
rownames(B) <- cpg_site_names
cpg_site_names <- G[[1]]
G <- G[, -1, with = FALSE]
G <- as.data.frame(G)
rownames(G) <- cpg_site_names
scaled_B <- t(scale(t(B)))
df <- as.data.frame(scaled_B)
scaled_G <- t(scale(t(G)))
df2 <- as.data.frame(scaled_G)
orderr <- rownames(df2)
df <- df[orderr,]
df <- na.omit(df)
orderr <- rownames(df)
df2 <- df2[orderr, ]

# Define clustering methods and distance measures
clustering_methods <- list(
  kmeans = function(data, k) kmeans(data, centers = k)$cluster,
  hclust = function(data, k) cutree(hclust(dist(data, method = "euclidean"), method = "ward.D2"), k = k)
)
distance_measures <- c("euclidean", "manhattan", "cosine")

# Function to determine the optimal number of clusters using Silhouette Method
optimal_clusters_silhouette <- function(data, clust_method, max_clusters = 10) {
  avg_sil <- sapply(2:max_clusters, function(k) {
    cl <- clust_method(data, k)
    ss <- silhouette(cl, dist(data))
    mean(ss[, 3])
  })
  plot(2:max_clusters, avg_sil, type = "b", xlab = "Number of Clusters (K)", ylab = "Average Silhouette Width")
  optimal_k <- which.max(avg_sil) + 1
  return(optimal_k)
}

# Loop through methods
for (dist_method in distance_measures) {
  d <- dist(df, method = dist_method)
  
  for (clust_name in names(clustering_methods)) {
    clust_method <- clustering_methods[[clust_name]]
    
    # Determine the optimal number of clusters
    num_clusters <- optimal_clusters_silhouette(df, clust_method)
    clusters <- clust_method(df, num_clusters)
    df_clustered <- df
    df_clustered$Cluster <- as.factor(clusters)
    
    df_clustered <- df_clustered[order(df_clustered$Cluster), ]
    df_order <- rownames(df_clustered)
    df2_ordered <- df2[df_order, ]
    
    for (i in 1:num_clusters) {
      train_df <- df_clustered[df_clustered$Cluster == i, ]
      ident <- rownames(train_df)
      test_df <- df2_ordered[ident, ]
      
      long_df_test <- test_df %>%
        rownames_to_column(var = "Sample") %>%
        pivot_longer(cols = -c(Sample), names_to = "Age", values_to = "Beta") %>%
        mutate(Age = as.numeric(sub("age_", "", Age)))
      
      long_df_train <- train_df %>%
        rownames_to_column(var = "Sample") %>%
        pivot_longer(cols = -c(Sample, Cluster), names_to = "Age", values_to = "Beta") %>%
        mutate(Age = as.numeric(sub("age_", "", Age)))
      
      fit <- lm(Beta ~ bs(Age, df = 5), data = long_df_train)
      fit2 <- lm(Beta ~ bs(Age, df = 5), data = long_df_test)
      long_df_test$Fitted <- predict(fit2)
      long_df_train$Fitted <- predict(fit)
      
      ggplot(long_df_test, aes(x = Age, y = Beta)) +
        geom_point(size=1) +
        geom_line(aes(y = Fitted), color = "blue") +
        labs(title = paste("Test Cluster", i, "-", dist_method, clust_name), x = "Age", y = "Scaled Beta Values") +
        theme_minimal()
      ggsave(
        paste0("test_cluster_", i, "_", dist_method, "_", clust_name, ".png"),
        plot = last_plot(),
        bg="white",
        scale = 1,
        dpi = 300
      )
      
      ggplot(long_df_train, aes(x = Age, y = Beta)) +
        geom_point(size=1) +
        geom_line(aes(y = Fitted), color = "blue") +
        labs(title = paste("Train Cluster", i, "-", dist_method, clust_name), x = "Age", y = "Scaled Beta Values") +
        theme_minimal()
      ggsave(
        paste0("train_cluster_", i, "_", dist_method, "_", clust_name, ".png"),
        plot = last_plot(),
        bg="white",
        scale = 1,
        dpi = 300
      )
    }
  }
}
