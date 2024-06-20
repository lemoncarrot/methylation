#take in matrix B with columns as ages and rows as cpg sites

library(dplyr)
library(data.table)
library(factoextra)
library(pheatmap)
library(cluster)


# Reformatting
B <- readr::read_csv("reduced_training_matrix_v2.csv") 
# Remove ages with no samples
B <- B[, colSums(is.na(B)) != nrow(B)]
# Remove rows with NA values
B <- na.omit(B)
# Reassign CpG site names to rownames and convert back to "data frame"
cpg_site_names <- B[[1]]
B <- B[, -1, with = FALSE]
B <- as.data.frame(B)
rownames(B) <- cpg_site_names

# Seed for subset1
set.seed(928)
# Random sampling for size reduction
sampled_indices <- sample(nrow(B), size = 20000)
subset1_B <- B[sampled_indices, ]
# Scaling
subset1_B_scaled <- t(scale(t(subset1_B)))

dissimilarity <- as.dist(1 - cor(t(subset1_B_scaled)))
hc <- hclust(dissimilarity, method = "ward.D2")
plot(hc, main = "Hierarchical Clustering with Ward.D2 Method")

sil_widths <- sapply(2:10, function(k) {
  clusters <- cutree(hc, k)
  silhouette_avg <- mean(silhouette(clusters, dissimilarity)[, "sil_width"])
  return(silhouette_avg)
})
plot(2:10, sil_widths, type = 'b', pch = 19, xlab = 'Number of clusters', ylab = 'Average silhouette width',
     main = 'Silhouette Analysis for Optimal Clusters')


# Choose the number of clusters with the highest average silhouette width
optimal_clusters <- which.max(sil_widths) + 1  # +1 because our k starts at 2
clusters <- cutree(hc, k = optimal_clusters)
subset1_B_scaled <- as.data.frame(subset1_B_scaled)
subset1_B_scaled$Cluster <- as.factor(clusters)

temp <- subset1_B_scaled[which(subset1_B_scaled$Cluster==2),]
# Remove last column (Cluster identifiers)
temp <- temp[, -ncol(temp)]

# Plot heatmap of cluster
pheatmap(temp,
         scale = "none",
         show_rownames = FALSE,
         show_colnames = TRUE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize_col = 6,
         angle_col = 45)

# Seed for representative CpG selection
set.seed(9135) #randomly select a seed each time
#randomly selected sub1clus1 - cg04671739, cg09125402, cg14341418
#randomly selected sub1clus2 - cg13355704, cg01164664, cg16587794
rep_sample_index <- sample(1:nrow(temp), 1, replace=F)
rep_sample_name <- rownames(temp)[rep_sample_index]
print(rep_sample_name)
cp_site <- temp[rep_sample_index, ] 
cp_site <- as.numeric(cp_site[1, ])
ages <- as.numeric(gsub("age_", "", colnames(temp)))

loess_fit <- loess(cp_site ~ ages, span = 0.6)  # You can adjust the 'span' for desired smoothness
plot(ages, cp_site, xlab = "Age", ylab = "Scaled Beta", main = paste0("Loess Fit for ", rep_sample_name),
     pch = 20, col = "blue", ylim = c(min(cp_site), max(cp_site)))
lines(ages, predict(loess_fit, data.frame(ages = ages)), col = "red", lwd = 2)

# Seed for subset1
set.seed(5145)
# Random sampling for size reduction
sampled_indices2 <- sample(nrow(B), size = 20000)
subset2_B <- B[sampled_indices2, ]
# Scaling
subset2_B_scaled <- t(scale(t(subset2_B)))

dissimilarity <- as.dist(1 - cor(t(subset2_B_scaled)))
hc <- hclust(dissimilarity, method = "ward.D2")
plot(hc, main = "Hierarchical Clustering with Ward.D2 Method")

#sil_widths <- sapply(2:10, function(k) {
#  clusters <- cutree(hc, k)
#  silhouette_avg <- mean(silhouette(clusters, dissimilarity)[, "sil_width"])
#  return(silhouette_avg)
#})
#plot(2:10, sil_widths, type = 'b', pch = 19, xlab = 'Number of clusters', ylab = 'Average silhouette width',
#     main = 'Silhouette Analysis for Optimal Clusters')


# Choose the number of clusters with the highest average silhouette width
optimal_clusters <- which.max(sil_widths) + 1  # +1 because our k starts at 2
clusters <- cutree(hc, k = optimal_clusters)
subset2_B_scaled <- as.data.frame(subset2_B_scaled)
subset2_B_scaled$Cluster <- as.factor(clusters)

temp2 <- subset2_B_scaled[which(subset2_B_scaled$Cluster==2),]
# Remove last column (Cluster identifiers)
temp2 <- temp2[, -ncol(temp2)]
pheatmap(temp2,
         scale = "none",
         show_rownames = FALSE,
         show_colnames = TRUE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize_col = 6,
         angle_col = 45)

# Seed for representative CpG selection
set.seed(5152) #randomly select a seed each time
#randomly selected sub2clus1 - cg10162251, cg16963622, cg06290884
#randomly selected sub2clus2 - cg04217778, cg07581978, cg04010091
rep_sample_index <- sample(1:nrow(temp2), 1, replace=F)
rep_sample_name <- rownames(temp2)[rep_sample_index]
print(rep_sample_name)
cp_site <- temp2[rep_sample_index, ] 
cp_site <- as.numeric(cp_site[1, ])
ages <- as.numeric(gsub("age_", "", colnames(temp2)))

loess_fit <- loess(cp_site ~ ages, span = 0.6)  # You can adjust the 'span' for desired smoothness
plot(ages, cp_site, xlab = "Age", ylab = "Scaled Beta", main = paste0("Loess Fit for ", rep_sample_name),
     pch = 20, col = "blue", ylim = c(min(cp_site), max(cp_site)))
lines(ages, predict(loess_fit, data.frame(ages = ages)), col = "red", lwd = 2)

# Validation in testing set
testing_B <- readr::read_csv("reduced_testing_matrix_v1.csv") 
# Remove ages with no samples
testing_B <- testing_B[, colSums(is.na(testing_B)) != nrow(testing_B)]
# Remove rows with NA values
testing_B <- na.omit(testing_B)
# Reassign CpG site names to rownames and convert back to "data frame"
cpg_site_names_testing <- testing_B[[1]]
testing_B <- testing_B[, -1, with = FALSE]
testing_B <- t(scale(t(testing_B)))
testing_B <- as.data.frame(testing_B)
rownames(testing_B) <- cpg_site_names_testing

site_to_validate <- "cg11084334"
#print(site_to_validate)
cp_site <- testing_B[site_to_validate, ]
cp_site <- as.numeric(cp_site[1, ])
#cp_site <- t(scale(t(cp_site)))
ages <- as.numeric(gsub("age_", "", colnames(testing_B)))
loess_fit <- loess(cp_site ~ ages, span = 0.6)  # You can adjust the 'span' for desired smoothness
plot(ages, cp_site, xlab = "Age", ylab = "Scaled Beta", main = paste0("Loess Fit for ", site_to_validate),
     pch = 20, col = "blue", ylim = c(min(cp_site), max(cp_site)))
lines(ages, predict(loess_fit, data.frame(ages = ages)), col = "red", lwd = 2)

set.seed(125)
# Random sampling for size reduction
sampled_indices <- sample(nrow(testing_B), size = 20000)
subset <- testing_B[sampled_indices, ]
pheatmap(subset,
         scale = "none",
         show_rownames = FALSE,
         show_colnames = TRUE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize_col = 6,
         angle_col = 45)
