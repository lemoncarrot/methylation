#take in matrix B with columns as ages and rows as cpg sites

library(dplyr)
library(data.table)
library(factoextra)
library(pheatmap)
library(cluster)
library(corrplot)
library(tidyr)
library(tibble)
library(sparcl)

# Reformatting
#B <- readr::read_csv("training_matrix_v3.csv") 
#B <- readr::read_csv("testing_matrix_v3.csv")
B <- readr::read_csv("6_20_training_matrix_v2.csv")
G <- readr::read_csv("6_20_testing_matrix_v2.csv")
# Remove ages with no samples
# reduced to 89 cols
B <- B[, colSums(is.na(B)) != nrow(B)]
G <- G[, colSums(is.na(G)) != nrow(G)]

# Remove rows with NA values
#reduced to 227963 rows, could add imputation
#97958 in testing
B <- na.omit(B)
G <- na.omit(G)

# Reassign CpG site names to rownames and convert back to "data frame"
cpg_site_names <- B[[1]]
B <- B[, -1, with = FALSE]
B <- as.data.frame(B)
rownames(B) <- cpg_site_names
#df <- B

cpg_site_names <- G[[1]]
G <- G[, -1, with = FALSE]
G <- as.data.frame(G)
rownames(G) <- cpg_site_names

scaled_B <- t(scale(t(B)))
df <- as.data.frame(scaled_B)

scaled_G <- t(scale(t(G)))
df2 <- as.data.frame(scaled_G)

orderr <- rownames(df2)
#df2 is test
df <- df[orderr,]

df <- na.omit(df)
orderr <- rownames(df)

df2 <- df2[orderr, ]

df_matrix <- as.matrix(df)
km.perm <- KMeansSparseCluster.permute(df_matrix,K=2,wbounds=seq(3,7,len=15),nperms=5)
plot(km.perm)
km.out <- KMeansSparseCluster(df_matrix, K=2,wbounds=km.perm$bestw)

# Seed for subset1
#set.seed(70)
# Random sampling for size reduction
#sampled_indices <- sample(nrow(B), size = 70000)
#subset1_B <- B[sampled_indices, ]
# Scaling
#subset1_B_scaled <- t(scale(t(subset1_B)))
#df <- as.data.frame(subset1_B_scaled)

#set.seed(123)

#elbow_method <- function(max_k) {
#  wss <- numeric(max_k)
#  for (k in 1:max_k) {
#    km <- kmeans(df, centers = k, nstart=25)
#    wss[k] <- sum(km$withinss)
#  }
#  plot(1:max_k, wss, type='b', xlab='Number of clusters', ylab='Within groups sum of squares')
#}

#elbow_method(10)

#optimal_clusters <- 3

#km <- kmeans(df, centers=optimal_clusters, nstart=25)
#df$cluster <- as.factor(km$cluster)
#df <- df[order(df$cluster), ]

#B <- B[final_order,]
#G <- G[final_order, ]

#png(filename="training2.png", width=10, height=8, units="in", res=300)
pheatmap(as.matrix(df[, -ncol(df)]),  # Exclude cluster column
         scale = "none",
         show_rownames = FALSE,
         show_colnames = TRUE,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize_col = 6,
         angle_col = 45,
         main="Testing",
         breaks = seq(-5, 5, length.out = 101)
         )
#dev.off()



#df1 <- df[which(df$cluster==1), ]

#age_values <- as.numeric(sub("age_", "", colnames(B)))

#df <- df[, -ncol(df)]
d <- dist(df, method = "euclidean")  # Compute the distance matrix
hc <- hclust(d, method = "ward.D2")  # Perform hierarchical clustering

# Plot the dendrogram
plot(hc, main = "Dendrogram", xlab = "Sample index", ylab = "Height")

num_clusters = 8 #change number
clusters <- cutree(hc, k = num_clusters)
df$Cluster <- as.factor(clusters)

df <- df[order(df$Cluster), ]
df_order <- rownames(df)
df2 <- df2[df_order, ]

library(splines)
library(dplyr)
library(tidyverse)

for (i in 1:num_clusters) {
  train_df <- df[df$Cluster == i, ]
  ident <- rownames(train_df)
  test_df <- df2[ident, ]
  
  long_df_test <- test_df %>%
    rownames_to_column(var = "Sample") %>%
    pivot_longer(cols = -c(Sample), names_to = "Age", values_to = "Beta") %>%
    mutate(Age = as.numeric(sub("age_", "", Age)))
  
  long_df_train <- train_df %>%
    rownames_to_column(var = "Sample") %>%
    pivot_longer(cols = -c(Sample, cluster), names_to = "Age", values_to = "Beta") %>%
    mutate(Age = as.numeric(sub("age_", "", Age)))
  
  #fit <- loess(Beta ~ Age, data = long_df, span = 0.3)
  fit <- lm(Beta ~ bs(Age, df = 5), data = long_df_train)
  fit2 <- lm(Beta ~ bs(Age, df = 5), data = long_df_test)
  long_df_test$Fitted <- predict(fit2)
  long_df_train$Fitted <- predict(fit)
  
  
  
  ggplot(long_df_test, aes(x = Age, y = Beta)) +
    geom_point(size=1) +
    geom_line(aes(y = Fitted), color = "blue") +
    labs(title = paste("Test Cluster", i), x = "Age", y = "Scaled Beta Values") +
    theme_minimal()
  ggsave(
    paste0("test_cluster_", i, ".png"),
    plot = last_plot(),
    bg="white",
    scale = 1,
    dpi = 300,
  )
  
  ggplot(long_df_train, aes(x = Age, y = Beta)) +
    geom_point(size=1) +
    geom_line(aes(y = Fitted), color = "blue") +
    labs(title = paste("Train Cluster", i), x = "Age", y = "Scaled Beta Values") +
    theme_minimal()
  ggsave(
    paste0("train_cluster_", i, ".png"),
    plot = last_plot(),
    bg="white",
    scale = 1,
    dpi = 300,
  )
}





# Reshape the data for ggplot
long_df <- df1 %>%
  rownames_to_column(var = "CpG") %>%
  pivot_longer(cols = -c(CpG, cluster), names_to = "Age", values_to = "Beta")

# Extract age information from the column names
long_df$Age <- as.numeric(sub("age_", "", long_df$Age))

library(npreg)

mod.ss <- ss(long_df$Beta, long_df$Age, nknots = 10)

loess_fit <- loess(Beta ~ Age, data = long_df, span = 0.3)


ggplot(long_df, aes(x=Age, y=Beta)) +
  geom_point(alpha=0.3) + 
  #geom_smooth(method="loess", se=FALSE)
  labs(x="Age", y="Scaled Beta value") + 
  theme_minimal() +
  xlim(0, 100)
ggsave(
  "filename",
  plot = last_plot(),
  scale = 1,
  dpi = 300,
)


