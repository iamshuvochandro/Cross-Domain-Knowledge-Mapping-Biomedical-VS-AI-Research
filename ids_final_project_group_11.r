if (!require("pacman")) install.packages("pacman")
pacman::p_load(tm, umap, ggplot2, dplyr, GGally, scales, tidyr, cluster, mclust, aricode)

ds1 <- read.csv("Data/Biomedical_Diabetics_DS1.csv")
ds2 <- read.csv("Data/AI_Diabetics_DS2.csv")

ds1$domain <- 0
ds2$domain <- 1
merged_data <- rbind(ds1, ds2)

corpus <- Corpus(VectorSource(merged_data$Abstract))
corpus <- tm_map(corpus, content_transformer(tolower))
corpus <- tm_map(corpus, removePunctuation)
corpus <- tm_map(corpus, removeNumbers)
corpus <- tm_map(corpus, removeWords, stopwords("english"))
corpus <- tm_map(corpus, stripWhitespace)
merged_data$Cleaned_Abstract <- sapply(corpus, as.character)
ds <- merged_data[!is.na(merged_data$Cleaned_Abstract) & merged_data$Cleaned_Abstract != "", ]
write.csv(ds, "ids_final_dataset_group_11.csv", row.names = FALSE)

corpus_final <- Corpus(VectorSource(ds$Cleaned_Abstract))
dtm <- DocumentTermMatrix(corpus_final, control = list(weighting = weightTfIdf))
tfidf_matrix <- as.matrix(dtm)

set.seed(123)
umap_results <- umap(tfidf_matrix)
umap_df <- as.data.frame(umap_results$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$domain <- ds$domain
umap_df$domain_name <- factor(umap_df$domain, labels = c("Biomedical", "AI"))

set.seed(123)
clusters <- kmeans(umap_df[, c("UMAP1", "UMAP2")], centers = 5)
umap_df$Cluster <- as.factor(clusters$cluster)
ds$Cluster <- as.factor(clusters$cluster)

cluster_table <- table(umap_df$Cluster, umap_df$domain)
cluster_pct <- prop.table(cluster_table, 1) * 100
mapping_summary <- as.data.frame(cluster_pct)
colnames(mapping_summary) <- c("Cluster", "Domain", "Percentage")
mapping_summary <- pivot_wider(mapping_summary, names_from = Domain, values_from = Percentage)
colnames(mapping_summary) <- c("Cluster", "Bio_Pct", "AI_Pct")
mapping_summary$Topic_Type <- ifelse(mapping_summary$Bio_Pct > 70, "Biomedical Specific",
                                     ifelse(mapping_summary$AI_Pct > 70, "AI Specific", 
                                            "Shared / Cross-Domain"))

get_top_keywords <- function(matrix, clusters, n = 10) {
  top_keywords <- list()
  for (i in unique(clusters)) {
    cluster_matrix <- matrix[clusters == i, , drop = FALSE]
    avg_tfidf <- colMeans(cluster_matrix)
    top_words <- sort(avg_tfidf, decreasing = TRUE)[1:n]
    top_keywords[[as.character(i)]] <- names(top_words)
  }
  return(top_keywords)
}
cluster_keywords <- get_top_keywords(tfidf_matrix, ds$Cluster)

sil <- silhouette(clusters$cluster, dist(umap_results$layout))
avg_sil <- mean(sil[, 3])
ari_score <- adjustedRandIndex(ds$domain, clusters$cluster)
nmi_score <- NMI(ds$domain, clusters$cluster)

print(mapping_summary)
print(cluster_keywords)
print(paste("Average Silhouette Score:", round(avg_sil, 4)))
print(paste("Adjusted Rand Index (ARI):", round(ari_score, 4)))
print(paste("Normalized Mutual Information (NMI):", round(nmi_score, 4)))

find_hull <- function(df) df[chull(df$UMAP1, df$UMAP2), ]
hulls <- umap_df %>% group_by(Cluster) %>% do(find_hull(.))

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(aes(shape = domain_name), size = 2, alpha = 0.6) +
  geom_polygon(data = hulls, aes(fill = Cluster), alpha = 0.1) +
  theme_minimal() +
  labs(title = "Cross-Domain Knowledge Map: Cluster Boundaries",
       shape = "Research Domain", color = "Research Topic (Cluster)")

ggplot(umap_df, aes(x = UMAP1, fill = Cluster)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  theme_minimal() +
  labs(title = "Distribution of UMAP Dimension 1", x = "UMAP1 Value", y = "Frequency")

ggplot(umap_df, aes(x = Cluster, fill = domain_name)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = percent) +
  theme_minimal() +
  labs(title = "Domain Composition of Clusters", y = "Percentage (%)", fill = "Domain")

ggplot(umap_df, aes(x = Cluster, y = UMAP2, fill = Cluster)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Variance of UMAP Dimension 2 by Cluster")

ggpairs(umap_df, columns = c("UMAP1", "UMAP2"), aes(color = Cluster, alpha = 0.5)) +
  theme_minimal() +
  labs(title = "Scatter Matrix of UMAP Dimensions")

ggplot(umap_df, aes(x = domain_name, y = UMAP1, fill = domain_name)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  theme_minimal() +
  labs(title = "Density Distribution of UMAP Dimensions by Domain")