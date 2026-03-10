############################################################
# 03 Technical replicate QC
############################################################

library(ggplot2)
library(ggrepel)

log_lfq_matrix = readRDS("../data/log_lfq_matrix.rds")
metadata = read.csv("../data/sample_metadata.csv")

# Correlation
cor_matrix =
  cor(log_lfq_matrix, use = "pairwise.complete.obs")

write.csv(cor_matrix,
          "../results/technical_correlation_matrix.csv")

# PCA
complete_idx = rowSums(is.na(log_lfq_matrix)) == 0

pca_matrix =
  log_lfq_matrix[complete_idx, ]

pca =
  prcomp(t(pca_matrix), scale = TRUE)

pca_df = as.data.frame(pca$x)

pca_df$Group = metadata$Group

ggplot(pca_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 4) +
  theme_minimal()

ggsave("../results/pca_quality_control.png")
