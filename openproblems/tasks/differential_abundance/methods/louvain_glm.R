# Milo
# @param sce SingleCellExperiment
# @return sce SingleCellExperiment

# Dependencies
library(SingleCellExperiment)
library(reshape2)
library(tidyverse)

# Make design matrix
sample_col <- "sample"
condition_col <- "condition"
design_df <- as.tibble(colData(sce)[c(sample_col, condition_col)]) %>%
  distinct() %>%
  dplyr::rename(sample = sample_col)
design <- formula(paste("Freq ~", condition_col, "+ offset(log(N_s))", collapse = " "))

condition_vec <- colData(sce)[[condition_col]]
sample_labels <- colData(sce)[[sample_col]]
clust.df <- data.frame("cell_id" = colnames(sce), "Louvain.Clust" = colData(sce)[["louvain"]])

clust.df$Sample <- sample_labels
clust.df$Condition <- condition_vec

louvain.count <- table(clust.df$Louvain.Clust, clust.df$Sample)
attributes(louvain.count)$class <- "matrix"

## Run test

df <- melt(louvain.count, varnames = c("cluster", "sample"), value.name = "Freq") %>%
  mutate(cluster = factor(cluster)) %>%
  left_join(design_df, by = "sample") %>%
  group_by(sample) %>%
  mutate(N_s = sum(Freq)) %>%
  ungroup() %>%
  group_by(cluster) %>%
  do(model = glm(design, data = ., family = "poisson"))

res_df <- t(sapply(df$model, function(x) summary(x)$coefficients[nrow(summary(x)$coefficients), ]))
res_df <- data.frame(res_df)
colnames(res_df) <- c("logFC", "Std. Error", "z value", "Pval")
louvain.res <- cbind(df, res_df) %>%
  mutate(FDR = p.adjust(Pval, method = "BH"))
rownames(louvain.res) <- louvain.res$cluster

clust.df$logFC <- louvain.res[clust.df$Louvain.Clust, "logFC"]
clust.df$FDR <- louvain.res[clust.df$Louvain.Clust, "FDR"]

## Convert cluster FC to single-cell probability
cell_p_cond1 <- exp(clust.df$logFC) / (exp(clust.df$logFC) + 1)
cell_p_cond2 <- 1 - cell_p_cond1
prob_mat <- cbind(cell_p_cond1, cell_p_cond2)
reducedDim(sce, "probability_estimate") <- as.matrix(prob_mat)

## Return
sce
