# Milo
# @param sce SingleCellExperiment
# @return sce SingleCellExperiment

# Dependencies
library(SingleCellExperiment)
library(miloR)
library(tidyverse)
library(Matrix)

# Make design matrix
sample_col = "sample"
condition_col = "condition"
design_df <- as.tibble(colData(sce)[c(sample_col, condition_col)]) %>%
    distinct() %>%
    column_to_rownames(sample_col)
design <- formula(paste('~', condition_col, collapse = ' '))

## Build graph neighbourhoods
k <- metadata(sce)$n_samples * 5
d <- 30
prop <- 0.1
milo <- Milo(sce)
milo <- buildGraph(milo, k=k, d=d, reduced.dim = "PCA")
milo <- makeNhoods(milo, prop = prop, k=k, d=d, reduced_dims = "PCA")

## Test DA
milo <- countCells(milo, meta.data = data.frame(colData(milo)), sample=sample_col)
milo <- calcNhoodDistance(milo, d=d, reduced.dim = "PCA")
DA_results <- testNhoods(milo, design = design, design.df = design_df)

## Convert nhood FC to single-cell probability
nhood_p <- exp(DA_results$logFC) / (exp(DA_results$logFC) + 1)
cell_p_cond1 <- milo@nhoods %*% nhood_p
cell_p_cond1 <- cell_p_cond1 / rowSums(milo@nhoods)
cell_p_cond2 <- 1 - cell_p_cond1
prob_mat <- cbind(cell_p_cond1, cell_p_cond2)
reducedDim(sce, "probability_estimate") <- as.matrix(prob_mat)

## Return
sce
