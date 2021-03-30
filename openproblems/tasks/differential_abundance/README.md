<!--- TODO: update --->

# Differential abundance (single cell condition likelihood estimation)

The goal of this task is to estimate the likelihood of a sample condition label given a cell's UMI count.

This task is calculated using experimental scRNA-seq data and then artificially creating sample labels over the dataset. The labels are created using a ground-truth probability function that is smooth over the dataset. Labels indicate `condition` and `replicate`. In the current task definition, `condition` probabilities are different for each condition across the data, but `replicate` assignments are uniform random across the dataset. We then concatenate the `condition` and `replicate` labels to create a `sample` label, e.g. `condition2.replicate3`.

Ground-truth probabilities are created using the `differential_abundance.dataset.utils.simulate_treatment()` function. Briefly, to create ground-truth probabilities that a given cell has a given condition label (i.e. was observed in a given condition), we take a convex combination of principal components of the data for each condition. These PC-combinations are then L1-normalized across samples to create a probability of each label given a cell (the probabilities of all the labels sums to 1 across conditions). Cells are then randomly assigned to each replicate.

## API

Datasets should contain the following attributes:

* `adata.obs["condition"]: array-like, shape=(n_obs,)`
  * indicates the condition for each cell. E.g. `condition4`
* `adata.obs["replicate"]: array-like, shape=(n_obs,)`
  * indicates the replicate for each cell. E.g. `replicate2`
* `adata.obs["sample"]: array-like, shape=(n_obs,)`
  * a concatenation of the condition and replicate for each cell. E.g. `condition3.replicate1`
* `adata.obsm["ground_truth_probability"]: array-like, shape=(n_obs,n_conditions)`
  * probability that each cell would be assigned a given condition label
* `adata.uns["conditions"]: array-like, shape=(n_conditions,)`
  * list of unique condition labels
* `adata.uns["replicates"]: array-like, shape=(n_replicates,)`
  * list of unique replicate labels
* `adata.uns["samples"]: array-like, shape=(n_conditions*n_replicates,)`
  * List of unique sample labels
* `adata.uns["n_conditions"]: int`
* `adata.uns["n_replicates"]: int`
* `adata.uns["n_samples"]: int`

Methods output:

* `adata.obs['probability_estimate']: array-like, shape=(n_obs, n_conditions)`
  * Estimate of the probability that each cell would have a given condition label

Metric definition:

Metrics should assess the similarity between `adata['ground_truth_probability']` and `adata.obs['probability_estimate']`.
