from ....tools.conversion import r_function
from ....tools.decorators import method
from ....tools.utils import check_version

import scanpy as sc
import pandas as pd

_louvain_glm = r_function("louvain_glm.R")


@method(
    method_name="Poisson GLM on Louvain clustering",
    paper_name="",
    paper_url="https://www.nature.com/articles/nature24489",
    paper_year=2017,
    code_url="",
    code_version=check_version("rpy2"),
    image="openproblems-r-extras",
)
def run_louvain_glm(adata):
    # Prepare anndata for anndata2ri conversion
    adata.obsm["ground_truth_probability"] = adata.obsm[
        "ground_truth_probability"
    ].values
    # Run louvain clustering
    k = adata.uns["n_samples"] * 5
    sc.pp.neighbors(adata, n_neighbors=k)
    sc.tl.louvain(adata)
    # Test with Poisson GLM
    adata = _louvain_glm(adata)
    # Reconvert probabilities to DataFrames
    adata.obsm["ground_truth_probability"] = pd.DataFrame(adata.obsm["ground_truth_probability"], index=adata.obs_names, columns=adata.uns["conditions"])
    adata.obsm["probability_estimate"] = pd.DataFrame(adata.obsm["probability_estimate"], index=adata.obs_names, columns=adata.uns["conditions"])
    return adata
