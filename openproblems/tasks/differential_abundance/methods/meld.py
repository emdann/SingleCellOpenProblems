from ....tools.decorators import method
from ....tools.utils import check_version

import meld
import numpy as np
import pandas as pd
import scanpy as sc


@method(
    method_name="MELD",
    paper_name="Quantifying the effect of experimental perturbations in "
    "single-cell data",
    paper_url="https://www.biorxiv.org/content/10.1101/532846v4",
    paper_year=2020,
    code_url="https://github.com/krishnaswamylab/MELD",
    code_version=check_version("meld"),
    # image="openproblems-template-image" # only if required
)
def run_meld(adata):
    # Library size normalize and sqrt transform adata
    result_dict = sc.pp.normalize_total(adata, target_sum=1e5, inplace=False)
    adata.layers["X_norm"] = result_dict["X"]
    adata.obs["norm_factor"] = result_dict["norm_factor"]
    adata.layers["X_norm_sqrt"] = adata.layers["X_norm"].sqrt()

    # Complete the result in-place
    meld_op = meld.MELD(verbose=False)
    adata.obsm["sample_densities"] = meld_op.fit_transform(
        adata.layers["X_norm_sqrt"], sample_labels=adata.obs["sample"]
    ).set_index(adata.obs_names)

    # Normalize the probability estimates for each condition per replicate
    adata.obsm["probability_estimate"] = pd.DataFrame(
        np.zeros(shape=(adata.n_obs, adata.uns["n_conditions"])),
        index=adata.obs_names,
        columns=adata.uns["conditions"],
    )

    for replicate in adata.uns["replicates"]:
        replicate_sample_mask = [
            replicate in sample for sample in adata.obsm["sample_densities"].columns
        ]
        curr_densities = adata.obsm["sample_densities"].loc[:, replicate_sample_mask]
        likelihoods = meld.utils.normalize_densities(curr_densities)
        likelihoods.columns = [col.split(".")[0] for col in likelihoods.columns]
        adata.obsm["probability_estimate"] += likelihoods

    # Calculate the average
    adata.obsm["probability_estimate"] /= adata.uns["n_replicates"]
    return adata
