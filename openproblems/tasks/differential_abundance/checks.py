import numpy as np


def check_dataset(adata):
    """Check that dataset output fits expected API."""
    assert "condition" in adata.obs
    assert "replicate" in adata.obs
    assert "sample" in adata.obs
    assert "ground_truth_probability" in adata.obsm

    assert "conditions" in adata.uns
    assert "replicates" in adata.uns
    assert "samples" in adata.uns
    assert "n_conditions" in adata.uns
    assert "n_replicates" in adata.uns
    assert "n_samples" in adata.uns

    # Ensure that ground truth probability is infact a probability distirbution
    assert np.all(np.min(adata.obsm["ground_truth_probability"]) >= 0)
    assert np.all(np.max(adata.obsm["ground_truth_probability"]) <= 1)
    assert np.allclose(np.sum(adata.obsm["ground_truth_probability"], axis=1), 1)

    # Check that sample labels are formatted properly
    for sample in adata.uns["samples"]:
        cond, rep = sample.split(".")
        assert cond in adata.uns["conditions"]
        assert rep in adata.uns["replicates"]

    return True


def check_method(adata):
    """Check that method output fits expected API."""
    # TODO: update
    assert "probability_estimate" in adata.obs
    return True
