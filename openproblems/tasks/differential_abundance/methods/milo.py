from ....tools.conversion import r_function
from ....tools.decorators import method
from ....tools.utils import check_version

_milo = r_function("milo.R")


@method(
    method_name="Milo",
    paper_name="Milo: differential abundance testing on single-cell data using k-NN graphs",
    paper_url="https://www.biorxiv.org/content/10.1101/2020.11.23.393769v1",
    paper_year=2020,
    code_url="https://github.com/MarioniLab/miloR",
    code_version=check_version("rpy2"),
    image="openproblems-r-extras",
)
def run_milo(adata):
    # Prepare anndata for anndata2ri conversion
    adata.obsm["ground_truth_probability"] = adata.obsm[
        "ground_truth_probability"
    ].values
    adata = _milo(adata)
    # Reconvert probabilities to DataFrames
    adata.obsm["ground_truth_probability"] = pd.DataFrame(
        adata.obsm["ground_truth_probability"],
        index=adata.obs_names,
        columns=adata.uns["conditions"],
    )
    adata.obsm["probability_estimate"] = pd.DataFrame(
        adata.obsm["probability_estimate"],
        index=adata.obs_names,
        columns=adata.uns["conditions"],
    )
    return adata
