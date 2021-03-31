from ....tools.conversion import r_function
from ....tools.decorators import method
from ....tools.utils import check_version

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
    adata = _louvain_glm(adata)
    return adata
