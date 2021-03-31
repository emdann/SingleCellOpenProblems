from ....tools.conversion import r_function
from ....tools.decorators import method


@method(
    method_name="DA-seq",
    paper_name="Detection of differentially abundant cell"
    "subpopulations discriminates biological states"
    "in scRNA-seq data"
    paper_url="https://www.biorxiv.org/content/10.1101/711929v3",
    paper_year=2020,
    code_url="https://github.com/KlugerLab/DAseq",
    #code_version="N/A",
    # image="openproblems-template-image" # only if required
)

_daseq = r_function("DAseq.R")

def DAseq(adata):
  adata = _daseq(adata)
  return adata
