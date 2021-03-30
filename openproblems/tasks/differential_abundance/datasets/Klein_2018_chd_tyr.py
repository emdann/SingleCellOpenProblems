from ....data.Klein_2018_zebrafish_embryo import load_zebrafish_chd_tyr
from ....tools.decorators import dataset
from .utils import simulate_treatment


@dataset("Chd/tyr CRISPR perturbation dataset")
def Klein_2018_chd_tyr_data(test=False):
    # Load UMI data
    adata = load_zebrafish_chd_tyr(test=test)
    # Simulate experiment as a combination of PC dimensions
    simulate_treatment(adata, n_conditions=3, n_replicates=2)

    return adata
