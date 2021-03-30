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

@dataset("Chd/tyr CRISPR perturbation dataset - multiple simulations")
def Klein_2018_chd_tyr_data_n_simulations(test=False, n_simulations=10):
    # Load UMI data
    adata = load_zebrafish_chd_tyr(test=test)
    # Simulate experimental conditions with N different seeds
    for i in range(n_simulations):
        seed = 0 + i
        if seed < 5:
            effect_size = 0.6
        else:
            effect_size = 0.8
        simulate_treatment(adata, seed=seed, n_conditions=2, n_replicates=3, effect_size=effect_size)
        obs_simulation = ["condition", "replicate", "sample"]
        adata.obs.columns = ['{x}_seed{s}'.format(x=x, s=seed) if x in obs_simulation else x for x in adata.obs.columns]
        adata.obsm["ground_truth_probability_seed{s}".format(s=seed)] = adata.obsm["ground_truth_probability"]
        uns_simulation = ["conditions", "replicates", "samples", "n_replicates", "n_conditions", "n_samples"]
        for u in uns_simulation:
            adata.uns['{u}_seed{s}'.format(u=u, s=seed)] = adata.uns[u]
        ## Save simulation params used for metrics
        adata.uns["DA_simulation_params_seed{s}".format(s=seed)] = {"max_effect_size":effect_size}
    return adata