from ....data.Wagner_2018_zebrafish_embryo_CRISPR import load_zebrafish_chd_tyr
from ....tools.decorators import dataset
from .utils import simulate_treatment


@dataset("Chd/tyr CRISPR perturbation dataset")
def Wagner_2018_two_condition(test=False):
    # Load UMI data
    adata = load_zebrafish_chd_tyr(test=test)
    # Simulate experiment as a combination of PC dimensions
    simulate_treatment(adata, n_conditions=3, n_replicates=2)

    return adata


@dataset("Chd/tyr CRISPR perturbation dataset - multiple simulations")
def Wagner_2018_chd_tyr_data_n_simulations(test=False, n_simulations=20):
    # Load UMI data
    adata = load_zebrafish_chd_tyr(test=test)

    # Simulate experimental conditions with N different seeds
    # Determine which seed corresponds to which effect size
    seeds = np.arange(n_simulations)
    effect_size = {}
    for seed in seeds:
        if seed < 10:
            effect_size[seed] = 0.6
        elif seed < 20:
            effect_size = 1

    for seed in range(seeds):
        curr_effect_size = effect_size[seed]
        simulate_treatment(
            adata,
            seed=seed,
            n_conditions=2,
            n_replicates=3,
            effect_size=curr_effect_size)
        obs_simulation = ["condition", "replicate", "sample"]
        adata.obs.columns = [
            "{x}_seed{s}".format(x=x, s=seed) if x in obs_simulation else x
            for x in adata.obs.columns
        ]
        adata.obsm["ground_truth_probability_seed{s}".format(s=seed)] = adata.obsm[
            "ground_truth_probability"
        ]
        uns_simulation = [
            "conditions",
            "replicates",
            "samples",
            "n_replicates",
            "n_conditions",
            "n_samples",
        ]
        for u in uns_simulation:
            adata.uns["{u}_seed{s}".format(u=u, s=seed)] = adata.uns[u]
        # Save simulation params used for metrics
        adata.uns["DA_simulation_params_seed{s}".format(s=seed)] = {
            "max_effect_size": effect_size
        }
    return adata
