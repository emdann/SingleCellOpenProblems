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
    # Save the map from seed to effect size in `uns`
    adata.uns['seed_info'] = {}
    for seed in seeds:
        if seed < 10:
            # This is a small effect size
            adata.uns['seed_info'][seed] = {'effect_size' : 0.6}
        elif seed < 20:
            # This is a large effect size
            adata.uns['seed_info'][seed] = {'effect_size' : 1 }

    # Iterate through seeds
    for seed in range(seeds):
        # Simulate a treatment effect across the dataset
        simulate_treatment(
            adata,
            seed=seed,
            n_conditions=2,
            n_replicates=3,
            effect_size=adata.uns['seed_info'][seed]['effect_size'],
            seed=seed)

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
        max_effect_size = (
            (adata.obsm["ground_truth_probability_seed{s}".format(s=seed)] - 0.5)
            .max()
            .max()
        )
        adata.uns["DA_simulation_params_seed{s}".format(s=seed)] = {
            "max_effect_size": max_effect_size
        }
    return adata
