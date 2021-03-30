from ....tools.decorators import metric

import sklearn.metrics


@metric(
    metric_name="MSE",
    maximize=False,
    # image="openproblems-template-image" # only if required
)
def mean_squared_error(adata):
    # TODO: update
    return sklearn.metrics.mean_squared_error(
        adata.obsm["ground_truth_probability"], adata.obsm["probability_estimate"]
    )
