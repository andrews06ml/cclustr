# cclustr 0.1.2

## Maintenance

* Updated the package vignette to use the `SynthDiabetes2` dataset from
  `mlbench` (>= 2.1-11), replacing `PimaIndiansDiabetes2`, which was
  removed from `mlbench` and is no longer available. Added
  `mlbench (>= 2.1-11)` to `Suggests` in DESCRIPTION.

# cclustr 0.1.1

* Initial release on CRAN.

## New features

* `as_mild_list()`: standardizes multiple imputation outputs into a 
  unified list of completed datasets.
* `cluster_imputations()`: applies clustering algorithms across all 
  imputed datasets.
* `consensus_clustering()`: builds a consensus partition via 
  co-assignment matrix.
* `validate_clustering()`: computes internal and stability validation 
  metrics.
* `choose_best_clustering()`: selects the optimal k via weighted rank 
  aggregation.
* `run_mi_clustering()`: convenience wrapper for the full pipeline.
* `plot_consensus_matrix()`: heatmap of the co-assignment matrix.
* `plot_consensus_dendrogram()`: dendrogram of the consensus solution.
* `plot_validation_metrics()`: multi-panel plot of validation metrics.

