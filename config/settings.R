# PUBLIC DEMONSTRATION CONFIGURATION — fully synthetic input data.
# Replace paths and parameters after reading docs/SCIENTIFIC_METHODS.md.
# This workflow is for independent two-group comparisons only.
settings <- list(
  abundance_file = "data/example/protein_abundance.csv",
  metadata_file = "data/example/sample_metadata.csv",
  contrasts_file = "data/example/contrasts.csv",
  output_dir = "results",
  pseudocount = 1,
  fdr_threshold = 0.05,
  absolute_log2fc_threshold = 0.58,
  imputation_quantile = 0.01,
  jitter_fraction = 0.01,
  random_seed = 65,
  rng_mode = "independent",
  min_replicates_per_condition = 2,
  heatmap_enabled = TRUE,
  heatmap_top_n = 50,
  heatmap_row_scaling = TRUE,
  heatmap_cluster_rows = TRUE,
  heatmap_cluster_columns = TRUE
)
