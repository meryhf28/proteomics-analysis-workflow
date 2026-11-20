# OPTIONAL PEAKS preparation settings. Supply a manifest yourself.
# All thresholds are USER-DEFINED, not universal proteomics recommendations.
peaks <- list(
  manifest_file = "data/peaks_manifest.csv",
  output_dir = "data/processed",
  min_unique = 2,
  min_score = 25,
  min_present_replicates = 2,
  # Example: "keratin|trypsin|contaminant" — adjust to the experiment.
  exclude_description_regex = ""
)
