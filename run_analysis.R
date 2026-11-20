# ==============================================================
# PROTEOMICS ANALYSIS WORKFLOW — RUN FROM REPOSITORY ROOT
# ==============================================================
# 1. Edit config/settings.R to point to your own input tables.
# 2. Prepare three CSVs as documented in README.md.
# 3. Install dependencies: Rscript scripts/install_dependencies.R
# 4. Run: Rscript run_analysis.R
# Every comparison produces Excel/CSV tables and a volcano PNG.
# This generic core starts from quantified protein abundance matrices.
# PEAKS filtering, FASTA retrieval and enrichment are separate extensions.
source("config/settings.R")
source("R/functions.R")
source("R/visualization.R")
required <- c("limma", "ggplot2", "openxlsx")
for (pkg in required) if (!requireNamespace(pkg, quietly=TRUE)) stop("Missing R package: ",pkg,". Run scripts/install_dependencies.R")
for (p in c(settings$abundance_file, settings$metadata_file, settings$contrasts_file)) if (!file.exists(p)) stop("Input file not found: ",p)
abundance <- utils::read.csv(settings$abundance_file, check.names=FALSE, stringsAsFactors=FALSE)
metadata <- utils::read.csv(settings$metadata_file, check.names=FALSE, stringsAsFactors=FALSE)
# Protect numeric sample IDs (5–16) from conversion to positional indices.
metadata$sample_id <- as.character(metadata$sample_id)
contrasts <- utils::read.csv(settings$contrasts_file, check.names=FALSE, stringsAsFactors=FALSE)
validate_tables(abundance, metadata, contrasts, settings$min_replicates_per_condition)
dir.create(settings$output_dir, recursive=TRUE, showWarnings=FALSE)
# The legacy_sequential option initializes RNG once across all contrasts;
# in this mode, contrast order changes the imputation random stream.
if (identical(settings$rng_mode, "legacy_sequential")) {
  set.seed(settings$random_seed)
  message("RNG: legacy sequential stream, seed ", settings$random_seed)
} else {
  message("RNG: independent deterministic seed per contrast")
}
summary <- vector("list", nrow(contrasts))
for (i in seq_len(nrow(contrasts))) {
  spec <- contrasts[i, ]
  message("Processing: ",spec$organism," — ",spec$treatment," vs ",spec$reference)
  result <- analyze_one_contrast(abundance, metadata, spec, settings)
  target <- save_outputs(result, spec, settings$output_dir, settings)
  if (isTRUE(settings$heatmap_enabled)) {
    # All significant proteins and the top N by adjusted p-value are plotted.
    for(mode in c("top_n","all")) {
      save_heatmap(result,target,paste(spec$organism,spec$treatment,"vs",spec$reference),
          mode=mode,n=settings$heatmap_top_n,scale_rows=settings$heatmap_row_scaling,
          cluster_rows=settings$heatmap_cluster_rows,
          cluster_columns=settings$heatmap_cluster_columns)
    }
  }
  summary[[i]] <- data.frame(organism=spec$organism, treatment=spec$treatment,
                             reference=spec$reference, tested=nrow(result$table),
                             significant=sum(result$table$significant), output=target)
}
utils::write.csv(do.call(rbind,summary),file.path(settings$output_dir,"analysis_summary.csv"),row.names=FALSE)
writeLines(capture.output(utils::sessionInfo()),file.path(settings$output_dir,"sessionInfo.txt"))
message("Completed. Results in: ", settings$output_dir)
