# General-purpose proteomics differential-abundance analysis.
# Core method: log2(x+1), column median centring, zeros treated as missing,
# per-column low-quantile + MAD-based Gaussian jitter, and limma.
# It does NOT claim equivalence to imputeLCMD::impute.MinDet.

validate_tables <- function(abundance, metadata, contrasts, min_reps) {
  required_abundance <- c("accession", "description")
  required_metadata <- c("sample_id", "organism", "condition")
  required_contrasts <- c("organism", "reference", "treatment")
  for (x in required_abundance) if (!x %in% names(abundance)) stop("Missing abundance column: ", x)
  for (x in required_metadata) if (!x %in% names(metadata)) stop("Missing metadata column: ", x)
  for (x in required_contrasts) if (!x %in% names(contrasts)) stop("Missing contrast column: ", x)
  if (anyNA(abundance$accession) || anyDuplicated(abundance$accession) || any(!nzchar(trimws(as.character(abundance$accession))))) stop("Protein accessions must be unique and nonmissing")
  if (anyNA(metadata[, required_metadata, drop=FALSE]) || anyDuplicated(metadata$sample_id)) stop("Metadata must contain unique, nonmissing sample IDs")
  if (any(!nzchar(trimws(as.character(metadata$sample_id))))) stop("Empty sample IDs are not permitted")
  if (nrow(contrasts) == 0L) stop("At least one contrast is required")
  if (anyNA(contrasts)) stop("Contrast definitions cannot contain missing values")
  if (any(!nzchar(trimws(as.character(contrasts$organism)))) || any(!nzchar(trimws(as.character(contrasts$reference)))) || any(!nzchar(trimws(as.character(contrasts$treatment))))) stop("Empty contrast labels are not allowed")
  if (any(contrasts$reference == contrasts$treatment)) stop("Reference and treatment must differ")
  if (anyDuplicated(paste(contrasts$organism, contrasts$reference, contrasts$treatment, sep="|"))) stop("Duplicate contrasts")
  if (!all(metadata$sample_id %in% names(abundance))) stop("Missing intensity columns: ", paste(setdiff(metadata$sample_id, names(abundance)), collapse=", "))
  if (!all(names(abundance) %in% c(required_abundance, metadata$sample_id))) stop("Unexpected abundance columns; only accession, description and sample IDs are supported")
  nums <- abundance[, metadata$sample_id, drop=FALSE]
  for (name in names(nums)) if (!is.numeric(nums[[name]]) || any(!is.finite(nums[[name]]) & !is.na(nums[[name]])) || any(nums[[name]] < 0, na.rm=TRUE)) stop("Abundances must be finite, nonnegative numeric values or NA: ", name)
  for (k in seq_len(nrow(contrasts))) {
    z <- contrasts[k, ]; md <- metadata[metadata$organism == z$organism & metadata$condition %in% c(z$reference, z$treatment), ]
    if (sum(md$condition == z$reference) < min_reps || sum(md$condition == z$treatment) < min_reps) stop("Insufficient replicates for ", z$organism, " / ", z$treatment, " vs ", z$reference)
  }
  invisible(TRUE)
}

low_quantile_impute <- function(x, groups, q = .01, jitter_fraction = .01) {
  # GROUPS select the columns, but the 1% quantile and jitter are calculated
  # independently within each sample.
  out <- x
  for (g in unique(groups)) for (j in which(groups == g)) {
    v <- out[, j]
    if (all(is.na(v))) next
    threshold <- as.numeric(stats::quantile(v, probs=q, na.rm=TRUE, names=FALSE))
    spread <- stats::mad(v, center=stats::median(v, na.rm=TRUE), na.rm=TRUE)
    missing <- is.na(v)
    if (any(missing)) v[missing] <- threshold + stats::rnorm(sum(missing), sd=max(jitter_fraction * spread, 1e-6))
    out[, j] <- v
  }
  out
}

analyze_one_contrast <- function(abundance, metadata, spec, settings) {
  md <- metadata[metadata$organism == spec$organism & metadata$condition %in% c(spec$reference, spec$treatment), , drop=FALSE]
  samples <- md$sample_id
  # Only analyse the requested organism and two conditions. Sample groupings
  # come from metadata, not hard-coded PEAKS run numbers or sample-name regex.
  x <- as.matrix(abundance[, samples, drop=FALSE]); storage.mode(x) <- "double"
  colnames(x) <- samples
  raw_zero <- !is.na(x) & x == 0
  transformed <- log2(x + settings$pseudocount)
  medians <- apply(transformed, 2, stats::median, na.rm=TRUE)
  if (any(!is.finite(medians))) stop("A sample has no observed intensities: ", paste(samples[!is.finite(medians)], collapse=", "))
  # Median centring removes a per-run location shift on the log scale.
  # Do not assume it corrects batch effects or changes in global proteome abundance.
  centered <- sweep(transformed, 2, medians, "-")
  centered[raw_zero] <- NA_real_
  missing_pre <- rowSums(is.na(centered))
  # Exclude proteins absent in every sample of this contrast; imputation is
  # impossible in that case and no statistical evidence exists.
  retain <- rowSums(!is.na(centered)) > 0
  if (sum(retain) < 2) stop("Too few proteins with measured abundance in contrast")
  # Reproducible RNG strategy:
  #   legacy_sequential: seed once in run_analysis.R, then use one continuous
  #       RNG stream across contrasts in the specified order. This mirrors
  #       an order-dependent sequence; use only for compatibility checks.
  #   independent: derive a deterministic seed per contrast, so reordering
  #       contrasts does not change the imputation of any given comparison.
  if (identical(settings$rng_mode, "independent")) {
    contrast_key <- paste(spec$organism, spec$reference, spec$treatment, sep="|")
    # Small bounded integer hash, not R string hashing or global RNG state.
    hash <- 0
    for (code in utf8ToInt(enc2utf8(contrast_key)))
      hash <- (hash * 131 + as.double(code)) %% 2147483646
    set.seed(as.integer((as.double(settings$random_seed) + hash) %% 2147483646 + 1))
  } else if (!identical(settings$rng_mode, "legacy_sequential")) {
    stop("Unknown rng_mode: ", settings$rng_mode)
  }
  imp <- low_quantile_impute(centered[retain, , drop=FALSE], md$condition, settings$imputation_quantile, settings$jitter_fraction)
  if (anyNA(imp)) stop("Unresolved missing values after imputation; inspect all-missing samples")
  rownames(imp) <- abundance$accession[retain]
  groups <- factor(md$condition, levels=c(spec$reference, spec$treatment))
  design <- stats::model.matrix(~ 0 + groups)
  colnames(design) <- c("reference", "treatment")
  if (qr(design)$rank < 2 || ncol(imp) <= 2) stop("Not enough samples to estimate differential abundance")
  # limma fits protein-specific linear models, then moderates residual variances
  # through empirical Bayes shrinkage. It does not replace biological design QC.
  fit <- limma::lmFit(imp, design)
  contrast <- limma::makeContrasts(treatment_minus_reference = treatment - reference, levels=design)
  fit <- limma::eBayes(limma::contrasts.fit(fit, contrast))
  tt <- limma::topTable(fit, coef="treatment_minus_reference", number=Inf, sort.by="none")
  tt$accession <- rownames(tt)
  # Preserve original accessions and prevent implicit reordering by merge().
  match_index <- match(abundance$accession[retain], tt$accession)
  if (anyNA(match_index)) stop("limma results missing protein accession(s)")
  output <- cbind(as.data.frame(abundance[retain, c("accession", "description")]),
                  tt[match_index, setdiff(names(tt), "accession"), drop=FALSE])
  output$log2FC <- output$logFC
  output$Pvalue <- output$P.Value
  output$FDR <- output$adj.P.Val
  output$significant <- !is.na(output$FDR) & output$FDR <= settings$fdr_threshold & abs(output$log2FC) >= settings$absolute_log2fc_threshold
  output$direction <- ifelse(!output$significant, "NS", ifelse(output$log2FC > 0, "Higher in treatment", "Higher in reference"))
  output$n_imputed <- missing_pre[retain][match(output$accession, abundance$accession[retain])] - rowSums(is.na(imp))[match(output$accession, rownames(imp))]
  output <- output[order(output$FDR, -abs(output$log2FC), na.last=TRUE), ]
  list(table=output, matrix=imp, metadata=md, omitted_accessions=abundance$accession[!retain])
}

save_outputs <- function(result, spec, output_dir, settings) {
  label <- paste(spec$organism, spec$treatment, "vs", spec$reference, sep="_")
  safe <- gsub("[^A-Za-z0-9_-]", "_", label)
  dir <- file.path(output_dir, safe); dir.create(dir, recursive=TRUE, showWarnings=FALSE)
  utils::write.csv(result$table, file.path(dir, "differential_abundance.csv"), row.names=FALSE)
  utils::write.csv(result$table[result$table$significant, , drop=FALSE], file.path(dir, "significant_proteins.csv"), row.names=FALSE)
  utils::write.csv(cbind(accession=rownames(result$matrix), as.data.frame(result$matrix)), file.path(dir,"imputed_matrix.csv"), row.names=FALSE)
  utils::write.csv(result$metadata, file.path(dir,"sample_metadata_used.csv"), row.names=FALSE)
  writeLines(result$omitted_accessions, file.path(dir,"omitted_all_missing_accessions.txt"))
  # Volcano display: x is the limma log2 fold change; y is -log10(BH FDR).
  # Points are classified using BOTH the FDR and fold-change thresholds.
  # A horizontal FDR line is therefore correctly aligned to the plotted metric.
  df <- result$table
  p <- make_volcano(df, spec, settings)
  # The horizontal guide represents the FDR cutoff on this adjusted-p axis.
  ggplot2::ggsave(file.path(dir,"volcano.png"), p, width=7, height=5, dpi=300)
  wb <- openxlsx::createWorkbook()
  for (sheet in c("All_proteins","Significant","Imputed_matrix")) openxlsx::addWorksheet(wb,sheet)
  openxlsx::writeData(wb,"All_proteins",result$table)
  openxlsx::writeData(wb,"Significant",result$table[result$table$significant, , drop=FALSE])
  openxlsx::writeData(wb,"Imputed_matrix",cbind(accession=rownames(result$matrix), as.data.frame(result$matrix)))
  openxlsx::saveWorkbook(wb,file.path(dir,"analysis_results.xlsx"),overwrite=TRUE)
  invisible(dir)
}
