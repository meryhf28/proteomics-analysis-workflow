# Parse-check every public .R script before running any tests.
code_files <- list.files(".", pattern="\\.R$", recursive=TRUE, full.names=TRUE)
for (f in code_files) tryCatch(parse(file=f), error=function(e) stop("Syntax failure in ", f, ": ", conditionMessage(e)))
message("PASS: R parser accepted ",length(code_files)," files.")
# Run from repository root: source("tests/run_tests.R")
# Only uses synthetic demonstration data; no user study files are present.
source("R/functions.R")
# The result exporter uses make_volcano(), defined in R/visualization.R.
# Load both modules explicitly so this test also works in a fresh R session.
source("R/visualization.R")
source("config/settings.R")
source("scripts/self_check.R")
stopifnot(is.function(make_volcano), is.function(save_heatmap))
if (!requireNamespace("limma", quietly=TRUE)) stop("Install limma before tests")
if (!requireNamespace("ggplot2", quietly=TRUE)) stop("Install ggplot2 before tests")
if (!requireNamespace("openxlsx", quietly=TRUE)) stop("Install openxlsx before tests")
a <- utils::read.csv(settings$abundance_file, check.names=FALSE, stringsAsFactors=FALSE)
m <- utils::read.csv(settings$metadata_file, check.names=FALSE, stringsAsFactors=FALSE)
c <- utils::read.csv(settings$contrasts_file, check.names=FALSE, stringsAsFactors=FALSE)
m$sample_id <- as.character(m$sample_id)
validate_tables(a, m, c, settings$min_replicates_per_condition)
stopifnot(identical(settings$rng_mode, "independent"))
# First-run determinism and contrast-order independence.
res <- lapply(seq_len(nrow(c)), function(i) analyze_one_contrast(a, m, c[i,], settings))
res_again <- lapply(rev(seq_len(nrow(c))), function(i) analyze_one_contrast(a, m, c[i,], settings))
for (i in seq_len(nrow(c))) {
  j <- match(i, rev(seq_len(nrow(c))))
  stopifnot(isTRUE(all.equal(res[[i]]$matrix, res_again[[j]]$matrix, tolerance=0)))
  stopifnot(isTRUE(all.equal(res[[i]]$table$log2FC, res_again[[j]]$table$log2FC, tolerance=0)))
  stopifnot(all(is.finite(res[[i]]$table$FDR)))
  stopifnot(all(res[[i]]$table$FDR >= 0 & res[[i]]$table$FDR <= 1))
  stopifnot(all(res[[i]]$table$significant == (res[[i]]$table$FDR <= settings$fdr_threshold & abs(res[[i]]$table$log2FC) >= settings$absolute_log2fc_threshold)))
  tmp <- tempfile("proteomics-test-")
  dir.create(tmp)
  folder <- save_outputs(res[[i]], c[i,], tmp, settings)
  stopifnot(file.exists(file.path(folder, "differential_abundance.csv")),
            file.exists(file.path(folder, "analysis_results.xlsx")),
            file.exists(file.path(folder, "volcano.png")))
  unlink(tmp, recursive=TRUE)
}
# Failure test: duplicate accessions must not be silently collapsed.
bad <- a
bad$accession[2] <- bad$accession[1]
stopifnot(inherits(try(validate_tables(bad,m,c,2),silent=TRUE),"try-error"))
message("PASS: synthetic end-to-end analysis, exports, repeatability, order-independence and input validation.")
