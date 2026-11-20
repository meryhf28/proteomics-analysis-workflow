# Complete command examples

## A. Synthetic demonstration (RStudio Console)

```r
getwd()                         # Confirm you are at repository root
source("scripts/install_dependencies.R")
source("tests/run_tests.R")     # Unit + synthetic end-to-end checks
source("tests/run_annotation_tests.R")  # Synthetic functional annotation checks
source("run_analysis.R")       # Writes synthetic results/
```

## B. Your independent two-condition experiment

1. Copy your own CSV files into a *private* directory outside GitHub.
2. Open `config/settings.R` and change `abundance_file`, `metadata_file`, and `contrasts_file` to their paths (forward `/` separators work in Windows).
3. Edit `output_dir`, and justify each filtering/imputation/significance parameter.
4. Execute `source("run_analysis.R")`; save a copy of `sessionInfo.txt` for audit.

## C. Starting with PEAKS exports

1. Create `data/peaks_manifest.csv` following `data/templates/peaks_manifest.csv`; every sample needs a matching PEAKS file.
2. Inspect `config/peaks_settings.R` and set organism-specific exclusions.
3. Execute `source("scripts/build_from_peaks.R")`.
4. Inspect `data/processed/identification_filter_audit.csv` before proceeding.
5. Point `config/settings.R` at `data/processed/protein_abundance.csv` and `data/processed/sample_metadata.csv`; create `data/processed/contrasts.csv`.
6. Run `source("run_analysis.R")`.

## D. FAQ

**Why no data in `data/raw`?** Unpublished study data were deliberately excluded. The synthetic example is sufficient to learn the interface.

**Does eggNOG run automatically?** No. You must perform the external annotation and curate category names.

**Can I reproduce the author's published figures?** Not from this public package alone; it has a reusable two-group workflow, not an exact released manuscript dataset.

**Does a green GitHub Actions badge prove statistical correctness?** No. It proves the named automated tests pass in that environment. Scientific appropriateness must be assessed separately.
