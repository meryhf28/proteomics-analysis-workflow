# Proteomics Analysis Workflow in R

**Version 0.9.0-rc4** · **Author: María Hernández-Fernández** · **License: MIT**

A configurable, reproducible R workflow for **label-free differential protein-abundance analysis**, with optional import of PEAKS protein tables, functional annotation integration, and publication-ready data visualizations. The workflow is designed for researchers analyzing independent two-group comparisons across biological datasets. It does not require a specific organism, protein identifier format, or treatment name.

The repository provides **synthetic demonstration data only**; no experimental datasets are distributed.

## Features

- Optional preparation of protein-abundance matrices from PEAKS CSV/XLSX exports using configurable identification and detection filters.
- Sample-wise log2 transformation, median centering, documented low-intensity missing-value imputation, and differential analysis with **limma**.
- Benjamini–Hochberg false discovery rate (FDR) adjustment and configurable significance thresholds.
- Reproducible comparisons with deterministic random-number initialization per contrast.
- Exports in CSV/XLSX format, volcano plots, sample-level heatmaps (top-N and all selected proteins), and analysis metadata.
- Optional accession-list and FASTA preparation for external annotation tools, import of eggNOG-mapper annotations, and integration of manually curated functional categories.
- Functional category summaries and plots, with optional combined panels. **Category summaries are descriptive, not statistical enrichment tests.**

## Requirements

R and the packages `limma` (Bioconductor), `ggplot2`, `openxlsx`, and `readxl`. RStudio is optional.

Install any missing dependencies from the repository root:

```r
source("scripts/install_dependencies.R")
```

## Quick start

Run the included synthetic example from the repository root:

```r
source("scripts/install_dependencies.R")
source("tests/run_tests.R")
source("tests/run_annotation_tests.R")
source("tests/run_visualization_tests.R")
source("run_analysis.R")
source("scripts/run_functional_analysis.R")  # Optional: curated synthetic annotations
```

The `results/` directory is generated locally and contains per-contrast tables and figures, an analysis summary, and `sessionInfo.txt`. Generated outputs are not part of the source-code release.

### Run the workflow on your own data

Prepare three input files:

| File | Required columns | Purpose |
|---|---|---|
| `protein_abundance.csv` | `accession`, `description`, one column per sample | Nonnegative raw protein intensities; missing values may be blank |
| `sample_metadata.csv` | `sample_id`, `organism`, `condition` | Links each sample column to a group and dataset |
| `contrasts.csv` | `organism`, `reference`, `treatment` | Defines independent two-group comparisons |

Sample IDs must match abundance-column names exactly. Positive log2 fold changes indicate higher abundance in the `treatment` condition relative to `reference`.

1. Use the example files in `data/example/` as templates and replace them with your own data.
2. Edit input paths and analytical parameters in `config/settings.R`.
3. Run `source("run_analysis.R")` from the repository root.
4. Inspect the exported tables, diagnostic plots, and `sessionInfo.txt` before interpreting results.

For optional PEAKS preprocessing, see `config/peaks_settings.R`, `data/templates/peaks_manifest.csv`, and `scripts/build_from_peaks.R`. PEAKS export headers, identification filters, and accession conventions must be checked for each dataset.

## Statistical approach and scope

The default workflow applies `log2(intensity + 1)` followed by per-sample median centering. Original zero intensities are treated as missing after transformation. Missing values are replaced by a **low-quantile estimate with MAD-scaled Gaussian jitter**; this is a custom procedure and **not** the canonical MinDet algorithm. It then fits independent, unpaired two-condition models using `limma`, adjusts p-values by the Benjamini–Hochberg method, and reports differential abundance. Default selection thresholds are **FDR ≤ 0.05** and **|log2FC| ≥ 0.58**. Volcano plots display `−log10(FDR)`.

These preprocessing and imputation assumptions are **not universally appropriate**. Researchers should justify them for their experimental design and missingness mechanism, and perform appropriate sensitivity checks. The supplied statistical design does not support paired samples, batch effects, longitudinal measurements, covariate adjustment, or multifactorial designs without modification.

Full details: [Scientific methods and assumptions](docs/SCIENTIFIC_METHODS.md).

## Functional annotation and visualization

Functional analysis is **optional** and requires user-provided annotation tables. The repository does **not** contact eggNOG-mapper or DAVID automatically and does **not** calculate GO/pathway enrichment. Category labels must be curated and identifier matching, annotation coverage, database versions, and any external enrichment background must be independently verified.

See [Annotation workflow](docs/ANNOTATION.md) and [Functional analysis and visualization](docs/FUNCTIONAL_WORKFLOW.md).

## Validation and reproducibility

The repository includes synthetic tests covering input validation, core analysis, repeatability, annotation integration, and visualization. GitHub Actions is configured to run these checks on repository updates; review the [Actions tab](https://github.com/meryhf28/proteomics-analysis-workflow/actions) for the status of the exact revision used.

The workflow has also been compared numerically with previously archived results from an originating analysis, as described in the [validation documentation](docs/VALIDATION.md). This does not establish independent biological validation for every organism, input format, or downstream annotation. Researchers should validate the workflow for their own data and reporting requirements.

## Documentation

- [Spanish user guide](docs/GUIA_ES.md)
- [Execution examples](docs/EXECUTION_EXAMPLES.md)
- [Scientific methods and assumptions](docs/SCIENTIFIC_METHODS.md)
- [Validation and limitations](docs/VALIDATION.md)
- [Functional annotation and plots](docs/FUNCTIONAL_WORKFLOW.md)

## License and citation

Copyright © 2026 María Hernández-Fernández. Distributed under the [MIT License](LICENSE).

Citation metadata is provided in [`CITATION.cff`](CITATION.cff). When a versioned archive and DOI become available, cite the archived software release corresponding to the version used in your analysis.
