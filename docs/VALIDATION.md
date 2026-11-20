# Software Validation and Reproducibility

## Scope

The workflow supports configurable label-free differential protein-abundance analysis, optional PEAKS preprocessing, functional annotation integration and plotting. Validation uses synthetic fixtures and automated checks. These tests assess implementation behavior, not biological validity for every dataset.

## Automated tests

- `tests/run_tests.R`: syntax, inputs, missing-value imputation, differential-abundance analysis, output generation, deterministic repeatability and contrast-order independence.
- `tests/run_annotation_tests.R`: import of example annotation files, curated category merging, retention of missing annotations and functional plots.
- `tests/run_visualization_tests.R`: heatmaps, scaling and clustering options, functional metrics and composite figures.
- `run_analysis.R` and `scripts/run_functional_analysis.R`: synthetic end-to-end runs.

GitHub Actions runs these checks using `.github/workflows/r-ci.yml`. See the repository's Actions tab for the current status; individual executions are evidence for the commit they tested, not all future changes.

## Reproducibility

The software records R session information, supplies synthetic examples and uses independent deterministic seeds for stochastic imputation per contrast. Results depend on software versions, selected parameters, and data quality. Archive analysis settings and input provenance for each study.

## Historical numerical comparison

An earlier research-specific implementation was numerically compared to archived analytical outputs for preprocessing, imputation and moderated statistics. Private experimental datasets are not distributed. This comparison does not independently validate all features or the generalized workflow across organisms and experimental designs.

## Known limitations

- Intended for independent two-condition comparisons; paired, longitudinal, batch-adjusted or multifactor designs are not supported as shipped.
- Low-quantile imputation with Gaussian jitter is a custom method; its missingness assumptions require dataset-specific review.
- Protein identification quality and biological interpretation require independent assessment.
- Functional category assignments require user curation and identifier mapping checks.
- Category counts are descriptive, not statistical enrichment analyses; external annotation/enrichment services are not called automatically.
- Synthetic tests do not establish universal statistical or biological correctness.

## Run locally

From the repository root:

```r
source("scripts/install_dependencies.R")
source("tests/run_tests.R")
source("tests/run_annotation_tests.R")
source("tests/run_visualization_tests.R")
source("run_analysis.R")
source("scripts/run_functional_analysis.R")
```

Before generating local results, optionally run:

```bash
python scripts/check_public_release.py
```

Review generated tables and graphics and confirm suitability for each experimental design before drawing scientific conclusions.
