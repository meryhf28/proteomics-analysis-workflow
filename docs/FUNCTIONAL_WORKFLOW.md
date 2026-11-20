# Complete visualization and functional annotation workflow

The downstream annotation and visualization modules support configurable group labels and identifiers without assuming a particular organism or treatment.

## Step A — statistical analysis and heatmaps

Run `source("run_analysis.R")` in RStudio from the repository root. The configured two-group contrasts each produce differential abundance results, FDR-based volcano plots, and two heatmaps:

- `heatmap_top_n.png`: the first N **significant** proteins ranked by FDR (ties by absolute log2FC), default N=50; if fewer than N are significant, plot those available.
- `heatmap_all.png`: all significant proteins. With fewer than one nonconstant protein, a `SKIPPED.txt` note is emitted.

Both use the sample-level imputed and median-centred log2 abundance matrix; each row is by default standardized to a z-score. **A row z-score depicts relative sample patterns, not an absolute protein abundance or log2 fold change.** Constant rows are excluded from row scaling. Hierarchical clustering uses Euclidean distance and complete linkage, independently on proteins and samples. Matrices underlying plots are exported as CSV. Heatmaps are generated with ggplot2; they are not fabricated from contrast-level fold changes.

Users can modify `heatmap_enabled`, `heatmap_top_n`, `heatmap_row_scaling`, `heatmap_cluster_rows`, and `heatmap_cluster_columns` in `config/settings.R`. The exported matrix retains accession/sample identifiers even when a large heatmap hides protein labels.

## Step B — external annotations (optional)

An external eggNOG result can be prepared using `scripts/import_eggnog.R`; a FASTA can be generated locally via `scripts/extract_fasta.R`. The workflow **does not submit data to external services**. DAVID tables and enrichment terms are NOT interchangeable with protein-level differential abundance significance; document organism-specific accession mapping, background, database versions and enrichment correction separately.

For functional barplots, create a curated mapping with columns `accession,category` and optional `block,annotation_source,annotation_description`. Each accession must be unique. Categories may be a published ontology, a manual classification or an agreed functional grouping, but must not be misrepresented as statistically enriched pathways. Store a manifest with `contrast_dir,mapping_csv`, one row per contrast, and set its path in `config/functional_settings.R`.

Then run `source("scripts/run_functional_analysis.R")` from the repository root, after statistical analysis.

Exports for each contrast include: annotated proteins, annotation coverage, a quantitative category summary, and **six** graph types: per-category counts, mean log2FC, negative-only mean, positive-only mean, negative minimum, and positive maximum. Zero observations are not silently substituted for absent positive or negative categories. Functional category counts represent proteins, not enrichment significance.

When two or more contrasts are provided, an optional three-panel image combines the first two volcano plots and the **first contrast's** functional mean-log2FC plot. The filename makes this explicit. The combined layout is for presentation only and does not imply direct statistical comparison of functional categories across contrasts.

## Step C — special-purpose tables

- `differential_abundance.csv` / `analysis_results.xlsx`: all tested proteins, sorted by FDR.
- `significant_proteins.csv`: proteins meeting BOTH FDR and absolute log2FC thresholds.
- `imputed_matrix.csv`: sample-level values for the per-contrast heatmaps.
- `annotated_significant_proteins.csv`: optional one-row-per-accession curated annotations.

A selection of 43 proteins from a specific manuscript was intentionally **not** hard-coded; researchers can filter any exported table according to clearly documented scientific criteria. Custom heatmap accession selections are supported through the `save_heatmap(..., mode="custom", custom_ids=...)` function.

## Important limitations

- A single, independent two-group limma model per contrast is supported; more elaborate paired, batch, repeated-measures or multifactor designs require a different statistical design matrix.
- Imputation uses an explicitly documented low-quantile Gaussian-jitter strategy **not** an established `imputeLCMD::impute.MinDet` call.
- DAVID gene-list and chart workbooks cannot generally be merged safely without accession-to-gene disambiguation. Input preparation is required.
- Generalized plotting routines have synthetic tests; independent validation is still required for user-supplied datasets.
