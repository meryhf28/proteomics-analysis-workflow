# Functional annotation: explicit, optional workflow

**There is no universal direct crosswalk from PEAKS accessions to DAVID/eggNOG terms.** Species-specific database accessions and annotation versions matter, and external services change.

## 1. Prepare protein identifiers

After `source("run_analysis.R")`, execute `source("scripts/prepare_annotation_inputs.R")` from the project root. Each contrast's results subfolder gains `IDs_for_external_annotation.txt`.

## 2. If necessary, export a local FASTA

From an RStudio Terminal (not the R Console), run:

```bash
Rscript scripts/extract_fasta.R path/to/reference.fasta path/to/IDs_for_external_annotation.txt path/to/output.fasta
```

Reference FASTA sequence and identifiers must match your protein database. Check the `.unmatched.txt` report. External services (DAVID/eggNOG) are used outside this repository.

## 3. Import an eggNOG result (optional)

For a tab-delimited eggNOG-mapper result with `#query` and `Description` columns:

```bash
Rscript scripts/import_eggnog.R path/to/emapper.annotations path/to/eggnog_mapping.csv
```

The intermediate mapping has `accession`, `annotation_description`, `category`, `annotation_source`. `category` remains blank: **automatic annotation is not equivalent to manual grouping into biological categories**. Inspect/query identifiers; in some exports, IDs are protein FASTA headers not accessions.

## 4. Curate functional categories and merge

Create `annotation_mapping.csv` with `accession,category` (one row per accession) and optional `description`/`source` columns. Explicitly record the category assignment criteria and annotation version. Then run:

```bash
Rscript scripts/merge_annotations.R path/to/significant_proteins.csv path/to/annotation_mapping.csv path/to/annotated.csv
Rscript scripts/plot_annotations.R path/to/annotated.csv path/to/annotation_summary
```

The plot is a **count of categorized significantly abundant proteins**, not a GO enrichment analysis. Labels not annotated are retained as missing and excluded from category counts; inspect coverage.

## 5. DAVID caution

DAVID gene lists and enrichment charts should be processed using explicit accession-to-gene mapping; gene symbols can be duplicated, multi-mapped, or absent. DAVID term p-values (and multiple-comparison adjustments) refer to enrichment, not protein-level limma p-values. The user must select a defensible background universe. This release does not parse arbitrary DAVID workbook exports or reconstruct enrichment calculations.
