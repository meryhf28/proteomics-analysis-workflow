# Statistical design, transformations and interpretation

## Scope
Independent biological replicate samples, two-condition comparisons within each organism. Input is protein-level intensity quantified using a compatible label-free pipeline. The model does not account for study blocks, subject pairing, batch, dose, time, random effects or nested experiments.

## Procedure
1. **PEAKS filter (optional).** A file per sample, with accession, unique peptide count, protein confidence score, protein area and description. Quality filter by organism and condition according to configurable thresholds. The default two-unique-peptides, score 25, and two-positive-samples thresholds are illustrative configuration defaults, not universal guidelines. Exclusion strings must be tailored to the organism and experimental design. Retain the filter audit.
2. **Assemble matrix.** Unique `accession` keys across organisms, all sample names taken from metadata. Intensities must be nonnegative, otherwise fail. Missing intensities may be blank.
3. **Transform.** For each abundance `x`, compute `log2(x + 1)` by default. The pseudocount is configurable and affects low intensities and zeros.
4. **Median centre.** Subtract the median of observed log-intensities separately in each sample. This approach assumes global changes can be removed by centring; it is not batch correction.
5. **Zero-as-missing.** Original input zeros become `NA` after centring. Zero is assumed to indicate nondetection, not a measured true zero.
6. **Impute.** For each sample, calculate the 0.01 quantile of observed centred values and replace missing with that quantile plus Gaussian jitter. Jitter SD = max(0.01 × sample MAD, 1e-6). This implementation is **not the standard deterministic MinDet method**. The randomness is controlled by a reproducible seed. Independent per-contrast random seeds are used by default to make results insensitive to contrast ordering.
7. **Model.** Fit a separate `limma::lmFit` two-condition, no-intercept design for each organism/contrast; compute treatment minus reference contrast; apply `limma::eBayes` variance moderation. `limma::topTable` gives unadjusted p-values and BH-FDR values for the single coefficient. The workflow does not pool FDR correction across comparisons.
8. **Classify.** Mark proteins significant when both FDR <= 0.05 and absolute log2FC >= 0.58 (editable). Note that the effect-size threshold is applied after inference, not via `limma::treat`.
9. **Display.** Volcano x-axis is log2FC and y-axis is -log10(BH FDR). Use the same FDR threshold for the horizontal guide. The category bar plot counts curated labels, not enrichment or overrepresentation.

## Important cautions
- High rates of missingness and imputing data not missing at random can induce false positives; inspect per-sample intensity and missingness before using the method.
- Proteins missing in all samples of a contrast are excluded, not imputed.
- `n_imputed` counts missing observations for each retained protein within a contrast.
- Labels such as `organism_A` and `treatment` are illustrative only.
- The example dataset is synthetic and cannot validate the assumptions of any biological study.
- Statistical reproducibility of an archived dataset does not in itself demonstrate that the method is the optimal analysis for new data.
