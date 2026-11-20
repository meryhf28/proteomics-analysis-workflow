# Changelog

## [0.9.0-rc4] - 2026-11-20

### Added
- Configurable sample-level heatmaps (Top N, all eligible proteins and custom selections).
- Curated functional category summaries, visualizations and optional composite figures.
- Additional synthetic annotation and plotting tests, plus CI coverage.

### Fixed
- Explicitly load visualization functions in the clean-session test runner.

### Validation
- Synthetic core, annotation, visualization and end-to-end tests completed in R 4.3.3 on Windows.
- GitHub Actions tests passed for the corresponding tested revision; the latest status is available in Actions.

### Statistical scope
- The two-condition limma model, preprocessing, imputation and selection thresholds were unchanged by this release candidate's visualization changes.

## [0.9.0-rc3] - 2026-11-20

- Generalized the workflow for configurable abundance tables, metadata and contrasts.
- Added standardized CSV/XLSX exports, independent random seeds and synthetic fixtures.
- Added documentation, functional annotation import utilities and automated tests.
- Aligned volcano y-axis with the BH-FDR selection threshold.

## Development history

The general-purpose software derives from an earlier research-specific R implementation. Private experimental inputs and archived study outputs are excluded. See [validation](docs/VALIDATION.md) for testing scope.
