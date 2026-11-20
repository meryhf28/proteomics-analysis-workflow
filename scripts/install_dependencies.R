# Install CRAN packages and Bioconductor limma, only if missing.
cran <- c("ggplot2", "openxlsx", "readxl")
missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) install.packages(missing, repos="https://cloud.r-project.org")
if (!requireNamespace("limma",quietly=TRUE)) {
  if (!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager",repos="https://cloud.r-project.org")
  BiocManager::install("limma", ask=FALSE, update=FALSE)
}
message("Dependencies are available.")
