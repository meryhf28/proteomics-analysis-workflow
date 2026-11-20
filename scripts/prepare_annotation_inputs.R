# Stage 2: create accession lists for external annotation services.
# Run this AFTER run_analysis.R. No API calls are made or simulated.
source('config/settings.R')
root <- settings$output_dir
s <- utils::read.csv(file.path(root,'analysis_summary.csv'),stringsAsFactors=FALSE)
clean_accession <- function(x) sub('^.*::','',as.character(x))
for (k in seq_len(nrow(s))) {
  dir <- s$output[k]
  p <- utils::read.csv(file.path(dir,'significant_proteins.csv'),stringsAsFactors=FALSE)
  ids <- unique(clean_accession(p$accession))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  writeLines(ids,file.path(dir,'IDs_for_external_annotation.txt'))
}
message('Identifier lists exported. External mapping is a MANUAL step (DAVID / eggNOG).')
