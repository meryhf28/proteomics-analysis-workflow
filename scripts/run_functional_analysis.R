# Run from the repository root with source("scripts/run_functional_analysis.R").
# This module is OPTIONAL: analysis results exist even without annotation.
# Supply curated accession-to-category CSVs for each contrast via the manifest
# specified in config/functional_settings.R. Do not pass raw DAVID enrichment
# p-values as protein-level FDR values.
source("R/visualization.R")
source("config/functional_settings.R")
source("config/settings.R")
if (!file.exists(functional_settings$manifest)) stop("No functional manifest found; see docs/FUNCTIONAL_WORKFLOW.md")
man <- utils::read.csv(functional_settings$manifest,stringsAsFactors=FALSE,check.names=FALSE)
if(!all(c("contrast_dir","mapping_csv") %in% names(man))) stop("Manifest needs contrast_dir,mapping_csv")
if(!nrow(man)) stop("Functional manifest is empty")
plots <- list()
for(i in seq_len(nrow(man))) {
  dir <- file.path(settings$output_dir,man$contrast_dir[i])
  sigfile <- file.path(dir,"significant_proteins.csv")
  mapfile <- man$mapping_csv[i]
  if(!file.exists(sigfile) || !file.exists(mapfile)) stop("Missing input for ",man$contrast_dir[i],": ",sigfile," or ",mapfile)
  a <- utils::read.csv(sigfile,stringsAsFactors=FALSE)
  b <- utils::read.csv(mapfile,stringsAsFactors=FALSE,check.names=FALSE)
  if(!all(c("accession","category") %in% names(b))) stop("Curated table must have accession,category")
  if(anyDuplicated(b$accession)) stop("Duplicate accession in annotation mapping: ",mapfile)
  idx <- match(a$accession,b$accession)
  annotated <- cbind(a,b[idx,setdiff(names(b),"accession"),drop=FALSE])
  if(!"category" %in% names(annotated)) stop("Missing functional category")
  folder <- file.path(dir,"functional")
  dir.create(folder,recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(annotated,file.path(folder,"annotated_significant_proteins.csv"),row.names=FALSE)
  utils::write.csv(data.frame(total=nrow(a),annotated=sum(!is.na(idx)),unmatched=nrow(a)-sum(!is.na(idx))),
                   file.path(folder,"annotation_coverage.csv"),row.names=FALSE)
  summary <- save_functional_outputs(annotated,folder,man$contrast_dir[i])
  if(!is.null(summary) && nrow(summary)) {
    plots[[man$contrast_dir[i]]] <- list(
      volcano=make_volcano(utils::read.csv(file.path(dir,"differential_abundance.csv")),
          list(organism=man$contrast_dir[i],treatment="Treatment",reference="Reference"),settings),
      functional=plot_functional_summary(summary,"mean_log2FC",man$contrast_dir[i]))
  }
  message("Functional coverage for ",man$contrast_dir[i],": ",sum(!is.na(idx)),"/",nrow(a))
}
if(length(plots)>=2L && isTRUE(functional_settings$combined_panels)) {
  keys<-names(plots)[seq_len(2)]
  # A third panel illustrates the first contrast's functional mean fold changes;
  # this is specified in the filename and is NOT a cross-organism comparison.
  out <- file.path(settings$output_dir,"combined_volcano_pair_and_first_functional.png")
  save_combined_panels(plots[[keys[1]]]$volcano,plots[[keys[2]]]$volcano,
                       plots[[keys[1]]]$functional,out)
  message("Combined panels: ",out)
}
message("Functional analysis complete (provided curated categories only).")
