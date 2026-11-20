# ================================================================
# OPTIONAL STAGE 0 — PEAKS protein identification -> abundance CSV
# ================================================================
# Protein identification: use the maximum unique
# peptides >= 2, maximum score >= 25, and intensity >0 in >=2
# replicates, evaluated WITHIN EACH ORGANISM across its two groups.
# Only organisms matching those manifest rows enter a contrast.
# The user must justify these thresholds in config/peaks_settings.R.
# Input manifest: path to one PEAKS export per biological sample.
# Accepted files: .xlsx or .csv, with explicitly mapped column names.
# No sample groups or species are hardcoded.
source('config/peaks_settings.R')
if (!requireNamespace('readxl',quietly=TRUE)) stop('Install readxl to parse PEAKS Excel exports.')
if (!file.exists(peaks$manifest_file)) stop('PEAKS manifest missing: ', peaks$manifest_file)
manifest <- utils::read.csv(peaks$manifest_file, stringsAsFactors=FALSE, check.names=FALSE)
# Keep numeric-looking sample IDs as column names, not positional indices.
manifest$sample_id <- as.character(manifest$sample_id)
needed <- c('sample_id','organism','condition','file')
if (!all(needed %in% names(manifest))) stop('Manifest must have: ', paste(needed,collapse=', '))
if (nrow(manifest)==0L || anyNA(manifest[,needed]) || anyDuplicated(manifest$sample_id)) stop('Manifest incomplete or duplicates')
if (any(!nzchar(trimws(manifest$sample_id))) || any(!nzchar(trimws(manifest$file)))) stop('Manifest has empty sample IDs or paths')
if (!all(file.exists(manifest$file))) stop('Missing PEAKS file(s): ',paste(manifest$file[!file.exists(manifest$file)],collapse=', '))
get_col <- function(nms, choices, mandatory=TRUE) {
  found <- choices[choices %in% nms]
  if (!length(found)) {if (mandatory) stop('No matching column found among: ', paste(choices,collapse=' / ')); return(NA_character_)}
  found[[1]]
}
clean_id <- function(x) {
  x <- trimws(as.character(x)); x <- sub('^(sp|tr)\\|','',x,ignore.case=TRUE)
  x <- sub('\\|.*$','',x); x <- sub(' .*$','',x); x <- toupper(x)
  x[!nzchar(x)] <- NA_character_; x
}
load_sample <- function(i) {
  f <- manifest$file[i]; sid <- manifest$sample_id[i]
  df <- if (grepl('\\.csv$',f,ignore.case=TRUE)) utils::read.csv(f,check.names=FALSE) else readxl::read_excel(f)
  nms <- names(df)
  ac <- get_col(nms,c('Accession','Protein','ID','Master Protein Accessions'))
  un <- get_col(nms,c('#Unique','Unique','Unique Peptides','Peptides Unique'))
  sc <- get_col(nms,c('-10lgP','minus10lgP','10lgP','Protein Score','Confidence','Score'))
  ar <- get_col(nms,c(paste0('Area Sample ',sid),'Area','Intensity','Abundance'))
  de <- get_col(nms,c('Description','Protein Description','Protein Desc'),FALSE)
  x <- data.frame(accession=clean_id(df[[ac]]),description=if(is.na(de)) rep('',nrow(df)) else as.character(df[[de]]),
                  unique_n=suppressWarnings(as.numeric(df[[un]])), score=suppressWarnings(as.numeric(df[[sc]])),
                  area=suppressWarnings(as.numeric(df[[ar]])), sample_id=sid,
                  organism=manifest$organism[i], condition=manifest$condition[i], stringsAsFactors=FALSE)
  x <- x[!is.na(x$accession),,drop=FALSE]
  if (anyDuplicated(x$accession)) {
    stop('Duplicate protein accession in ',f,'. Resolve ambiguity in PEAKS before using this workflow.')
  }
  if (any(x$area < 0,na.rm=TRUE)) stop('Negative protein intensity in ',f)
  x
}
long <- do.call(rbind,lapply(seq_len(nrow(manifest)),load_sample))
if (length(peaks$exclude_description_regex) && nzchar(peaks$exclude_description_regex)) {
  long <- long[is.na(long$description) | !grepl(peaks$exclude_description_regex,long$description,ignore.case=TRUE,perl=TRUE),,drop=FALSE]
}
if (!nrow(long)) stop('No protein rows left after exclusion filter')
# For each organism, maxima are computed WITHIN each condition.
# A protein is retained if it reaches the unique+score criteria in
# at least one group AND the replicate-presence criterion in at
# least one group; these need not be the same group (original logic).
key <- unique(long[,c('organism','accession'),drop=FALSE])
keep <- rep(FALSE,nrow(key))
for (k in seq_len(nrow(key))) {
  x <- long[long$organism==key$organism[k] & long$accession==key$accession[k],,drop=FALSE]
  cond <- unique(x$condition)
  id_ok <- any(vapply(cond,function(g) {
    y<-x[x$condition==g,,drop=FALSE]
    if (!any(is.finite(y$unique_n)) || !any(is.finite(y$score))) return(FALSE)
    max(y$unique_n,na.rm=TRUE)>=peaks$min_unique && max(y$score,na.rm=TRUE)>=peaks$min_score
  },logical(1)))
  presence_ok <- any(vapply(cond,function(g) sum(x$condition==g & !is.na(x$area) & x$area>0)>=peaks$min_present_replicates,logical(1)))
  keep[k] <- id_ok && presence_ok
}
key$keep <- keep
filtered <- merge(long,key,by=c('organism','accession'),all.x=TRUE,sort=FALSE)
filtered <- filtered[filtered$keep,,drop=FALSE]
if(!nrow(filtered)) stop('No proteins passed quality filters. Inspect settings and raw PEAKS data.')
# A protein ID may recur in distinct organisms. For use by the generic
# downstream pipeline, we prefix IDs with organism to remain unique.
filtered$global_id <- paste(filtered$organism,filtered$accession,sep='::')
ids <- unique(filtered$global_id)
samples <- manifest$sample_id
wide <- data.frame(accession=ids,description=rep('',length(ids)),stringsAsFactors=FALSE)
for (j in seq_along(ids)) {
  x <- filtered[filtered$global_id==ids[j],,drop=FALSE]
  desc <- x$description[!is.na(x$description) & nzchar(x$description)]
  if (length(desc)) wide$description[j] <- desc[[1]]
}
for (sid in samples) {
  values <- filtered[filtered$sample_id==sid,c('global_id','area'),drop=FALSE]
  wide[[sid]] <- values$area[match(ids,values$global_id)]
  # Other-organism proteins stay NA, not 0.
}
dir.create(peaks$output_dir,recursive=TRUE,showWarnings=FALSE)
utils::write.csv(wide,file.path(peaks$output_dir,'protein_abundance.csv'),row.names=FALSE,na='')
utils::write.csv(manifest[,c('sample_id','organism','condition')],file.path(peaks$output_dir,'sample_metadata.csv'),row.names=FALSE)
utils::write.csv(key,file.path(peaks$output_dir,'identification_filter_audit.csv'),row.names=FALSE)
message('PEAKS filtering completed. Configure config/settings.R to use the generated files.')
