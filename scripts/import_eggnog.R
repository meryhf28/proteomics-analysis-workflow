# Optional eggNOG-mapper annotation import.
# Run in a terminal: Rscript scripts/import_eggnog.R input_annotations.tsv output_mapping.csv
# Output is an INTERMEDIATE mapping; categories are not inferred from descriptions.
# eggNOG-mapper commonly provides '#query' / 'query' and 'Description' fields.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2L) stop("Usage: Rscript scripts/import_eggnog.R input.tsv output_mapping.csv")
if (!file.exists(args[1])) stop("Input does not exist: ",args[1])
lines <- readLines(args[1], warn=FALSE)
# Find '#query' header, preserving it; ignore other commented lines.
header <- which(grepl("^#?query\\t", lines, ignore.case=TRUE))
if (length(header)) {
  data <- paste(lines[c(header[1], which(!startsWith(lines,"#") & nzchar(lines)))],collapse="\n")
  x <- utils::read.delim(text=data, check.names=FALSE, stringsAsFactors=FALSE, quote="")
} else {
  x <- utils::read.delim(args[1], check.names=FALSE, stringsAsFactors=FALSE, comment.char="#", quote="")
}
key <- intersect(c("#query","query","Query","query_name"),names(x))
if (!length(key)) stop("Could not find query column; inspect eggNOG export format")
desc <- intersect(c("Description","description"),names(x))
accession <- trimws(as.character(x[[key[1]]]))
out <- data.frame(accession=accession,
  annotation_description=if(length(desc)) as.character(x[[desc[1]]]) else rep("",nrow(x)),
  category=rep(NA_character_,nrow(x)),
  annotation_source="eggNOG-mapper", stringsAsFactors=FALSE)
out <- out[!is.na(out$accession) & nzchar(out$accession),,drop=FALSE]
if (anyDuplicated(out$accession)) stop("Duplicate query accessions; resolve before merging")
dir.create(dirname(args[2]), recursive=TRUE, showWarnings=FALSE)
utils::write.csv(out,args[2],row.names=FALSE,na="")
message("Exported ",nrow(out)," accession-level rows. Curate 'category' before making category plots.")
