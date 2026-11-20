# Stage 3: optional sequence retrieval from an EXISTING local FASTA.
# Usage: Rscript scripts/extract_fasta.R path/to/reference.fasta path/to/IDs.txt path/to/output.fasta
# Assumption: accession occurs as a FASTA token (e.g. >sp|P12345|...)
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=3L) stop('Usage: Rscript scripts/extract_fasta.R input.fasta ids.txt output.fasta')
if (!all(file.exists(args[1:2]))) stop('FASTA or ID list not found')
ids <- trimws(readLines(args[2],warn=FALSE)); ids<-unique(ids[nzchar(ids)])
lines <- readLines(args[1],warn=FALSE)
starts <- which(startsWith(lines,'>'))
if (!length(starts)) stop('FASTA input contains no records')
ends <- c(starts[-1]-1L,length(lines))
extract_id <- function(header) {
  h<-sub('^>','',header)
  tok<-strsplit(h,'[[:space:]]+')[[1]][1]
  parts<-strsplit(tok,'\\|')[[1]]
  if(length(parts)>=3L && tolower(parts[1]) %in% c('sp','tr')) parts[2] else tok
}
ref <- vapply(lines[starts],extract_id,character(1))
if(anyDuplicated(ref)) warning('Duplicate accessions in FASTA; first record is used')
selected<-match(ids,ref); missing<-ids[is.na(selected)]; selected<-unique(selected[!is.na(selected)])
if(length(selected)) {
  rows<-unlist(Map(seq.int,starts[selected],ends[selected]),use.names=FALSE)
  writeLines(lines[rows],args[3])
} else writeLines(character(0),args[3])
writeLines(missing,paste0(args[3],'.unmatched.txt'))
message('Found ',length(selected),' accessions; not found ',length(missing))
