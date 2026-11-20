# Stage 4: join a USER-CURATED accession-to-annotation table.
# Minimal input: accession,category; extra columns are preserved.
# CSVs from DAVID/eggNOG may need manual header and ID mapping.
# This script DOES NOT claim to parse arbitrary proprietary exports.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L) stop('Usage: Rscript scripts/merge_annotations.R significant_proteins.csv annotation_mapping.csv output.csv')
if (!all(file.exists(args[1:2]))) stop('Missing significant proteins or annotation mapping file')
a <- utils::read.csv(args[1],check.names=FALSE,stringsAsFactors=FALSE)
b <- utils::read.csv(args[2],check.names=FALSE,stringsAsFactors=FALSE)
if(!'accession' %in% names(a) || !all(c('accession','category') %in% names(b))) stop('Required mapping columns accession,category')
if(anyDuplicated(b$accession)) stop('One row per accession required; aggregate multi-term mappings explicitly first')
# Match by accession without silently dropping proteins lacking annotations.
idx <- match(a$accession,b$accession)
c <- cbind(a, b[idx,setdiff(names(b),'accession'),drop=FALSE])
dir.create(dirname(args[3]), recursive=TRUE, showWarnings=FALSE)
utils::write.csv(c,args[3],row.names=FALSE)
message('Annotation coverage: ',sum(!is.na(idx)),' / ',nrow(a),' proteins')
