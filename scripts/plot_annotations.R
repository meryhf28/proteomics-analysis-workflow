# Stage 5: generic annotations bar plot, optional exploratory heatmap.
# Usage: Rscript scripts/plot_annotations.R annotated.csv output_prefix
# Input from merge_annotations.R contains 'category' and 'log2FC'.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L) stop('Usage: Rscript scripts/plot_annotations.R annotated.csv output_prefix')
if(!requireNamespace('ggplot2',quietly=TRUE)) stop('Install ggplot2')
x <- utils::read.csv(args[1],stringsAsFactors=FALSE)
if(!all(c('category','log2FC') %in% names(x))) stop('Expected category,log2FC')
x <- x[!is.na(x$category) & nzchar(x$category) & is.finite(x$log2FC),,drop=FALSE]
if(!nrow(x)) stop('No valid annotated proteins to plot')
counts <- as.data.frame(table(x$category),stringsAsFactors=FALSE)
names(counts)<-c('category','protein_count')
dir.create(dirname(args[2]), recursive=TRUE, showWarnings=FALSE)
utils::write.csv(counts,paste0(args[2],'_category_counts.csv'),row.names=FALSE)
p <- ggplot2::ggplot(counts,ggplot2::aes(x=stats::reorder(category,protein_count),y=protein_count)) +
  ggplot2::geom_col(fill='#35688C') + ggplot2::coord_flip() +
  ggplot2::labs(x='Functional category',y='Number of differentially abundant proteins',title='Functional annotation overview') +
  ggplot2::theme_bw(base_size=11)
ggplot2::ggsave(paste0(args[2],'_functional_categories.png'),p,width=8,height=max(4,0.3*nrow(counts)+1),dpi=300)
message('Annotation barplot generated. Categories depend on the curated annotation map.')
