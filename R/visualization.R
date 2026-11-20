# Reusable visualization and export functions for any two-group protein contrast.
# These functions use the differential abundance table AND the sample-level
# imputed abundance matrix; fold changes alone are NOT suitable for sample heatmaps.

select_heatmap_ids <- function(table, mode="top_n", n=50L, custom_ids=NULL) {
  if (!all(c("accession", "FDR", "log2FC", "significant") %in% names(table))) stop("Incomplete statistical table")
  sig <- table[!is.na(table$significant) & table$significant, , drop=FALSE]
  sig <- sig[order(sig$FDR, -abs(sig$log2FC), sig$accession), , drop=FALSE]
  if (mode=="all") return(unique(as.character(sig$accession)))
  if (mode=="top_n") {
    if (length(n)!=1 || !is.finite(n) || n<1) stop("top_n must be a positive integer")
    return(head(unique(as.character(sig$accession)), as.integer(n)))
  }
  if (mode=="custom") {
    if (is.null(custom_ids)) stop("Provide custom_ids for custom heatmaps")
    # Curated selections must be drawn from the significant set; unmatched IDs
    # are explicitly reported to prevent silently plotting irrelevant entries.
    missing <- setdiff(unique(custom_ids), sig$accession)
    if (length(missing)) warning("Custom IDs not significant or absent: ",paste(missing,collapse=", "))
    return(intersect(unique(as.character(custom_ids)),as.character(sig$accession)))
  }
  stop("Unknown heatmap selection mode: ", mode)
}

build_heatmap_matrix <- function(result, ids, scale_rows=TRUE, cluster_rows=TRUE,
                                 cluster_columns=TRUE) {
  x <- result$matrix
  if (!is.matrix(x) || !is.numeric(x) || is.null(rownames(x))) stop("Numeric imputed protein-by-sample matrix required")
  ids <- intersect(ids, rownames(x))
  if (!length(ids)) return(NULL)
  m <- x[ids, ,drop=FALSE]
  # Rows with zero across-sample variance have no defined z-score.
  # Drop them when row scaling is enabled; report what was removed.
  if (scale_rows) {
    sds <- apply(m,1,stats::sd)
    keep <- is.finite(sds) & sds>0
    if (any(!keep)) warning(sum(!keep)," proteins omitted due to zero/undefined variance in row z-scores")
    m <- m[keep,,drop=FALSE]
    if (!nrow(m)) return(NULL)
    means <- rowMeans(m)
    m <- (m-means)/sds[keep]
  }
  if (cluster_rows && nrow(m)>=2L) {
    idx <- stats::hclust(stats::dist(m),method="complete")$order
    m <- m[idx,,drop=FALSE]
  }
  if (cluster_columns && ncol(m)>=2L) {
    idx <- stats::hclust(stats::dist(t(m)),method="complete")$order
    m <- m[,idx,drop=FALSE]
  }
  m
}

save_heatmap <- function(result, dir, label, mode="top_n", n=50L,
                         scale_rows=TRUE, cluster_rows=TRUE,
                         cluster_columns=TRUE, custom_ids=NULL) {
  ids <- select_heatmap_ids(result$table,mode,n,custom_ids)
  mat <- build_heatmap_matrix(result, ids, scale_rows,cluster_rows,cluster_columns)
  path_prefix <- file.path(dir,paste0("heatmap_",mode))
  if (is.null(mat)) {
    writeLines("No eligible nonconstant protein rows; heatmap not generated.",paste0(path_prefix,"_SKIPPED.txt"))
    return(invisible(NULL))
  }
  # ggplot2 tiles do not require pheatmap: clustering order is computed above.
  df <- data.frame(
    accession=rep(rownames(mat),times=ncol(mat)),
    sample_id=rep(colnames(mat),each=nrow(mat)),
    value=as.vector(mat), stringsAsFactors=FALSE)
  df$accession <- factor(df$accession,levels=rev(rownames(mat)))
  df$sample_id <- factor(df$sample_id,levels=colnames(mat))
  p <- ggplot2::ggplot(df,ggplot2::aes(x=sample_id,y=accession,fill=value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low="#3867A6",mid="white",high="#BC494B",midpoint=0,
        name=if(scale_rows) "Row z-score" else "Log2 centered") +
    ggplot2::labs(x="Sample", y="Protein",title=paste(label,"—",mode,"significant proteins")) +
    ggplot2::theme_minimal(base_size=10) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1),panel.grid=ggplot2::element_blank())
  # For large heatmaps, suppress row labels but preserve IDs in CSV.
  if(nrow(mat)>35) p <- p+ggplot2::theme(axis.text.y=ggplot2::element_blank(),axis.ticks.y=ggplot2::element_blank())
  utils::write.csv(cbind(accession=rownames(mat),as.data.frame(mat)),paste0(path_prefix,"_matrix.csv"),row.names=FALSE)
  ggplot2::ggsave(paste0(path_prefix,".png"),p,width=8,height=max(4,min(15,nrow(mat)*.16+2)),dpi=300)
  invisible(mat)
}

make_volcano <- function(df, spec, settings) {
  df$significance_label <- ifelse(df$significant,"Significant","Not significant")
  df$significance_label <- factor(df$significance_label,levels=c("Not significant","Significant"))
  ggplot2::ggplot(df,ggplot2::aes(x=log2FC,y=-log10(pmax(FDR,.Machine$double.xmin)),colour=significance_label)) +
    ggplot2::geom_point(alpha=.75,size=1.6,na.rm=TRUE) +
    ggplot2::geom_vline(xintercept=c(-settings$absolute_log2fc_threshold,settings$absolute_log2fc_threshold),linetype="dashed") +
    ggplot2::geom_hline(yintercept=-log10(settings$fdr_threshold),linetype="dotted") +
    ggplot2::scale_colour_manual(values=c("Not significant"="#EC817A","Significant"="#26B7C0"),drop=FALSE) +
    ggplot2::labs(title=paste(spec$organism,":",spec$treatment,"vs",spec$reference),
                  x="log2 fold change",y="-log10(BH-adjusted p-value)",colour="Protein classification") +
    ggplot2::theme_bw(base_size=11)
}

# Functional analysis is not GO enrichment. Category labels must be supplied
# by the user and are never guessed from free-text descriptions.
functional_summary <- function(annotated, exclude_unknown=FALSE) {
  if (!all(c("accession","category","log2FC") %in% names(annotated))) stop("Expected accession, category and log2FC")
  a <- annotated
  a$category <- trimws(as.character(a$category))
  valid <- !is.na(a$category) & nzchar(a$category) & is.finite(a$log2FC)
  if (exclude_unknown) valid <- valid & !tolower(a$category) %in% c("unknown function","uncharacterized","unknown")
  a <- a[valid,,drop=FALSE]
  if (!nrow(a)) return(data.frame())
  # Every accession must have one curated category for these summaries;
  # duplicated accessions inflate counts and are rejected.
  if(anyDuplicated(a$accession)) stop("One category per accession required for summaries; resolve duplicates explicitly")
  splitcat <- split(a$log2FC,a$category)
  do.call(rbind,lapply(names(splitcat),function(cat){
    v <- splitcat[[cat]]
    data.frame(category=cat,n=length(v),mean_log2FC=mean(v),
       n_negative=sum(v<0),n_positive=sum(v>0),
       mean_negative=if(any(v<0)) mean(v[v<0]) else NA_real_,
       mean_positive=if(any(v>0)) mean(v[v>0]) else NA_real_,
       min_negative=if(any(v<0)) min(v[v<0]) else NA_real_,
       max_positive=if(any(v>0)) max(v[v>0]) else NA_real_,
       stringsAsFactors=FALSE)
  }))
}

plot_functional_summary <- function(summary, metric="n", title="Functional categories") {
  available <- c("n","mean_log2FC","mean_negative","mean_positive","min_negative","max_positive")
  if(!metric %in% available) stop("Unknown functional metric: ",metric)
  if(!nrow(summary)) stop("No annotated proteins for functional plots")
  df <- summary[is.finite(summary[[metric]]),,drop=FALSE]
  if(!nrow(df)) stop("No finite values for this functional summary metric")
  df$category <- factor(df$category, levels=df$category[order(df[[metric]])])
  ggplot2::ggplot(df,ggplot2::aes(x=category,y=.data[[metric]]))+
    ggplot2::geom_col(fill="#35688C")+ggplot2::coord_flip()+
    ggplot2::labs(x="Functional category",y=metric,title=title)+ggplot2::theme_bw(base_size=11)
}

save_functional_outputs <- function(annotated, dir, label="Functional annotation") {
  dir.create(dir,recursive=TRUE,showWarnings=FALSE)
  s <- functional_summary(annotated)
  utils::write.csv(s,file.path(dir,"functional_summary.csv"),row.names=FALSE)
  if(!nrow(s)) {writeLines("No curated functional categories available",file.path(dir,"NO_CURATED_CATEGORIES.txt"));return(invisible(NULL))}
  for(metric in c("n","mean_log2FC","mean_negative","mean_positive","min_negative","max_positive")) {
    # Empty negative/positive facets are expected for some datasets.
    if(!any(is.finite(s[[metric]]))) next
    plot <- plot_functional_summary(s,metric,paste(label,"—",metric))
    ggplot2::ggsave(file.path(dir,paste0("functional_",metric,".png")),plot,
                   width=8,height=max(4,.3*nrow(s)+1),dpi=300)
  }
  invisible(s)
}

# Combined figure: two volcano panels and a category summary, all drawn from
# numerical outputs. No PNG image manipulation or hard-coded organisms.
# R's grid package lays out ggplot grobs without an additional dependency.
save_combined_panels <- function(plot_a,plot_b,plot_c,path) {
  stopifnot(inherits(plot_a,"ggplot"),inherits(plot_b,"ggplot"),inherits(plot_c,"ggplot"))
  grDevices::png(path,width=3200,height=1700,res=250)
  on.exit(grDevices::dev.off(),add=TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,3,widths=grid::unit(c(1,1,1.3),"null"))))
  for (i in seq_len(3)) {
    p<-list(plot_a,plot_b,plot_c)[[i]]
    grid::pushViewport(grid::viewport(layout.pos.row=1,layout.pos.col=i))
    grid::grid.draw(ggplot2::ggplotGrob(p))
    grid::popViewport()
  }
  grid::popViewport()
  invisible(path)
}
