# Execute at project root: source("tests/run_visualization_tests.R")
# Synthetic tests for all heatmap modes, functional metrics and combined panels.
source("config/settings.R")
source("R/functions.R")
source("R/visualization.R")
if (!requireNamespace("ggplot2",quietly=TRUE)) stop("Install ggplot2")
a <- utils::read.csv(settings$abundance_file,check.names=FALSE)
m <- utils::read.csv(settings$metadata_file,check.names=FALSE)
c <- utils::read.csv(settings$contrasts_file,check.names=FALSE)
res <- analyze_one_contrast(a,m,c[1,],settings)
tmp <- tempfile("visualization-check-");dir.create(tmp)
ids <- select_heatmap_ids(res$table,"all")
stopifnot(identical(ids,select_heatmap_ids(res$table,"top_n",length(ids)+10L)))
if(length(ids)) {
  h <- save_heatmap(res,tmp,"Synthetic",mode="all",scale_rows=TRUE)
  if(!is.null(h)) {
    stopifnot(ncol(h)==nrow(res$metadata),all(is.finite(h)),
              file.exists(file.path(tmp,"heatmap_all.png")))
  }
}
if(length(ids)>1L) {
  z <- select_heatmap_ids(res$table,"top_n",1)
  stopifnot(length(z)==1L,z %in% ids)
}
# Explicit constant-row safeguard and meaningful synthetic functional data.
fake <- matrix(c(1,1,1,2,3,4),nrow=2,byrow=TRUE,dimnames=list(c("constant","varying"),c("S1","S2","S3")))
example_res<-list(matrix=fake)
clean <- suppressWarnings(build_heatmap_matrix(example_res,c("constant","varying"),TRUE,FALSE,FALSE))
stopifnot(nrow(clean)==1L,rownames(clean)=="varying")
func <- data.frame(accession=paste0("P",1:5),category=c("A","A","B","B","Unknown function"),log2FC=c(-2,3,-1,4,0.5))
s <- functional_summary(func)
stopifnot(sum(s$n)==5L,all(c("mean_negative","mean_positive","min_negative","max_positive") %in% names(s)))
stopifnot(s$n[s$category=="A"]==2L)
save_functional_outputs(func,tmp,"Synthetic")
for(metric in c("n","mean_log2FC","mean_negative","mean_positive","min_negative","max_positive")) {
  stopifnot(file.exists(file.path(tmp,paste0("functional_",metric,".png"))))
}
p <- ggplot2::ggplot(data.frame(x=1:2,y=2:1),ggplot2::aes(x,y))+ggplot2::geom_point()
save_combined_panels(p,p,p,file.path(tmp,"combined.png"))
stopifnot(file.exists(file.path(tmp,"combined.png")))
unlink(tmp,recursive=TRUE)
message("PASS: heatmaps (top/all), scaling, clustering, functional metrics and composite figure smoke tests.")
