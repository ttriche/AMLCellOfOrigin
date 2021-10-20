#' Cluster samples on an NMF embedding
#'
#' Plots the sample embeddings matrix of an NMF model, grouping based on the colnames.
#'
#' @param nmf_model object of class \code{\link{RcppML::nmf}}.
#' @param ... arguments to \code{\link{Seurat::FindClusters}}
#' @export
#' @import Seurat
#' @importFrom uwot umap
#' @return ggplot2 object of UMAP coordinates with basic styling, you can use other ggplot2 parameters to style the plot further
#'
cluster_nmf <- function(nmf_model, data, resolution = 1){
  rownames(data) <- paste0("feature", 1:nrow(data))
  colnames(data) <- paste0("sample", 1:ncol(data))
  A <- suppressWarnings(CreateSeuratObject(data))
  A@reductions$nmf <- CreateDimReducObject(embeddings = t(nmf_model$h), loadings = nmf_model$w, stdev = nmf_model$d, key = "nmf_", assay = "RNA")
  rownames(A@reductions$nmf@cell.embeddings) <- colnames(A@assays$RNA@counts)
  A <- FindNeighbors(A, k.param = 10, reduction = "nmf", dims = 1:ncol(nmf_model$w))
  A <- FindClusters(A, resolution = resolution, verbose = FALSE)
  return(as.numeric(as.vector(A@meta.data$seurat_clusters)) + 1)
}
