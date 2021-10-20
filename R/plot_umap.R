#' Plot UMAP reduction of an NMF model
#'
#' Plots the sample embeddings matrix of an NMF model, grouping based on the colnames.
#'
#' @param nmf_model object of class \code{\link{RcppML::nmf}}.
#' @param ... arguments to \code{\link{uwot::umap}}
#' @export
#' @import ggplot2
#' @importFrom uwot umap
#' @return ggplot2 object of UMAP coordinates with basic styling, you can use other ggplot2 parameters to style the plot further
#'
plot_umap <- function(nmf_model, ...){
  set.seed(123)
  u <- uwot::umap(t(nmf_model$h), ...)
  df <- data.frame("umap1" = u[, 1], "umap2" = u[, 2], "group" = colnames(nmf_model$h))
  df$size <- 1
  for(i in 1:nrow(df))
    if(df$group[i] != "AML sample") df$size[i] = 1.8
  ggplot(df, aes(x = umap1, y = umap2, color = group)) + geom_point(size = df$size) + theme_void()
}
