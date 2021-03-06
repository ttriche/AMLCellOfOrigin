% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_nmf.R
\name{cluster_nmf}
\alias{cluster_nmf}
\title{Cluster samples on an NMF embedding}
\usage{
cluster_nmf(nmf_model, data, resolution = 1)
}
\arguments{
\item{nmf_model}{object of class \code{\link{RcppML::nmf}}.}

\item{...}{arguments to \code{\link{Seurat::FindClusters}}}
}
\value{
ggplot2 object of UMAP coordinates with basic styling, you can use other ggplot2 parameters to style the plot further
}
\description{
Plots the sample embeddings matrix of an NMF model, grouping based on the colnames.
}
