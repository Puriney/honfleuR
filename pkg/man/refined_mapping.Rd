% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/map_cell_location.R
\name{refined_mapping}
\alias{refined_mapping}
\title{Quantitative refinement of spatial inferences}
\usage{
refined_mapping(object, genes.use, cells.num, bins = 64)
}
\arguments{
\item{object}{Seurat object}

\item{genes.use}{Genes to use to drive the refinement procedure.}

\item{cells.num}{Number of centroids candiates to be used to estimate mean
and variance parameters of for each cell and each bin. Default: 2n where n is
number of \code{genes.use}.}

\item{bins}{Number of bins. Default: 64.}
}
\value{
Seurat object, where mapping probabilities for each bin are stored
in object@final.prob
}
\description{
Refines the initial mapping with more complex models that allow gene
expression to vary quantitatively across bins (instead of 'on' or 'off'),
and that also considers the covariance structure between genes.
}
\details{
Full details given in spatial mapping manuscript.
}
\note{
Speed-up twin of \code{refined.mapping}
}

