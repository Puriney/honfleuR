
set.ifnull <- function(x,y) {
  if (is.null(x)) {
    x <- y
  }
  return(x)
}

#' @title Distance Calculation
#'
#' @description Compute distance between query vector and all rows (or columns)
#' of data matrix in 1 v.s. N fashion. Alternatively pairwise compute distance
#' among rows (or columns).
#'
#' @param x A matrix or data.frame.
#' @param v A query vector.
#' @param struct A character indicating structure of the \code{x}. String
#' "\code{fn}" means feature by sample; "\code{nf}" means sample by feature.
#' Default: '\code{fn}'.
#' @param method Distance metric (Default: euclidean).
#' @return A vector of distance or pairwise distance matrix.
#' @note Default structure of matrix is feature by sample, which is opposite to
#' conventions. But this is a good compromise so that friendly to vectorized
#' computation and smaller memory usage.
#' @export
calc_dist <- function(x, v, struct = 'fn', method = c('euclidean')) {
  ## x: matrix. F x N
  ## v: vector. 1 x F
  if (struct == 'nf'){
    x <- t(x)
  }
  if (nrow(x) != length(v)) {
    stop("Dimension imcompatible")
  }
  if (method == 'euclidean') {
      return(sqrt(colSums((x - v) ^ 2)))
  }
}


#
# Map cell to proper location
#

#' @title For given bin, fetch cloeset cell-centroids.
#'
#' @description For given bin, find closest cell-centroids in the space.
#' Distance metric is L2-distance.
#'
#' @param b A scalar indicating the bin index.
#' @param centroids A data.frame recording all the cell-centroids candidates
#' with cell names as rows and coordinates as columns.
#' @param cells.num A scalar indicating the number of cell-centroids to be
#' fetched.
#' @param bins A scalar indicating total bins number. Default: 64.
#' @param compatible A logic indicating whether running function compatible with
#' original \code{seura} package function. If \code{TRUE}, this function will
#' actually yieid \code{cells.num + 1} closest centroids. If \code{FALSE},
#' exactly \code{cells.num} centroids will be reported. Default: TRUE.
#'
#' @return A character vector of cells names.
#' @note Please notice the original \code{fetch.closest} function in
#' \code{seurat} package actually
#' report \code{2n+1} cells where \code{n} is the number of landmarked genes,
#' while in paper it claimed \code{2n}.
#' @export
fetch_closest <- function(b, centroids, cells.num, bins = 64,
                          compatible = TRUE) {
  bins.base <- sqrt(bins)
  b.x <- (b - 1) %%  bins.base + 1
  b.y <- (b - 1) %/% bins.base + 1
  b.pos <- c(b.x, b.y)
  b.dists <- calc_dist(x = centroids, v = b.pos, struct = 'nf')
  b.dists.sort <- sort(b.dists)
  if (compatible) cells.num <- cells.num + 1
  out <- names(b.dists.sort)[1:cells.num]
  return(out)
}

#' @title Calculate spatial centroid of given cell.
#'
#' @description
#' For given one cell, calculate spatial centroid, i.e. center of mass, of the
#' spatial probability map. [helper function]
#'
#' @param probs A vector of probability of cell originating from all bins. Length
#' is equal to \code{bins}.
#' @param bins A scalar. Number of bins (default: 64).
#' @export
calc_cell_centroid <- function(probs, bins = 64) {
  xbins <- ybins <- sqrt(bins)
  xidx <- rep(1:xbins, times = ybins)
  yidx <- rep(1:ybins, each = xbins)
  x <- sum(xidx * probs)
  y <- sum(yidx * probs)
  return(c(x, y))
}


partition <- function(x, pSize){
  n <- length(x)
  pNum <- n %/% pSize
  g <- rep(1:pNum, each = pSize)
  out <- n %% pSize
  if (out != 0){
    if(out < pSize / 2){
      g <- c(g, rep(pNum, out))
    } else{
      g <- c(g, rep(pNum + 1, out))
    }
  }
  l <- split(x = x, f = g)
  return(l)
}
