



#
# Map cell to proper location
#

#' @export
fetch.closest=function(bin,all.centroids,num.cell) {
  bin.y=(bin-1)%/%8+1
  bin.x=(bin-1)%%8+1
  all.centroids=rbind(all.centroids,c(bin.x,bin.y))
  all.dist=as.matrix(dist(all.centroids))
  return(names(sort(all.dist[nrow(all.dist),]))[2:(num.cell+2)])
}

#' @title Calculate spatial centroid of given cell. 
#'
#' @description
#' For given one cell, calculate spatial centroid, i.e. center of mass, of the 
#' spatial probability map. [helper function]
#'
#' @param probs A vector of probability of cell originating from all bins. Length
#' is equal to \code{bins}.
#' @param bins A scalar. Number of bins (default: 64). It should be the product
#' of \code{xbins} and \code{ybins}.
#' @export
calc_cell_centroid <- function(probs, bins = 64) {
  xbins <- ybins <- sqrt(bins)
  xidx <- rep(1:xbins, times = ybins)
  yidx <- rep(1:ybins, each = xbins)
  x <- sum(xidx * probs)
  y <- sum(yidx * probs)
  return(c(x, y))
}
