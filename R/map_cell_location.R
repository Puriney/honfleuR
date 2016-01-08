
# Internal, not documented for now
#' @export
setGeneric("map.cell2",
  function(object, cell.name, do.plot = FALSE, safe.use = TRUE,
           text.val = NULL, do.rev = FALSE) standardGeneric("map.cell2")
)
setMethod("map.cell2", "seurat",
  function(object, cell.name, do.plot = FALSE, safe.use = TRUE,
           text.val = NULL, do.rev = FALSE) {
    insitu.matrix <- object@insitu.matrix
    insitu.genes <- colnames(insitu.matrix)
    insitu.genes <- intersect(insitu.genes, rownames(object@imputed))

    insitu.use <- insitu.matrix[, insitu.genes]
    imputed.use <- object@imputed[insitu.genes, ]

    safe_fxn  <- sum
    if (safe.use) {
      safe_fxn <- log_add
    }

    all.avail.states <- lapply(insitu.matrix, unique)
    all.avail.states <- unlist(lapply(insitu.genes,
                                      function(x) paste(x,
                                                        all.avail.states[[x]],
                                                        "post", sep =".")
                                      )
                              )

    missing.states <- setdiff(all.avail.states, colnames(object@mix.probs))
    if (length(missing.states) > 0){
      stop(paste0("Error: ", all.avail.states[missing.states], " is missing from the mixture fits"))
    }

    all.avail.probs <- object@mix.probs[cell.name, all.avail.states] # yellow-page
    bins.avail.probs <- sapply(insitu.genes,
                               function(x) {
                                 g.stat.idx <- paste(x, insitu.matrix[ ,x],
                                                     "post", sep = ".")
                                 log(as.numeric(all.avail.probs[ ,g.stat.idx]))
                               })
    bins.avail.probs <- data.frame(bins.avail.probs)
    tmp <- apply(bins.avail.probs, 2, log_add)
    bins.avail.probs.scale <- bins.avail.probs - matrix(rep(tmp, 64), nrow = 64,
                                                        byrow = TRUE)

    THREHOLD <- -9.2
    bins.avail.probs.scale[bins.avail.probs.scale < THREHOLD] <- THREHOLD

    cell.bins.prob <- exp(apply(bins.avail.probs.scale, 1, safe_fxn))
    cell.bins.prob <- cell.bins.prob / sum(cell.bins.prob)

    #
    # plot
    #
    if (do.plot) {
      par(mfrow = c(1, 2))
      txt.matrix <- matrix(rep("", 64), nrow = 8)
      if (!is.null(text.val)) {
        txt.matrix[text.val] <- "X"
      }
      if (do.rev) {
        bins.avail.probs.scale <- bins.avail.probs.scale[unlist(lapply(0:7, function(x) seq(1, 57, 8) + x)), ]
      }
      aheatmap(matrix(cell.bins.prob, nrow = 8), Rowv = NA, Colv = NA,
               txt = txt.matrix, col = bwCols)
      aheatmap(bins.avail.probs.scale, Rowv = NA, Colv = NA)
      rp()
    }
    ## output
    ## cell to 64 bins probs - numeric vector
    return(cell.bins.prob)

  }
)

#' Infer spatial origins for single cells
#'
#' Probabilistically maps single cells based on (imputed) gene expression
#' estimates, a set of mixture models, and an in situ spatial reference map.
#'
#'
#' @param object Seurat object
#' @param cells.use Which cells to map
#' @return Seurat object, where mapping probabilities for each bin are stored
#' in object@@final.prob
#' @export
setGeneric("initial.mapping",
  function(object,cells.use=NULL) standardGeneric("initial.mapping")
)
#' @export
setMethod("initial.mapping", "seurat",
  function(object, cells.use = NULL) {
    cells.use <- set.ifnull(cells.use, colnames(object@data))
    every.prob <- sapply(cells.use,
                         function(x) map.cell2(object, x, FALSE, FALSE))
    object@final.prob <- data.frame(every.prob)
    rownames(object@final.prob) <- paste("bin.",
                                         rownames(object@final.prob), sep="")
    return(object)
  }
)


slimdmvnorm=function (x, mean = rep(0, p), sigma = diag(p), log = FALSE) 
{
    x <- matrix(x, ncol = length(x))
  p <- ncol(x)
    dec <- tryCatch(chol(sigma), error = function(e) e)
    tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
      rss <- colSums(tmp^2)
      logretval <- -sum(log(diag(dec))) - 0.5 * p * log(2 * pi) - 0.5 * rss
        names(logretval) <- rownames(x)
        logretval
}

#' Quantitative refinement of spatial inferences
#'
#' Refines the initial mapping with more complex models that allow gene
#' expression to vary quantitatively across bins (instead of 'on' or 'off'),
#' and that also considers the covariance structure between genes.
#'
#' Full details given in spatial mapping manuscript.
#'
#' @param object Seurat object
#' @param genes.use Genes to use to drive the refinement procedure.
#' @param cells.num Number of centroids candiates to be used to estimate mean
#' and variance parameters of for each cell and each bin. Default: 2n where n is
#' number of \code{genes.use}.
#' @param bins Number of bins. Default: 64.
#' @return Seurat object, where mapping probabilities for each bin are stored
#' in object@@final.prob
#' @import fpc
#' @export
setGeneric("refined_mapping",
           function(object, genes.use, cells.num, bins = 64)
             standardGeneric("refined_mapping"))
#' @export
setMethod("refined_mapping", "seurat",
  function(object, genes.use, cells.num, bins = 64) {
    genes.use <- intersect(genes.use, rownames(zf@imputed))
    if (missingArg(cells.num)) {
      cells.num <- 2 * length(genes.use)
    }
    cells.name <- colnames(zf@data)
    centroids.pos <- t(sapply(cells.name,
                          function(x) calc_cell_centroid(zf@final.prob[, x])
                          ))
    bins.centroids <- sapply(1:bins, function(b) fetch_closest(b, centroids.pos,
                                                               cells.num))
    ####
    all.mu <- sapply(genes.use,
                     function(gene) sapply(1:64,
                                           function(bin) mean(as.numeric(zf@imputed[gene, fetch.closest(bin, centroids.pos, 2*length(genes.use))] ))))

    all.cov <- list()
    for (x in 1:64) {
      all.cov[[x]] <- cov(t(zf@imputed[genes.use, fetch.closest(x, centroids.pos, 2* length(genes.use))]))
    }

    mv.probs <- sapply(colnames(zf@data),
                       function(my.cell) sapply(1:64 , function(bin) slimdmvnorm(as.numeric(zf@imputed[genes.use, my.cell]), 
                                                                                 as.numeric(all.mu[bin, genes.use]), 
                                                                                 all.cov[[bin]]
                                                                                 )
                       )
                       )
    mv.final <- exp(sweep(mv.probs, 2, apply(my.probs, 2, log_add)))
    object@final.prob <- data.frame(mv.final)
    return(object)



  }
)
