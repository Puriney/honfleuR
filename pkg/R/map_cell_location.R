
# Internal, not documented for now

#' @title Helper function of \code{initial_mapping}.
#' @description Map one cell to bins.
#' @export
setGeneric("map_cell",
  function(object, cell.name, gene.exclude = NULL,
           do.plot = FALSE, safe.use = TRUE,
           text.val = NULL, do.rev = FALSE) standardGeneric("map_cell")
)
setMethod("map_cell", "seurat",
  function(object, cell.name, gene.exclude = NULL,
           do.plot = FALSE, safe.use = TRUE,
           text.val = NULL, do.rev = FALSE) {
    insitu.matrix <- object@insitu.matrix
    insitu.genes <- colnames(insitu.matrix)
    insitu.genes <- intersect(insitu.genes, rownames(object@imputed))
    insitu.genes <- setdiff(insitu.genes, gene.exclude)

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
#' @param gene.exclude Gene to be excluded
#' @return Seurat object, where mapping probabilities for each bin are stored
#' in \code{object@@final.prob}
#' @note Speed-up twin of \code{initial.mapping}.
#' @export
setGeneric("initial_mapping",
  function(object, cells.use=NULL, gene.exclude = NULL)
    standardGeneric("initial_mapping")
)
setMethod("initial_mapping", "seurat",
  function(object, cells.use = NULL, gene.exclude = NULL) {
    cells.use <- set.ifnull(cells.use, colnames(object@data))
    # cells.use <- ifelse(is.null(cells.use), colnames(object@data), cells.use)
    every.prob <- sapply(cells.use,
                         function(x) map_cell(object, x,
                                              gene.exclude = gene.exclude,
                                              FALSE, FALSE))
    object@final.prob <- data.frame(every.prob)
    rownames(object@final.prob) <- paste("bin.",
                                         rownames(object@final.prob), sep="")
    return(object)
  }
)

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
#' @importFrom dplyr summarise group_by
#' @importFrom tidyr spread
#' @importFrom mixtools logdmvnorm
#' @note Speed-up twin of \code{refined.mapping}
#' @export
setGeneric("refined_mapping",
           function(object, genes.use, cells.num, bins = 64)
             standardGeneric("refined_mapping"))
#' @export
setMethod("refined_mapping", "seurat",
  function(object, genes.use, cells.num, bins = 64) {
    genes.use <- intersect(genes.use, rownames(object@imputed))
    if (missingArg(cells.num)) {
      cells.num <- 2 * length(genes.use)
    }
    cells.name <- colnames(object@data)
    centroids.pos <- t(sapply(cells.name,
                          function(x) calc_cell_centroid(object@final.prob[, x])
                          ))
    bins.centroids <- t(sapply(1:bins, function(b) fetch_closest(b, centroids.pos,
                                                               cells.num)))
    ## estimate mean of imputed expression value of landmarked genes in all
    ## bins
    permu.centrds <- c(t(bins.centroids)) ## matrix to vector
    cells.num <- ncol(bins.centroids) ## cells.num is 2n+1 due to compatible modal
    permu.bins <- rep(paste0("bin.", seq_len(bins)), each = cells.num)
    gbcenexpr.df <- data.frame(permu.bins, permu.centrds,
                               stringsAsFactors = FALSE)
    temp.rep.idx <- rep(seq_len(length(permu.centrds)), length(genes.use))
    gbcenexpr.df <- gbcenexpr.df[temp.rep.idx, ]
    permu.genes <- rep(genes.use, each = length(permu.centrds))
    gbcenexpr.df <- cbind(permu.genes, gbcenexpr.df, stringsAsFactors = FALSE)
    gbcenexpr.df$expr <- unlist(sapply(seq_len(nrow(gbcenexpr.df)),
                              function(i) {
                                r = gbcenexpr.df$permu.genes[i]
                                c = gbcenexpr.df$permu.centrds[i]
                                as.numeric(object@imputed[r, c])
                              }))
    gb.mu <- summarise(group_by(gbcenexpr.df, permu.genes, permu.bins),
                       mu = mean(expr))
    gb.mu <- as.data.frame(gb.mu)
    gb.mu <- spread(gb.mu, permu.bins, mu)
    rownames(gb.mu) <- gb.mu$permu.genes
    gb.mu <- gb.mu[, -1]
    gb.mu <- gb.mu[genes.use, paste0('bin.', seq_len(bins))]

    gb.cov <- lapply(seq_len(bins), function(b)
                     cov(t(object@imputed[genes.use, bins.centroids[b, ]])))

    ## estimate the density for multivariate normal distribution
    imputed.expr <- t(object@imputed[genes.use, cells.name])
    mvnorm.logden <- t(sapply(seq_len(bins), function(b) {
                            b.mv.mu = gb.mu[, b]
                            b.mv.cov = gb.cov[[b]]
                            logdmvnorm(y = imputed.expr,
                                       mu = b.mv.mu, sigma = b.mv.cov)
                     }))

    ## substract the log_add
    mv.final <- exp(sweep(mvnorm.logden, 2, apply(mvnorm.logden, 2, log_add)))
    object@final.prob <- data.frame(mv.final)
    return(object)
  }
)
