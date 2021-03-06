#' @title Build mixture models of gene expression
#'
#' @description Models the imputed gene expression values as a mixture of
#' gaussian distributions. For a two-state model, estimates the probability that
#' a given cell is in the 'on' or 'off' state for any gene. Followed by a greedy
#' k-means step where cells are allowed to flip states based on the overall
#' structure of the data (see Manuscript for details)
#'
#' @param object Seurat object
#' @param gene Gene to fit
#' @param do.k Number of modes for the mixture model (default is 2).  When
#' \code{do.k} is set greater than 2, it cannot co-exist with parameter
#' \code{start.pct}.
#' @param num.iter Number of 'greedy k-means' iterations (default is 1)
#' @param do.plot Logic. Plot mixture model results if TRUE.
#' @param genes.use Genes to use in the greedy k-means step. (Default: landmark
#' genes. See manuscript for details)
#' @param start.pct (Optional) Initial estimates of the percentage of cells in
#' the 'on' state (usually estimated from the in situ map). Cannot co-exist with
#' \code{do.k}.  When \code{start.pct} is used, it is assumed that only 2 states
#' exisits and in the same time \code{do.k} is forced to be 2.
#' @return A Seurat object, where the posterior of each cell being in the 'on'
#' or 'off' state for each gene is stored in \code{object@@mix.probs}
#' @note This is the speed-up and debugged twin of \code{fit.gene.k} of
#' \link{seurat}.
#' @export
### @importFrom mixtools normalmixEM [not used at all]
setGeneric("fit_gene_k",
  function(object, gene, do.k = 2, num.iter = 1, do.plot = FALSE,
           genes.use = NULL, start.pct = NULL) standardGeneric("fit_gene_k"))
setMethod("fit_gene_k", "seurat",
  function(object, gene, do.k = 2, num.iter = 1, do.plot = FALSE,
           genes.use = NULL, start.pct = NULL) {
    data <- object@imputed
    data.use <- data[gene, ]
    names(data.use) <- colnames(data.use)
    scale.data <- t(scale(t(object@imputed)))
    genes.use <- set.ifnull(genes.use, rownames(scale.data))
    genes.use <- intersect(genes.use, rownames(scale.data))
    scale.data <- scale.data[genes.use, ]

    data.cut <- as.numeric(data.use[gene, ])
    cell.ident <- as.numeric(cut(data.cut, do.k))
    if (!(is.null(start.pct))) {
      if(!(is.null(do.k)) && do.k > 2) {
        stop("Parameter do.k with greater than 2 cannot coexist with parameter start.pct")
      }
      cell.ident <- rep(1, length(data.cut))
      cell.ident[data.cut > quantile(data.cut, 1-start.pct)] <- 2
    }
    cell.ident <- order(tapply(as.numeric(data.use),cell.ident,mean))[cell.ident]
    ident.table <- table(cell.ident)
    if (num.iter > 0) {
      for(i2 in 1:num.iter) {
        cell.ident <- iter.k.fit.fast(scale.data, cell.ident, data.use, do.k)
        ident.table <- table(cell.ident)
      }
    }
    data.use.t <- as.data.frame(t(data.use))
    data.use.t <- cbind(data.use.t, cell_ident = cell.ident)
    data.use.t.melt <- melt(data.use.t, id.vars = "cell_ident")
    kmodal.mu <- dcast(data.use.t.melt, cell_ident ~ variable, mean)
    kmodal.sd <- dcast(data.use.t.melt, cell_ident ~ variable, sd)
    kmodal.norm_factor <- as.numeric(ident.table / sum(ident.table))
    raw.probs <- sapply(1:do.k, function(k) {
      factor.k <- as.numeric(kmodal.norm_factor[k])
      mean.k <- kmodal.mu[k, gene]
      sd.k <- kmodal.sd[k, gene]
      prob.k <- (factor.k * dnorm(data.use.t[, gene], mean = mean.k, sd = sd.k))
      return(prob.k)
    })
    norm.probs <- as.data.frame(raw.probs / rowSums(raw.probs))
    #colnames(norm.probs)=unlist(lapply(1:do.k,function(x)paste(gene,x-1,"post",sep=".")))
    colnames(norm.probs) <- paste(gene, (1:do.k)-1, "post", sep = ".")
    row.names(norm.probs) <- row.names(data.use.t)
    norm.probs <- cbind(norm.probs, cell.ident)
    colnames(norm.probs)[ncol(norm.probs)] <- paste0(gene, ".ident")

    new.mix.probs <- data.frame(minusc(object@mix.probs,
                                       paste(gene, ".",sep = "")),
                                row.names = rownames(object@mix.probs))
    colnames(new.mix.probs)[1] <- "nGene"
    object@mix.probs <- cbind(new.mix.probs,norm.probs)

    if (do.plot) {
      nCol=2
      num.row=floor((do.k+1)/nCol-1e-5)+1
      hist(as.numeric(data.use),probability = TRUE,ylim=c(0,1),xlab=gene,main=gene);
      for(k in 1:do.k) {
        # hist(as.numeric(data.use),probability = TRUE,ylim=c(0,1),xlab=gene,main=gene);
        factor.k <- as.numeric(kmodal.norm_factor[k])
        mean.k <- kmodal.mu[k, gene]
        sd.k <- kmodal.sd[k, gene]
        lines(seq(-10,10,0.01),
              factor.k * dnorm(seq(-10,10,0.01), mean.k, sd.k),
              col = k + 1, lwd = 2);
      }
    }
    return(object)
  }
)

iter.k.fit.fast <- function(scale.data, cell.ident, data.use, do.k) {
  # cell.ident.K <- sort(unique(cell.ident))
  cell.ident.K <- 1:(do.k)
  means.all <- sapply(cell.ident.K, function(x)
                      rowMeans(scale.data[, cell.ident == x]))
  all.dist <- lapply(cell.ident.K, function(i)
                     calc_dist(scale.data, means.all[, i]))
  all.dist <- matrix(unlist(all.dist), ncol = length(cell.ident.K))
  cell.ident <- apply(all.dist, 1, which.min)
  cell.ident <- order(tapply(as.numeric(data.use), cell.ident, mean))[cell.ident]
  return(cell.ident)
}
