empP <- function(x, nullval) {
    return(length(which(nullval > x)) / length(nullval))
}


#' Calculate p values for PCs given the rotated intact matrix and rotated
#' permuted matrix.
#'
#' @param x Subtset or entire \code{object@@pca.x}
#' @param xp \code{fake.pcVals}
#' @return p-value matrix with genes as rows and PCs on columns
calc_PCpval <- function(x, xp){
  g.names <- rownames(x)
  pc.names <- colnames(x)
  x2 <- c(x) ## matrix-to-vector in column-wise
  idx <- c(col(x))
  pcpval <- sapply(1:length(x2),
                   function(i) {
                     token <- abs(x2[i])
                     pc.idx <- idx[i]
                     null.val <- abs(xp[, pc.idx])
                     empP(token, null.val)
                   })
  pcpval <- matrix(pcpval, length(g.names), length(pc.names))
  rownames(pcpval) <- g.names
  colnames(pcpval) <- pc.names
  return(as.data.frame(pcpval))
}


#' Determine statistical significance of PCA scores.
#'
#' Randomly permutes a subset of data, and calculates projected PCA scores for
#' these 'random' genes. Then compares the PCA scores for the 'random' genes
#' with the observed PCA scores to determine statistical signifance. End result
#' is a p-value for each gene's association with each principal component.
#'
#'
#' @param object Seurat object
#' @param num.pc Number of PCs to compute significance for
#' @param num.replicate Number of replicate samplings to perform
#' @param prop.freq Proportion of the data to randomly permute for each
#' replicate
#' @return Returns a Seurat object where object@@jackStraw2.empP represents
#' p-values for each gene in the PCA analysis. If project.pca is subsequently
#' run, object@@jackStraw2.empP.full then represents p-values for all genes.
#' @references Inspired by Chung et al, Bioinformatics (2014)
#' @export
setGeneric("jackStraw2", function(object, num.pc = 30, num.replicate = 100,
  prop.freq = 0.01, do.print = FALSE) standardGeneric("jackStraw2"))
setMethod("jackStraw2", "seurat",
  function(object, num.pc = 30, num.replicate = 100, prop.freq = 0.01,
           do.print = FALSE) {
    num.pc <- min(num.pc, ncol(object@pca.rot))
    pc.genes <- rownames(object@pca.x)
    data.use <- object@scale.data[pc.genes, ]
    S.MIN <- 2.01
    prop.freq <- max(prop.freq, (S.MIN / length(pc.genes)))
    md.x <- as.matrix(object@pca.x)
    # md.rot <- as.matrix(object@pca.rot)

    fake.pcVals.raw <- sapply(1:num.replicate,
                              function(x) jackRandom2(scaled.data  = data.use,
                                                      prop.use     = prop.freq,
                                                      r1.use       = 1,
                                                      r2.use       = num.pc,
                                                      seed.use     = x),
                              simplify = FALSE)
    # fake.pcVals <- sapply(1:num.pc,function(x)as.numeric(unlist(lapply(1:num.replicate,function(y)fake.pcVals.raw[[y]][,x]))))
    fake.pcVals <- do.call("rbind", fake.pcVals.raw)

    # set.seed(1234)
    # fake.pcVals <- replicate(num.replicate,
                             # jackRandom2(scaled.data  = data.use,
                                         # prop.use     = prop.freq,
                                         # r1.use       = 1,
                                         # r2.use       = num.pc),
                             # simplify = FALSE)
    # fake.pcVals <- do.call("rbind", fake.pcVals)

    rownames(fake.pcVals) <- 1:nrow(fake.pcVals)
    object@jackStraw.fakePC <- as.data.frame(fake.pcVals)
    object@jackStraw.empP <- calc_PCpval(md.x[, 1:num.pc], fake.pcVals)
    return(object)
  }
)

## almost identical with original jackRandom
jackRandom2 <- function(scaled.data, prop.use = 0.01, r1.use = 1, r2.use = 5,
                        seed.use) {
  if(!missingArg(seed.use)){
    set.seed(seed.use)
  }
  rand.genes <- sample(rownames(scaled.data), nrow(scaled.data) * prop.use)
  data.mod <- scaled.data
  data.mod[rand.genes, ] <- shuffleMatRow2(scaled.data[rand.genes, ])
  fake.pca <- prcomp(data.mod)
  fake.x <- fake.pca$x
  return(fake.x[rand.genes, r1.use:r2.use])
}

## almost identical as original shuffleMatRow
shuffleMatRow2 <- function(x) {
  x2 <- t(x)
  ind <- order(c(col(x2)), runif(length(x2)))
  x2 <- matrix(x2[ind], nrow=nrow(x), ncol=ncol(x), byrow=TRUE)
  return(x2)
}
