#
# Imputation wrapper function including various imputations schemes
# ---- Implement Me ----
#

#' PLSR scheme to perform imputation
#'
#' @param x data.frame of predictors(columns): observed gene expression.
#' @param y Vector of response, observed expression of given landmarked gene.
#' @param y.name Chractor. Name of given landmarked gene.
#' @param do.print Logic. Whether print gene name to screen.
#' @param validation Charactor. Scheme for validation. Default: CV.
#' @import pls
#' @export
plsr_preds_expr <- function(x, y, y.name, do.print = FALSE, validation = "CV"){
  if (class(x) != "data.frame") x <- as.data.frame(x)
  y <- as.numeric(y)
  data <- cbind(x, y)
  colnames(data)[ncol(data)] <- y.name
  f <- as.formula(paste0(y.name, " ~ ."))
  plsr.modal <- plsr(formula = f, data = data, validation = validation)
  plsr.rmsep <- RMSEP(plsr.modal, estimate = validation)
  k <- which.min(plsr.rmsep$val) - 1
  k <- ifelse(k == 0, 1, k) ## excluding intercept item
  plrs.fits <- predict(plsr.modal, x, ncomp = k) ## array: Cells x 1 x 1
  plrs.fits <- plrs.fits[, 1, 1]
  if (do.print) print(y.name)
  return(plrs.fits)
}


#' @title Impute expression of landmarked genes based on other genes via LASSO
#' linear modal
#'
#' @param x Matrix of predictors(columns), observed genes expression.
#' @param y Vector of response, observed expression of given landmarked gene.
#' @param s.use Numeric. See \link{predict.lars}'s parameter \code{s}.
#' @param y.name Chractor. Name of the given landmarked gene.
#' @param do.print Logic. Whether print \code{y.name} on screen.
#' @param use.gram Logic. See \link{lars}'s parameter \code{use.Gram}
#' @return Returns the imputed landmarked gene expession via LASSO modal.
#' @import lars
#' @export
lasso_preds_expr <- function(x, y, s.use = 20, y.name = NULL,
                            do.print = FALSE, use.gram = TRUE) {
  if (class(x) != "matrix") x <- as.matrix(x)
  lasso.modal <- lars(x, as.numeric(y), type = "lasso", max.steps = s.use * 2,
                      use.Gram = use.gram)
  #lasso.fits = predict.lars(lasso.modal,x,type = "fit",s = min(s.use,max(lasso.modal$df)))$fit
  lasso.fits <- predict.lars(lasso.modal, x, type = "fit", s = s.use)$fit
  if (do.print) print(y.name)
  return(lasso.fits)
}


#' @title Tilling LASSO linear modal to impute expression of landmark genes
#'
#' @param x Matrix of predictors(columns), observed genes expression.
#' @param y Vector of response, observed expression of given landmarked gene.
#' @param s.use Numeric. See \link{predict.lars}'s parameter \code{s}.
#' @param y.name Chractor. Name of the given landmarked gene.
#' @param do.print Logic. Whether print \code{y.name} on screen.
#' @param use.gram Logic. See \link{lars}'s parameter \code{use.Gram}
#' @param w Window size to set the training data size relative to input
#' \code{x}. Default: 0.8.
#' @return Returns the imputed landmarked gene expession via LASSO modal.
#' @import lars
#' @export
tilling_lasso_preds_expr <- function(x, y, s.use = 20, y.name = NULL,
                             do.print = FALSE, use.gram = TRUE, w = 0.8) {
  if(nrow(x) < 10) stop("Size of sequenced cells is less than 10")
  x <- as.matrix(x)
  y <- as.numeric(y)

  set.seed(1234)
  n <- nrow(x)
  s <- ceiling(n * (1 - w))
  idx <- 1:n
  ridx <- sample(idx, n, replace = FALSE)
#   if (all.equal(n %% s, 0)){
#     ridx.p <- partition(ridx, sep = rep(s, n/s))
#   } else {
#     ridx.p <- partition(ridx, sep = c(rep(s, n %/% s), n %% s))
#   }
  ridx.p <- partition(ridx, s)
  yhat <- rep(NA, nrow(x))
  parts <- seq_len(length(ridx.p))
  yhat[ridx] <- unlist(lapply(parts,
                   function(p) {
                     idx.un <- ridx.p[[p]]
                     idx.tr <- setdiff(ridx, idx.un)
                     x.tr <- x[idx.tr, ]
                     y.tr <- y[idx.tr]
                     model.p <- lars(x.tr, y.tr, type = "lasso",
                                     max.steps = s.use * 2, use.Gram = use.gram)
                     x.un <- x[idx.un, ]
                     o <- predict.lars(model.p, x.un, type = "fit", s = s.use)$fit
                     as.numeric(o)

                   }))
  if(do.print) print(y.name)
  return(yhat)
}


#' @title Impute gene expression and fill in Seurat object.
#'
#' @description Impute the expression values of import genes, e.g. landmark,
#' genes with other genes as predictors via L1-constrained linear models, i.e.
#' lasso (default), or partial least squares regression (PLSR).
#'
#' @param object Seurat object
#' @param genes.use A vector of genes (predictors) that can be used for
#' building the imputation models.
#' @param genes.fit A vector of response genes to be imputed.
#' @param scheme Imputation strategies: "lasso", "plsr".
#' @param s.use Maximum number of steps taken by the algorithm (lower values
#' indicate a greater degree of smoothing). Only used for "lasso" scheme.
#' @param do.print Print progress (output the name of each gene after it has
#' been imputed).
#' @param gram The use.gram argument passed to \link{lars} of lars package. See
#' \link{lars} for parameter details. Only used for "lasso" scheme.
#' @return Returns a Seurat object where the imputed values have been added to
#' \code{object@@data}. If \code{genes.fit} have already existed in object,
#' their values are imputed again via the selected \code{scheme} and refilled;
#' otherwise they are imputed and added to \code{object@@data}.
#' @note This extends the function of original \code{addImputedScore}.
#' @import lars
#' @import pls
#' @export
setGeneric("fill_imputed_expr",
  function(object, genes.use = NULL, genes.fit = NULL, scheme = "lasso",
           s.use = 20, do.print = FALSE, gram = TRUE)
  standardGeneric("fill_imputed_expr"))
setMethod("fill_imputed_expr", "seurat",
  function(object, genes.use = NULL, genes.fit = NULL, scheme = "lasso",
           s.use = 20, do.print = FALSE, gram = TRUE) {
    genes.use <- set.ifnull(genes.use, object@var.genes)
    genes.fit <- set.ifnull(genes.fit, object@var.genes)
    genes.use <- intersect(genes.use, rownames(object@data))
    genes.fit <- intersect(genes.fit, rownames(object@data))

    fitted.expr <- sapply(genes.fit, function(g) {
                           x = t(object@data[setdiff(genes.use, g), ])
                           y = object@data[g, ]
                           if (scheme == "lasso"){
                             lasso_preds_expr(x, y, s.use = s.use, y.name = g,
                                              do.print, gram)
                           } else if (scheme == "plsr"){
                             plsr_preds_expr(x, y, y.name =g, do.print)
                           } else {}
                         })
    fitted.expr <- as.data.frame(t(fitted.expr))
    # print(fitted.expr)

    genes.old <- intersect(genes.fit, rownames(object@imputed))
    genes.new <- setdiff(genes.fit, rownames(object@imputed))
    if (length(genes.old) > 0) {
      object@imputed[genes.old, ] <- fitted.expr[genes.old, ]
    }

    object@imputed <- rbind(object@imputed, fitted.expr[genes.new, ])
    return(object)
  }
)


