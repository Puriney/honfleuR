#
# Imputation wrapper function including various imputations schemes
# ---- Implement Me ----
#

#' PLSR scheme to perform imputation
#'
#' @param x data.frame of predictors(rows): observed gene expression.
#' @param y Vector of response, observed expression of given landmarked gene.
#' @param y.name Chractor. Name of given landmarked gene. 
#' @param do.print Logic. Whether print gene name to screen.
#' @param validation Charactor. Scheme for validation. Default: CV.
#' @import pls
#' @export
plsr_preds_expr <- function(x, y, y.name, do.print = FALSE, validation = "CV"){
  data <- as.data.frame(t(rbind(x, y)))
  x <- as.data.frame(t(x))
  f <- as.formula(paste0(y.name, " ~ ."))
  plsr.modal <- plsr(formula = f, data = data, validation = validation)
  plsr.rmsep <- RMSEP(plsr.modal, estimate = validation)
  k <- which.min(plsr.rmsep$val)
  plrs.fits <- predict(plsr.modal, x, ncomp = k)
  if (do.print) print(y.name)
  return(plrs.fits)
}
#' @title Impute expression of landmarked genes based on other genes via LASSO
#' linear modal
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
  lasso.modal <- lars(x, as.numeric(y), type = "lasso", max.steps = s.use * 2,
                      use.Gram = use.gram)
  #lasso.fits = predict.lars(lasso.modal,x,type = "fit",s = min(s.use,max(lasso.modal$df)))$fit
  lasso.fits <- predict.lars(lasso.modal, x, type = "fit", s = s.use)$fit
  if (do.print) print(y.name)
  return(lasso.fits)
}

#' @title Calculate imputed expression values
#'
#' @description Uses L1-constrained linear models (LASSO) to impute single cell
#' gene expression values.
#'
#' @param object Seurat object
#' @param genes.use A vector of genes (predictors) that can be used for
#' building the LASSO models.
#' @param genes.fit A vector of genes to impute values for
#' @param scheme Imputation strategies: "lasso", "plsr".
#' @param s.use Maximum number of steps taken by the algorithm (lower values
#' indicate a greater degree of smoothing)
#' @param do.print Print progress (output the name of each gene after it has
#' been imputed).
#' @param gram The use.gram argument passed to lars
#' @return Returns a Seurat object where the imputed values have been added to
#' object@@data
#' @import lars
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
                             plsr_preds_expr(t(x), y, y.name =g, do.print)
                           } else {}
                         })
    fitted.expr <- as.data.frame(t(lasso.fits))

    genes.old <- intersect(genes.fit, rownames(object@imputed))
    genes.new <- setdiff(genes.fit, rownames(object@imputed))
    if (length(genes.old) > 0) {
      object@imputed[genes.old, ] <- fitted.expr[genes.old, ]
    }

    object@imputed <- rbind(object@imputed, fitted.expr[genes.new, ])
    return(object)
  }
)
