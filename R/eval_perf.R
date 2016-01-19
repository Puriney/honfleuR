
#' @title Evaluate Seurat package performance
#'
#' @description Compare the inferred and observed in-situ patterns of landmark
#' genes to evaluate Seurat performance. Refer to paper 'Evaluating Seurat's
#' perforamce' for details.
#' @param object The seurat object.
#' @param genes.eval The landmark genes to be used to evaluted performance.
#' @export
setGeneric("eval_seurat",
  function(object, genes.eval) standardGeneric("eval_seurat"))
setMethod("eval_seurat", "seurat",
  function(object, genes.eval) {

    ## intinal mapping excluding specific landmark gene


    ## refined mapping excluding same specific landmark gene

    ## generating roc object and AUC value

    ## generating boxplot and violin plot



  }
)

