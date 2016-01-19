eval_seurat_innter <- function(object, g){
  ## intinal mapping excluding specific landmark gene
  object <- initial_mapping(zf, gene.exclude = g)
  ## refined mapping excluding same specific landmark gene
  num.pc <- 3
  num.genes <- 6
  genes.use <- pcTopGenes(object, pc.use = 1:num.pc, num.genes = num.genes,
                          use.full = TRUE, do.balanced = TRUE)
  # genes.use     <- setdiff(genes.use, g)
  object <- refined.mapping(zf, genes.use)
  ## generating roc object and AUC value
  probs.use <- object@final.prob
  data.use  <- exp(object@data) - 1
  bins.num  <- nrow(object@insitu.matrix)
  insilico.vector <- unlist(lapply(1:bins.num, function(b) {
                                     prob = as.numeric(probs.use[b, ])
                                     obsv = as.numeric(data.use[g, ])
                                     sum(prob * obsv)
                                    }))
  probs.total <- rowSums(probs.use)
  probs.total[probs.total < 0] <- 0

  pred <- insilico.vector / probs.total
  raw  <- object@insitu.matrix[, g]
  g.roc.obj <- roc(raw, pred)
  g.roc.obj
}

plot_roc_list <- function(l, main = " "){
  nms <- names(l)
  n   <- length(nms)
  mycolor <- brewer.pal(n, 'Set1')
  i <- 1
  auc.vals <- rep(0, n)
  for (obj in l) {
    if (i == 1){
      plot.roc(obj, col = mycolor[i], asp = 1, main = main,
               cex = 1.5, cex.lab = 1.5, cex.axis = 1.5)
    } else {
      lines.roc(obj, col = mycolor[i])
    }
    auc.vals[i] <- as.numeric(obj@auc)
    i <- i + 1
  }
  auc.vals <- format(auc.vals, digits = 3)
  legend.labs <- paste0(nms, " (", auc.vals, ")")
  legend("bottomright", legend = legend.labs, col = mycolor[1:n], lwd = 2)
  0
}
#' @title Evaluate Seurat package performance
#'
#' @description Compare the inferred and observed in-situ patterns of landmark
#' genes to evaluate Seurat performance. Refer to paper 'Evaluating Seurat's
#' perforamce' for details.
#' @param object The seurat object.
#' @param genes.eval The landmark genes to be used to evaluted performance.
#' @import pROC
#' @import RColorBrewer
#' @import vioplot
#' @export
setGeneric("eval_seurat",
  function(object, genes.eval, dir) standardGeneric("eval_seurat"))
setMethod("eval_seurat", "seurat",
  function(object, genes.eval, dir) {
    genes.roc.obj <- lapply(genes.eval, function(x)
                            eval_seurat_innter(object, x))
    names(genes.roc.obj) <- genes.eval
    genes.auc.val <- sapply(genes.roc.obj, function(o) o$auc)
    auc.dec.idx <- order(genes.auc.val, na.last = TRUE, decreasing = TRUE)
    genes.roc.obj <- genes.roc.obj[auc.dec.idx]

    genes.eval.default <- c("ADMP", "OSR1", "CDX4", "SOX3", "CHD", "SZL")

    genes.roc.default <- genes.roc.obj[genes.eval %in% genes.eval.default]
    pdf(paste0(dir, '/eval_seurat1', 7, 7))
    plot_roc_list(genes.roc.default)
    dev.off()

    if (length(genes.eval) > 6){
      pdf(paste0(dir, '/eval_seurat2', 7, 7))
      plot_roc_list(genes.roc.obj[1:6], main = "Top 6 genes")
      dev.off()
    } else {
      pdf(paste0(dir, '/eval_seurat2', 7, 7))
      plot_roc_list(genes.roc.obj, main = "ROC for all available landmark genes")
      dev.off()
    }

    pdf(paste0(dir, '/eval_seurat3', 3, 7))
    vioplot(as.numeric(genes.auc.val), names = "", main = "Prediction ROC",
            col = "tomato", rectCol = "grey")
    abline(h = median(as.numeric(genes.auc.val), lty = 2))
    text(x = 1, y = median(as.numeric(genes.auc.val)),
         median(as.numeric(genes.auc.val)))
    dev.off()

    return(as.data.frame(genes.auc.val))
  }
)

