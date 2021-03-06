#' Run seurat analysis pipeline (optional: excluding specific landmark gene)
#'
#' @param object Object of seurat. Load the object before imputation.
#' @param g Specific landmark gene (upper case). If not setting, all landmark genes from in
#' situ data are used. Otherwise, \code{g} is excluded for generating spatial
#' reference map
#' @param scheme Imputation strategy. Default: lasso. Other options are: plsr,
#' tlasso
#' @param insitu Intact in situ staining data
#' @return roc object
#' @import XLConnect
#' @export
eval_seurat_innter <- function(object, g = NULL, scheme = "lasso", insitu){
  if (is.null(g) || missingArg(g)) {
    cat(">> Run Seurat entire analysis\n")
  } else {
    cat(">> Evaluating performance on gene ", g , "\n")
  }
  ## read-in insitu and load into object
  full.insitu.g <- toupper(getSheets(insitu))
  insitu.matrix <- data.frame(sapply(1:length(full.insitu.g), function(i)
                                    as.numeric(as.matrix(insitu[i][2:9, 2:9]))))
  colnames(insitu.matrix) <- toupper(full.insitu.g)
  g.staining <- insitu.matrix[, g]
  insitu.genes <- setdiff(full.insitu.g, g)
  object@insitu.matrix <- insitu.matrix[, insitu.genes]

  ## imputation
  genes.sig <- pca.sig.genes(object, pcs.use = c(1,2,3), pval.cut = 1e-2,
                             use.full = TRUE)
  predictor.genes <- unique(c(genes.sig, object@var.genes))
  object <- fill_imputed_expr(object, genes.use = predictor.genes,
                              genes.fit = insitu.genes, scheme = scheme,
                              do.print = FALSE, s.use = 40, gram = FALSE)
  ## Part-2 ends here

  ## Part-3 starts
  ## bi-modal model construction
  for (i in rev(insitu.genes)) {
    object <- fit_gene_k(object, gene = i, num.iter = 1, do.k = 2,
                        start.pct = mean(object@insitu.matrix[, i]),
                        do.plot = FALSE)
  }
  ## intinal mapping
  object <- initial_mapping(object)
  ## refined mapping
  num.pc <- 3
  num.genes <- 6
  genes.use <- pcTopGenes(object, pc.use = 1:num.pc, num.genes = num.genes,
                          use.full = TRUE, do.balanced = TRUE)
  new.imputed <- genes.use[!genes.use %in% rownames(object@imputed)]
  predictor.use <- unique(c(object@var.genes,
                            pca.sig.genes(object, pcs.use = c(1,2,3),
                                          pval.cut = 1e-2, use.full = TRUE)
                            ))
  object <- fill_imputed_expr(object, genes.use = predictor.use,
                              genes.fit = new.imputed,
                              scheme = scheme,
                              do.print = FALSE, s.use = 40, gram = FALSE)
  object <- refined_mapping(object, genes.use)
  return(eval_roc_gene(object, g, g.staining))
}


#' Generate ROC for specific landmarkd gene
#'
#' @param object seurat object
#' @param g gene name
#' @param g.staining Value of gene in in situ staining data
#' @return roc object
#' @import pROC
eval_roc_gene <- function(object, g, g.staining){
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
  raw  <- g.staining
  g.roc.obj <- roc(raw, pred)
  g.roc.obj
}


plot_roc_list <- function(l, main = " "){
  nms <- names(l)
  n   <- length(nms)
  mycolor <- brewer.pal(ifelse(n < 3, 3, n), 'Set1')
  i <- 1
  auc.vals <- rep(0, n)
  for (obj in l) {
    if (i == 1){
      plot.roc(obj, col = mycolor[i], asp = 1, main = main,
               cex = 1.5, cex.lab = 1.5, cex.axis = 1.5)
    } else {
      lines.roc(obj, col = mycolor[i])
    }
    auc.vals[i] <- as.numeric(obj$auc)
    i <- i + 1
  }
  auc.vals <- format(auc.vals, digits = 3)
  legend.labs <- paste0(nms, " (", auc.vals, ")")
  legend("bottomright", legend = legend.labs, col = mycolor[1:i], lwd = 2)
  0
}

#' @title Evaluate Seurat package performance
#'
#' @description Compare the inferred and observed in-situ patterns of landmark
#' genes to evaluate Seurat performance. Refer to paper 'Evaluating Seurat's
#' perforamce' for details.
#' @param object The seurat object.
#' @param genes.eval The landmark genes to be used to evaluted performance.
#' @param dir Directory to save the ROC curve plot and violin plot.
#' @param parallel If "TRUE", \link{foreach} is used to run parallel.
#' @param scheme Imputation strategies.
#' @param insitu.path Path to in situ staining data.
#' @import pROC
#' @import RColorBrewer
#' @import vioplot
#' @import doMC
#' @import foreach
#' @export
setGeneric("eval_seurat",
  function(object, genes.eval, dir, parallel = TRUE,
           scheme = c("lasso", "plsr"), insitu.path) standardGeneric("eval_seurat"))
setMethod("eval_seurat", "seurat",
  function(object, genes.eval, dir, parallel = TRUE,
           scheme = c("lasso", "plsr"), insitu.path) {
    genes.roc.obj <- list()
    genes.eval.default <- c("ADMP", "OSR1", "CDX4", "SOX3", "CHD", "SZL")
    if (missingArg(genes.eval)){
      genes.eval <- genes.eval.default
      cat(genes.eval, " to be evaluated\n")
    }
    if (missingArg(dir)){
      dir <- getwd()
    } else {
      if (!dir.exists(dir)) dir.create(dir)
    }

    wb <- loadWorkbook(insitu.path, create = FALSE)
    if (!parallel) {
      genes.roc.obj <- lapply(genes.eval, function(x)
                              eval_seurat_innter(object, x, scheme = scheme,
                                                 insitu = wb))
    } else {
      registerDoMC(2)
      genes.roc.obj <- foreach(x = genes.eval) %dopar% {
                        eval_seurat_innter(object, x, scheme = scheme,
                                           insitu = wb)}
    }
    names(genes.roc.obj) <- genes.eval
    genes.auc.val <- sapply(genes.roc.obj, function(o) o$auc)
    auc.dec.idx <- order(genes.auc.val, na.last = TRUE, decreasing = TRUE)
    genes.roc.obj <- genes.roc.obj[auc.dec.idx]

    genes.roc.default <- genes.roc.obj[names(genes.roc.obj) %in% genes.eval.default]
    pdf(paste0(dir, '/eval_seurat1.pdf'), 7, 7)
    plot_roc_list(genes.roc.default)
    dev.off()

    if (length(genes.eval) > 6){
      pdf(paste0(dir, '/eval_seurat2.pdf'), 7, 7)
      plot_roc_list(genes.roc.obj[1:6], main = "Top 6 genes")
      dev.off()

      pdf(paste0(dir, '/eval_seurat3.pdf'), 3, 7)
      vioplot(as.numeric(genes.auc.val), names = "", ylim = c(0, 1),
              col = "tomato", rectCol = "grey")
      abline(h = median(as.numeric(genes.auc.val)), lty = 2)
      text(x = 1, y = median(as.numeric(genes.auc.val)),
           format(median(as.numeric(genes.auc.val)), digits = 3))
      dev.off()
    } else {
      pdf(paste0(dir, '/eval_seurat2.pdf'), 7, 7)
      plot_roc_list(genes.roc.obj, main = "ROC for all available landmark genes")
      dev.off()
    }
    # return(as.data.frame(genes.auc.val))
    return(genes.roc.obj)
  }
)


# #Internal, not documented for now
# setGeneric("calc.insitu", function(object,gene,do.plot=FALSE,do.write=FALSE,write.dir="~/window/insitu/",col.use=bwCols,do.norm=FALSE,cells.use=NULL,do.return=FALSE,probs.min=0,do.log=FALSE, use.imputed=FALSE, bleach.use=0) standardGeneric("calc.insitu"))
# setMethod("calc.insitu", "seurat",
          # function(object,gene,do.plot=FALSE,do.write=FALSE,write.dir="~/window/insitu/",col.use=bwCols,do.norm=FALSE,cells.use=NULL,do.return=FALSE,probs.min=0,do.log=FALSE,use.imputed=FALSE,bleach.use=0) {
            # cells.use=set.ifnull(cells.use,colnames(object@final.prob))
            # probs.use=object@final.prob
            # data.use=exp(object@data)-1
            # if (use.imputed) data.use=exp(object@imputed)-1
            # cells.use=cells.use[cells.use%in%colnames(probs.use)]; cells.use=cells.use[cells.use%in%colnames(data.use)]
            # #insilico.stain=matrix(unlist(lapply(1:64,function(x) sum(probs.use[x,]*data.use[gene,]))),nrow=8,ncol=8)
            # insilico.vector=unlist(lapply(1:64,function(x) sum(as.numeric(probs.use[x,cells.use])*as.numeric(data.use[gene,cells.use]))))
            # probs.total=apply(probs.use,1,sum)
            # probs.total[probs.total<probs.min]=probs.min
            # insilico.stain=(matrix(insilico.vector/probs.total,nrow=8,ncol=8))
            # if (do.log) insilico.stain=log(insilico.stain+1)
            # if (bleach.use > 0) {
              # insilico.stain=insilico.stain-bleach.use
              # insilico.stain=minmax(insilico.stain,min=0,max=1e6)
            # }
            # if (do.norm) insilico.stain=(insilico.stain-min(insilico.stain))/(max(insilico.stain)-min(insilico.stain))
            # title.use=gene
            # if (gene %in% colnames(object@insitu.matrix)) {
              # pred.use=prediction(insilico.vector/probs.total,object@insitu.matrix[,gene],0:1)
              # perf.use=performance(pred.use,"auc")
              # auc.use=round(perf.use@y.values[[1]],3)
              # title.use=paste(gene,sep=" ")
              # cat(gene, " ", auc.use, "\n")
            # }
            # if (do.write) {
              # write.table(insilico.stain,paste(write.dir,gene,".txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
            # }
            # if (do.plot) {
              # aheatmap(insilico.stain,Rowv=NA,Colv=NA,col=col.use, main=title.use)
            # }
            # if (do.return) {
              # return(as.vector(insilico.stain))
            # }
            # return(object)
          # }
# )


