eval_seurat_innter <- function(object, g){
  cat(">> Evaluating gene ", g , "\n")
  ## intinal mapping excluding specific landmark gene
  object <- initial_mapping(object, gene.exclude = g)
  ## refined mapping excluding same specific landmark gene
  num.pc <- 3
  num.genes <- 6
  genes.use <- pcTopGenes(object, pc.use = 1:num.pc, num.genes = num.genes,
                          use.full = TRUE, do.balanced = TRUE)
  new.imputed <- genes.use[!genes.use %in% rownames(object@imputed)]
  predictor.use <- unique(c(object@var.genes,
                            pca.sig.genes(object, pcs.use = c(1,2,3),
                                          pval.cut = 1e-2, use.full = TRUE)
                            ))
  object <- addImputedScore(object, genes.use = predictor.use,
                            genes.fit = new.imputed,
                            do.print = FALSE, s.use = 40, gram = FALSE)
  # genes.use     <- setdiff(genes.use, g)
  object <- refined.mapping(object, genes.use)
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
#' @param parallel If "TRUE", \link{foreach} is used to run parallel.
#' @import pROC
#' @import RColorBrewer
#' @import vioplot
#' @import doMC
#' @import foreach
#' @export
setGeneric("eval_seurat",
  function(object, genes.eval, dir, parallel = TRUE) standardGeneric("eval_seurat"))
setMethod("eval_seurat", "seurat",
  function(object, genes.eval, dir, parallel = TRUE) {
    genes.roc.obj <- list()
    genes.eval.default <- c("ADMP", "OSR1", "CDX4", "SOX3", "CHD", "SZL")
    if (missingArg(genes.eval)){
      genes.eval <- genes.eval.default
      cat(genes.eval, " to be evaluated\n")
    }
    if (missingArg(dir)){
      dir <- getwd()
    }
    if (!parallel) {
      genes.roc.obj <- lapply(genes.eval, function(x)
                              eval_seurat_innter(object, x))
    } else {
      registerDoMC(2)
      genes.roc.obj <- foreach(x = genes.eval) %dopar% {
                        eval_seurat_innter(object, x)}
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


#Internal, not documented for now
setGeneric("calc.insitu", function(object,gene,do.plot=FALSE,do.write=FALSE,write.dir="~/window/insitu/",col.use=bwCols,do.norm=FALSE,cells.use=NULL,do.return=FALSE,probs.min=0,do.log=FALSE, use.imputed=FALSE, bleach.use=0) standardGeneric("calc.insitu"))
setMethod("calc.insitu", "seurat",
          function(object,gene,do.plot=FALSE,do.write=FALSE,write.dir="~/window/insitu/",col.use=bwCols,do.norm=FALSE,cells.use=NULL,do.return=FALSE,probs.min=0,do.log=FALSE,use.imputed=FALSE,bleach.use=0) {
            cells.use=set.ifnull(cells.use,colnames(object@final.prob))
            probs.use=object@final.prob
            data.use=exp(object@data)-1
            if (use.imputed) data.use=exp(object@imputed)-1
            cells.use=cells.use[cells.use%in%colnames(probs.use)]; cells.use=cells.use[cells.use%in%colnames(data.use)]
            #insilico.stain=matrix(unlist(lapply(1:64,function(x) sum(probs.use[x,]*data.use[gene,]))),nrow=8,ncol=8)
            insilico.vector=unlist(lapply(1:64,function(x) sum(as.numeric(probs.use[x,cells.use])*as.numeric(data.use[gene,cells.use]))))
            probs.total=apply(probs.use,1,sum)
            probs.total[probs.total<probs.min]=probs.min
            insilico.stain=(matrix(insilico.vector/probs.total,nrow=8,ncol=8))
            if (do.log) insilico.stain=log(insilico.stain+1)
            if (bleach.use > 0) {
              insilico.stain=insilico.stain-bleach.use
              insilico.stain=minmax(insilico.stain,min=0,max=1e6)
            }
            if (do.norm) insilico.stain=(insilico.stain-min(insilico.stain))/(max(insilico.stain)-min(insilico.stain))
            title.use=gene
            if (gene %in% colnames(object@insitu.matrix)) {
              pred.use=prediction(insilico.vector/probs.total,object@insitu.matrix[,gene],0:1)
              perf.use=performance(pred.use,"auc")
              auc.use=round(perf.use@y.values[[1]],3)
              title.use=paste(gene,sep=" ")
              cat(gene, " ", auc.use, "\n")
            }
            if (do.write) {
              write.table(insilico.stain,paste(write.dir,gene,".txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
            }
            if (do.plot) {
              aheatmap(insilico.stain,Rowv=NA,Colv=NA,col=col.use, main=title.use)
            }
            if (do.return) {
              return(as.vector(insilico.stain))
            }
            return(object)
          }
)


