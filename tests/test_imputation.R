source("tests/depends.R")
library(pls)
library(lars)

source("R/impute_expr.R")
load('/Users/yunyan/Yun_Codes/Projects/seurat_dev/output_part2.Robj')

genes.sig <- pca.sig.genes(zf,pcs.use = c(1,2,3), pval.cut = 1e-2, use.full = TRUE)
insitu.genes <- colnames(zf@insitu.matrix)

imputed.dt0 <- zf@imputed
lasso.genes.use <- unique(c(genes.sig,zf@var.genes))

cat("## Testing fill_imputed_expr - lasso: correct? \n")
zf1 <- fill_imputed_expr(zf, genes.use = lasso.genes.use, genes.fit = insitu.genes, scheme = "lasso", do.print=FALSE, s.use=40, gram=FALSE)
imputed.dt1 <- zf1@imputed

if (all.equal(imputed.dt0, imputed.dt1)){
  cat(">> Pass lasso scheme works correctly as original seurat\n")
} else {
  cat("xx Fail !!!")
}
cat("## Testing fill_imputed_expr - plsr: works?")
zf2 <- fill_imputed_expr(zf, genes.use = lasso.genes.use, genes.fit = insitu.genes, scheme = "plsr", do.print=FALSE)

cat("## Demonstrate the benifit of imputation")
genePlot(zf,"MIXL1","OSR1",use.imputed = FALSE,col="grey",cex.use = 1)
genePlot(zf,"MIXL1","OSR1",use.imputed = TRUE,col="black",cex.use = 1)
genePlot(zf1,"MIXL1","OSR1",use.imputed = TRUE,col="blue",cex.use = 1)
genePlot(zf2,"MIXL1","OSR1",use.imputed = TRUE,col="brown",cex.use = 1)
