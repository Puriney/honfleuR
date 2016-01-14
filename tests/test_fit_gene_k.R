source('tests/depends.R')

## source functions to be tested
source('R/bimodal_distribution.R')

## load required object
load('/Users/yunyan/Yun_Codes/Projects/seurat_dev/output_part2.Robj')

## Begin testing
outLog <- paste0(test_res_dir, '/bimodal.log')

insitu.genes <- colnames(zf@insitu.matrix)
g <- insitu.genes[1]

sink(outLog)
cat("## Testing bimodal estimationg ... ", mytime, "\n")
zf0 <- fit.gene.k(zf, g, do.plot=FALSE, do.k = 2, start.pct=mean(zf@insitu.matrix[, g]), num.iter = 1)
zf1 <- fit_gene_k(zf, g, do.plot=FALSE, do.k = 2, start.pct=mean(zf@insitu.matrix[, g]), num.iter = 1)

if (all.equal(zf0@mix.probs, zf1@mix.probs)){
  cat(">> Passing: bimodal distributions are correctly estimated ",
      mytime, "\n")
  cat("## Evaluating performance improved ... ", mytime, "\n")
  set.seed(1234)
  genes <- sample(insitu.genes, floor(length(insitu.genes) * 0.2), replace = F)
  cat("-- using random insitu genes for evaluation: ", genes, "\n")
  res <- microbenchmark(
    seurat = for(i in rev(genes)) zf=fit.gene.k(zf,i,do.plot=FALSE,do.k = 2,start.pct=mean(zf@insitu.matrix[,i]),num.iter = 1),
    seurat2 = for(i in rev(genes)) zf=fit_gene_k(zf,i,do.plot=FALSE,do.k = 2,start.pct=mean(zf@insitu.matrix[,i]),num.iter = 1),
    times = 5
  )
  print(res)
  cat("-- generating pdf", mytime, "\n")
  pdf(paste0(test_res_dir, '/bimodal.pdf'), 7, 7)
  p <- autoplot(res) +
    ggtitle('Estimating landmarked genes bimodal')
  print(p)
  dev.off()
} else{
  warning("xx Fail: estimated bimodal distribution is different from paper results ", mytime)
}

cat("## Finish testing ", mytime, "\n")
sink()
