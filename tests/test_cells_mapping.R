## pkg
library(Seurat)
library(mixtools)
library(dplyr)
library(tidyr)
library(microbenchmark)

if (!exists("run_perf")){
  assign("run_perf", FALSE)
}

## src the codes to be tested
source("R/map_cell_location.R")
source("R/utility.R")

outLog <- "./tests/fig_cells_mapping/cells_mapping.log"
theme_set(theme_gray(base_size = 22))

##
## Load seurat object with fitted
##
load("/Users/yunyan/Yun_Codes/Projects/seurat_dev/output_part3_justFitted.Robj")

##
## Part-1: Test initial mapping
##

## 1.1 :Test just one cell mapping: `map_cell` func
onecell <- colnames(zf@data)[100]
onecell.prob0 <- map.cell(zf, onecell, FALSE, FALSE)
onecell.prob1 <- map_cell(zf, onecell, FALSE, FALSE)

sink(outLog)
if (all.equal(onecell.prob0, onecell.prob1)) {
  cat(">> Pass: map_cell works exactly same as map.cell \n")
} else {
  warning("xx Fail: map_cell works exactly same as map.cell")
}
sink()


## 1.2 :Test `initial_mapping` func
zf0 <- initial.mapping(zf)
zf1 <- initial_mapping(zf)
cells.name <- colnames(zf@data)
cells.test <- sample(cells.name, length(cells.name) * 0.2,
                     replace = FALSE)
sink(file = outLog, append = T)
if (all.equal(zf0@final.prob, zf1@final.prob)){
  cat(">> Pass: initial_mapping works\n")
  if (run_perf){
    cat("-- Estimate Performance Improvement ...\n")
    res <- microbenchmark(
      seurat = initial.mapping(zf, cells.use = cells.test),
      serrat2 = initial_mapping(zf, cells.use = cells.test),
      times = 10
    )
    print(res)
    pdf("./tests/fig_cells_mapping/init_mapping.pdf", 7, 7)
    p <- autoplot(res) +
      ggtitle(paste0('Initial mapping on random ',
                     length(cells.test), " cells"))
    print(p)
    dev.off()
  }
  zf <- zf1
} else {
  warning("xx Fail: initial_mapping fail")
}
sink()
##
## Part-2 Refined Mapping
##

## Preliminary settings
num.pc=3; num.genes=6
genes.use=pcTopGenes(zf,pc.use = 1:num.pc,num.genes = num.genes,use.full = TRUE,do.balanced = TRUE)
new.imputed=genes.use[!genes.use%in%rownames(zf@imputed)]
lasso.genes.use=unique(c(zf@var.genes,pca.sig.genes(zf,pcs.use = c(1,2,3), pval.cut = 1e-2, use.full = TRUE)))
zf <- addImputedScore(zf, genes.use=lasso.genes.use,genes.fit=new.imputed, do.print=FALSE, s.use=40, gram=FALSE)

zf0 <- refined.mapping(zf, genes.use)
zf1 <- refined_mapping(zf, genes.use)
sink(outLog, append = T)
if (all.equal(zf0@final.prob, zf1@final.prob)){
  cat(">> Pass: refined_mapping works\n")
  if (run_perf){
    cat("-- Estimate Performance Improvement ...\n")
    res <- microbenchmark(
      seurat = refined.mapping(zf, genes.use),
      serrat2 = refined_mapping(zf, genes.use),
      times = 5
    )
    print(res)
    pdf("./tests/fig_cells_mapping/refined_mapping.pdf", 7, 7)
    p <- autoplot(res) +
      ggtitle('Refined mapping')
    print(p)
    dev.off()
  }
  zf <- zf1
} else {
  warning("xx Fail: refined_mapping fail")
}
sink()