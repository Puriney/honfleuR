source('tests/depends.R')
library("doMC")
library("foreach")

source("R/map_cell_location.R")
source("R/eval_perf.R")

load('/Users/yunyan/Yun_Codes/Projects/seurat_dev/output_part3.Robj')

genes.eval <- colnames(zf@insitu.matrix)

roc.list.original <- eval_seurat(zf, genes.eval = genes.eval,
                            dir = test_res_dir, parallel = T)
save(roc.list.orignal, file = paste0(test_res_dir, "/roc_list_original.Rdata"))
