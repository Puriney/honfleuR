## pkg
library(Seurat)
library(mixtools)
library(dplyr)
library(tidyr)
library(microbenchmark)
library(pROC)
library(RColorBrewer)
library(vioplot)

mytime <- format(Sys.time(), "%Y-%m-%d %X %Z")
source('pkg/R/utility.R')

test_res_dir <- 'tests/fig_log'
if (!dir.exists(test_res_dir)){
  dir.create(test_res_dir)
}

theme_set(theme_bw(base_size = 22))
