## pkg
library(Seurat)
library(mixtools)
library(dplyr)
library(tidyr)
library(microbenchmark)

mytime <- format(Sys.time(), "%Y-%m-%d %X %Z")
source('R/utility.R')

test_res_dir <- 'tests/fig_log'
if (!dir.exists(test_res_dir)){
  dir.create(test_res_dir)
}

theme_set(theme_gray(base_size = 22))