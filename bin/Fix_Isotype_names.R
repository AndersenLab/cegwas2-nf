#!/usr/bin/env Rscript
library(tidyverse)
library(cegwas2)
library(data.table)


for(i in grep(".tsv", list.files(), value = T)){
  
  temp_file <- data.table::fread(i, header = T)
  
  trt <- colnames(temp_file)[2]

  pr_trt <- cegwas2::process_phenotypes(temp_file)

  system(glue::glue("rm {i}"))
  
  readr::write_tsv(pr_trt, glue::glue("{trt}.tsv"))
}