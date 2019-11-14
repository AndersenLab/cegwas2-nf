#!/usr/bin/env Rscript
library(tidyverse)
library(cegwas2)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

# load trait file
traits <- readr::read_tsv(args[1])

# fix strain names
if(args[2] == "fix"){
fixed_names <- cegwas2::process_phenotypes(traits)
} else {
  fixed_names <- traits
}
  


for(i in 1:(ncol(fixed_names)-1)){
  t_df <- fixed_names[,c(1,i+1)]
  trait_name <- colnames(t_df)[2]
  write.table(t_df, 
              file = glue::glue("pr_{trait_name}.tsv"),
              quote = F, col.names = T, row.names = F, sep="\t")
}

write.table(fixed_names$strain, 
            file = glue::glue("Phenotyped_Strains.txt"),
            quote = F, col.names = F, row.names = F)
