#!/usr/bin/env Rscript
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

# load bulk input file
pheno <- readr::read_tsv(args[1])

# iterate through each trait and make a ped output
for(i in 1:(ncol(pheno)-1)){
  # subset bulk input
  t_df <- pheno[,c(1,i+1)]
  # extract strain name
  trait_name <- colnames(t_df)[2]
  # format df to ped
  trait_to_ped <- t_df %>%
    dplyr::mutate(Fam = "elegans", Sample = strain, Paternal = 0, Maternal = 0, Sex = 2)%>%
    dplyr::select(Fam:Sex, trait_name)
  # convert NA to -9
  trait_to_ped[is.na(trait_to_ped)] <- -9
  # save ped
  write.table(trait_to_ped, 
              file = glue::glue("{trait_name}.ped"),
              quote = F, col.names = F, row.names = F)
}