library(Rcpp)
library(tidyverse)
library(correlateR)
library(RSpectra)

args <- commandArgs(trailingOnly = TRUE)

phenotyped_strain_snps <- readr::read_tsv(args[1]) %>%
  na.omit()

analysis_chromosome <- args[2]

chrom_geno <- list()
for(chrom in 1:length(unique(phenotyped_strain_snps$CHROM))){
  t_chrom <- unique(phenotyped_strain_snps$CHROM)[chrom]
  t_df <- dplyr::filter(phenotyped_strain_snps, CHROM == t_chrom)
  t_df <- t_df[,5:ncol(t_df)]
  
  keepMarkers <- data.frame(
    MAF = apply(t_df,
                MARGIN = 1,
                FUN = function(x){
                  x[x==-1] <- 0
                  return(sum(x, na.rm = T)/length(x))}
    )) 
  
  t_df <- dplyr::mutate(t_df, MAF = keepMarkers$MAF) %>%
    dplyr::filter(MAF >= 0.05, MAF <= 0.95) %>%
    dplyr::select(-MAF)
  
  chrom_geno[[chrom]] <- t_df
}

lapply(chrom_geno, nrow)

eigenvalues <- list()
for(chrom in 1:length(chrom_geno)){
  
  snpcor <- correlateR::cor(t(chrom_geno[[chrom]]))
  
  snpeigen <- eigs_sym(snpcor, 100, opts = list(retvec = FALSE))
  
  snpeigen$values[snpeigen$values>1] <- 1
  
  eigenvalues[[chrom]] <- sum(snpeigen$values)
}

evals <- unlist(eigenvalues)

independent_tests <- data.frame(independent_tests = sum(evals))

readr::write_csv(independent_tests, path = glue::glue("{analysis_chromosome}_independent_snvs.csv"))
