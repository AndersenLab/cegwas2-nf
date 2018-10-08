library(Rcpp)
library(tidyverse)
library(correlateR)
library(RSpectra)

args <- commandArgs(trailingOnly = TRUE)

phenotyped_strain_snps <- readr::read_tsv(args[1]) %>%
  na.omit()

analysis_chromosome <- args[2]

strain_count <- (ncol(phenotyped_strain_snps) - 4)

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
  
  # initialize parameters for while loop
  test_eigenvalues <- 100
  significant_eigenvalues <- test_eigenvalues
  eigen_increment <- 100
  
  while (significant_eigenvalues == test_eigenvalues) {
    
    # first test 200 eigenvalues
    test_eigenvalues <- test_eigenvalues + eigen_increment
    
    print(glue::glue("Testing {test_eigenvalues} Eigenvalues"))
    
    snpcor <- correlateR::cor(t(chrom_geno[[chrom]]))
    snpeigen <- eigs_sym(snpcor, test_eigenvalues, opts = list(retvec = FALSE))
    snpeigen$values[snpeigen$values>1] <- 1
    
    significant_eigenvalues <- sum(snpeigen$values)
    
    if(significant_eigenvalues == test_eigenvalues){
      print(glue::glue("Significant Eigenvalues - {significant_eigenvalues} Exceeds or is Equal to the Number Tested {test_eigenvalues}\n Testing {eigen_increment} More"))
    }
    
  } # end while loop
  
  eigenvalues[[chrom]] <- significant_eigenvalues
} # end chrom loop 

evals <- unlist(eigenvalues)

independent_tests <- data.frame(independent_tests = sum(evals))

readr::write_csv(independent_tests, path = glue::glue("{analysis_chromosome}_independent_snvs.csv"))
