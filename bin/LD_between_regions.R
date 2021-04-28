library(genetics) 
library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
gm <- read.table(file = args[1], header = T)
processed_mapping <- read.delim(args[2], stringsAsFactors=FALSE)
TRAIT <- args[3]

snp_df <- processed_mapping %>% na.omit()

ld_snps <- dplyr::filter(gm, CHROM %in% snp_df$CHROM, POS %in% snp_df$POS)


if ( nrow(ld_snps) > 1 ) {
    
    ld_snps <- data.frame(snp_id = paste(ld_snps$CHROM, ld_snps$POS,
                                         sep = "_"), data.frame(ld_snps)[, 5:ncol(ld_snps)])
    
    sn <- list()
    
    for (i in 1:nrow(ld_snps)) {
        sn[[i]] <- genetics::genotype(as.character(gsub(1, "T/T",
                                                        gsub(-1, "A/A", ld_snps[i, 4:ncol(ld_snps)]))))
    }
    
    test <- data.frame(sn)
    colnames(test) <- (ld_snps$snp_id)
    ldcalc <- t(genetics::LD(test)[[4]])^2
    diag(ldcalc) <- 1
    
    write.table(ldcalc, paste0(TRAIT, "_LD_between_QTL_regions.tsv"), quote=F, row.names = T, col.names = NA, sep="\t")
}
