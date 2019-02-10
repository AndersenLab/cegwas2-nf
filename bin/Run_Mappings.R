#!/usr/bin/env Rscript
library(tidyverse)
library(rrBLUP)
library(ggbeeswarm)

# argument information
# 1 - Genetoype matrix
# 2 - Phenotype data 
# 3 - Number of cores of parallel processing
# 4 - P3D, Boolean, F = EMMA, T = EMMAx
# 5 - Number of independent tests from Eigen Decomposition of SNV correlation matrix
# 6 - String, What significance threshold to use for defining QTL, 
#       BF = Bonferroni, 
#       EIGEN = Defined by Number of independent tests from Eigen Decomposition of SNV correlation matrix, 
#       or a user defined number
# 7 - If two QTL are less than this distance from each other, combine the QTL into one, (DEFAULT = 1000)
# 8 - Number of SNVs to the left and right of the peak marker used to define the QTL confidence interval, (DEFAULT = 150)

# load arguments
args <- commandArgs(trailingOnly = TRUE)

# define the number of independent tests used for EIGEN threshold
independent_tests <- as.numeric(args[5])
independent_test_cutoff <- -log10(0.05/independent_tests)

# Define number of cores available for parallel processing
cores_avail <- as.numeric(args[3])

# load genotype matrix
genotype_matrix <- readr::read_tsv(args[1]) %>%
  na.omit()

# load phenotpe data
phenotype_data <- readr::read_tsv(args[2]) %>%
  na.omit() %>%
  as.data.frame()

# generate kinship matrix
kinship_matrix <- rrBLUP::A.mat(t(genotype_matrix[,5:ncol(genotype_matrix)]), n.core = cores_avail)

# extract trait name from phenotype data
trait_name <- colnames(phenotype_data)[2]

# define method for setting significance threshold
significance_threshold <- args[6]

# set significance threshold
if(significance_threshold == "BF"){
  QTL_cutoff <- NA
} else if(significance_threshold == "EIGEN"){
  QTL_cutoff <- independent_test_cutoff
} else {
  QTL_cutoff <- as.numeric(args[6])
}


# mapping function
gwa_mapping <- function (data, 
                         cores = cores_avail, 
                         kin_matrix = kinship_matrix, 
                         snpset = genotype_matrix, 
                         min.MAF = 0.05,
                         p3d = args[4]) {
  x <- data

  y <- snpset %>% dplyr::mutate(marker = paste0(CHROM, "_", POS)) %>% 
    dplyr::select(marker, everything(), -REF, -ALT) %>% 
    as.data.frame()

  kin <- as.matrix(kin_matrix)

  pmap <- rrBLUP::GWAS(pheno = x, 
                       geno = y, 
                       K = kin, 
                       min.MAF = min.MAF, 
                       n.core = cores, 
                       P3D = as.logical(p3d), 
                       plot = FALSE)

  return(pmap)
}

# plotting function
manplot_edit <- function(plot_df, 
                         bf_line_color = "gray",
                         eigen_line_color = "gray",
                         eigen_cutoff = independent_test_cutoff) {
  plot_traits <- unique(plot_df$trait)
  plots <- lapply(plot_traits, function(i) {
    plot_df %>%
      dplyr::filter(trait == i,
                    CHROM != "MtDNA") %>%
      dplyr::distinct(marker, .keep_all = T) %>%
      dplyr::mutate(EIGEN_CUTOFF = eigen_cutoff) %>%
      dplyr::mutate(EIGEN_SIG = ifelse(log10p > BF, "1", 
                                       ifelse(log10p > EIGEN_CUTOFF, "2", "0")) )%>%
      ggplot2::ggplot(.) +
      ggplot2::aes(x = POS/1e6, y = log10p) +
      ggplot2::scale_color_manual(values = c("0" = "black", 
                                             "1" = "red",
                                             "2" = "hotpink3")) +
      ggplot2::geom_rect(ggplot2::aes(xmin = startPOS/1e6, 
                                      xmax = endPOS/1e6, 
                                      ymin = 0, 
                                      ymax = Inf, 
                                      fill = "blue"), 
                         color = "blue",fill = "cyan",linetype = 2, 
                         alpha=.3)+
      ggplot2::geom_hline(ggplot2::aes(yintercept = BF),
                          color = bf_line_color, 
                          alpha = .75,  
                          size = 1) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = EIGEN_CUTOFF),
                          color = eigen_line_color, 
                          alpha = .75,  
                          size = 1,
                          linetype = 2) +
      ggplot2::geom_point( ggplot2::aes(color= factor(EIGEN_SIG)) ) +
      ggplot2::facet_grid( . ~ CHROM, scales = "free_x" , space = "free_x") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16),
                     axis.text.y = ggplot2::element_text(size = 16),
                     axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                     axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                     strip.text.x = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                     strip.text.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                     plot.title = ggplot2::element_text(size = 24, face = "bold", vjust = 1), 
                     panel.background = ggplot2::element_rect(color = "black",size = 1.2),
                     legend.position = "none",
                     strip.background = element_blank()) +
      ggplot2::labs(x = "Genomic Position (Mb)",
                    y = expression(-log[10](italic(p))),
                    title = plot_traits)
  })
  plots
}

# process mapping function
process_mapping_df <- function (mapping_df, 
                                phenotype_df, 
                                CI_size = as.numeric(args[8]), 
                                snp_grouping = as.numeric(args[7]), 
                                BF = NA,
                                geno = genotype_matrix) {
  pheno <- phenotype_df 
  
  pheno$trait <- colnames(phenotype_df)[2]
  
  colnames(pheno) <- c("strain", "value", "trait")
  
  if (!is.na(BF) & BF < 1) {
    message("Taking -log10 of BF.")
    BF <- -log10(BF)
  }
  
  colnames(mapping_df) <- c("marker", "CHROM", "POS", "log10p")
  
  mapping_df <- mapping_df %>% 
    dplyr::mutate(trait = colnames(phenotype_df)[2]) %>%
    dplyr::group_by(trait) %>% 
    dplyr::filter(log10p != 0) %>% 
    dplyr::mutate(BF = ifelse(is.na(BF), -log10(0.05/sum(log10p > 0, na.rm = T)), BF)) %>% 
    dplyr::mutate(aboveBF = ifelse(log10p >= BF, 1, 0))
  
  Processed <- mapping_df %>% 
    dplyr::filter(sum(aboveBF, na.rm = T) > 0) %>% 
    dplyr::ungroup()
  
  snpsForVE <- Processed %>% 
    dplyr::filter(aboveBF == 1) %>% 
    dplyr::select(marker, trait)
  
  snpsForVE$trait <- as.character(snpsForVE$trait)
  
  if (nrow(snpsForVE) > 0) {
    
    row.names(pheno) <- gsub("-", "\\.", row.names(pheno))
    
    pheno$trait <- gsub("-", "\\.", pheno$trait)
    
    rawTr <- pheno %>%
      dplyr::left_join(., snpsForVE, by = "trait")
    
    rawTr$marker <- as.character(rawTr$marker)
    rawTr$strain <- as.character(rawTr$strain)
    
    snp_df <- geno %>% 
      dplyr::select(-REF, -ALT)
    
    gINFO <- snp_df %>% 
      dplyr::mutate(marker = paste(CHROM, POS, sep = "_")) %>% 
      dplyr::filter(marker %in% snpsForVE$marker) %>% 
      tidyr::gather(strain, allele, -marker, -CHROM, -POS)
    
    gINFO$marker <- as.character(gINFO$marker)
    gINFO <- suppressWarnings(data.frame(gINFO) %>% 
                                dplyr::left_join(., snpsForVE, by = "marker") %>%
                                dplyr::left_join(rawTr, ., by = c("trait", "strain", "marker")))
    
    cors <- gINFO %>% dplyr::group_by(trait, marker) %>% 
      dplyr::mutate(var.exp = cor(value, allele, use = "pairwise.complete.obs", 
                                  method = "spearman")^2)
    
    CORmaps <- Processed %>% 
      dplyr::left_join(., cors, by = c("trait", "marker", "CHROM", "POS"), copy = TRUE)
    processed_mapping_df <- Processed
    correlation_df <- CORmaps
    phenotypes <- as.character(unique(processed_mapping_df$trait))
    intervals <- list()
    for (i in 1:length(phenotypes)) {
      PeakDF <- processed_mapping_df %>% 
        dplyr::filter(trait ==  phenotypes[i]) %>% 
        dplyr::group_by(CHROM, trait) %>% 
        dplyr::mutate(index = 1:n()) %>% 
        dplyr::mutate(peaks = cumsum(aboveBF)) %>% 
        dplyr::filter(aboveBF == 1) %>% 
        dplyr::group_by(CHROM, trait) %>% 
        dplyr::mutate(nBF = n()) %>% dplyr::group_by(CHROM, trait) %>% 
        dplyr::arrange(CHROM, POS)
      
      SNPindex <- processed_mapping_df %>% 
        dplyr::filter(trait == phenotypes[i]) %>% 
        dplyr::group_by(CHROM, trait) %>% 
        dplyr::mutate(index = 1:n()) %>%
        dplyr::distinct(CHROM, POS, .keep_all = T) %>% 
        dplyr::select(CHROM, POS, index) %>% 
        dplyr::filter(POS == min(POS) | POS == max(POS))
      
      findPks <- PeakDF %>% 
        dplyr::filter(trait == phenotypes[i]) %>% 
        dplyr::group_by(CHROM) %>% 
        dplyr::arrange(CHROM, POS)
      
      if (findPks$nBF == 1 & length(unique(findPks$CHROM)) ==  1) {
        findPks$pID <- 1
        findPks <- findPks %>% 
          dplyr::group_by(CHROM, pID, trait) %>%
          dplyr::mutate(start = min(index) - CI_size, end = max(index) + CI_size)
        
        for (k in 1:nrow(findPks)) {
          tSNPs <- SNPindex %>% dplyr::filter(CHROM == findPks$CHROM[k])
          if (findPks$start[k] < min(tSNPs$index)) {
            findPks$start[k] <- min(tSNPs$index)
          }
          else if (findPks$end[k] > max(tSNPs$index)) {
            findPks$end[k] <- max(tSNPs$index)
          }
        }
        intervals[[i]] <- findPks %>% dplyr::ungroup()
      }
      else {
        findPks$pID <- 1
        for (j in 2:nrow(findPks)) {
          findPks$pID[j] <- ifelse(abs(findPks$index[j] - findPks$index[j - 1]) < snp_grouping & findPks$CHROM[j] == findPks$CHROM[j - 1],
                                   findPks$pID[j - 1],
                                   findPks$pID[j - 1] + 1)
        }
        findPks <- findPks %>% 
          dplyr::group_by(CHROM, pID, trait) %>% 
          dplyr::mutate(start = min(index) - CI_size, end = max(index) + CI_size)
        
        for (k in 1:nrow(findPks)) {
          tSNPs <- SNPindex %>% dplyr::filter(CHROM == findPks$CHROM[k])
          
          if (findPks$start[k] < min(tSNPs$index)) {
            findPks$start[k] <- min(tSNPs$index)
          }
          else if (findPks$end[k] > max(tSNPs$index)) {
            findPks$end[k] <- max(tSNPs$index)
          }
        }
      }
      intervals[[i]] <- findPks %>% dplyr::ungroup()
    }
    intervalDF <- data.table::rbindlist(intervals)
    peak_df <- intervalDF
    peak_list <- intervals
    Pos_Index_Reference <- processed_mapping_df %>%
      dplyr::group_by(CHROM, trait) %>% 
      dplyr::mutate(index = 1:n()) %>% 
      dplyr::mutate(peaks = cumsum(aboveBF)) %>% 
      dplyr::select(trait, CHROM, POS, index) %>% 
      dplyr::filter(index %in% c(unique(peak_df$start), unique(peak_df$end))) %>% 
      dplyr::ungroup()
    
    Pos_Index_Reference$trait <- as.character(Pos_Index_Reference$trait)
    interval_positions <- list()
    for (i in 1:length(peak_list)) {
      print(paste(100 * signif(i/length(peak_list), 3), 
                  "%", sep = ""))
      peak_list[[i]]$trait <- as.character(peak_list[[i]]$trait)
      peak_list[[i]] <- peak_list[[i]] %>% dplyr::arrange(desc(log10p)) %>% 
        dplyr::distinct(pID, .keep_all = T)
      trait_i <- unique(peak_list[[i]]$trait)
      index_i <- c(peak_list[[i]]$start, peak_list[[i]]$end)
      CHROM_i <- peak_list[[i]]$CHROM
      PKpos <- data.frame(Pos_Index_Reference) %>% 
        dplyr::filter(trait == trait_i & index %in% index_i & CHROM %in% CHROM_i) %>% 
        dplyr::left_join(., peak_list[[i]], by = c("trait",  "CHROM")) %>% 
        dplyr::mutate(issues = ifelse(start ==  index.x | end == index.x, 1, 0)) %>% 
        dplyr::filter(issues != 0) %>% 
        dplyr::select(trait, CHROM, POS.x, POS.y, pID, log10p, index.x, index.y, start, end) %>% 
        dplyr::group_by(CHROM, pID) %>% 
        dplyr::mutate(startPOS = min(POS.x),  peakPOS = POS.y, endPOS = max(POS.x)) %>%
        dplyr::distinct(trait, CHROM, pID, peakPOS, .keep_all = T) %>% dplyr::select(trait, CHROM, POS = POS.y, startPOS, peakPOS, endPOS, peak_id = pID)
      interval_positions[[i]] <- PKpos
    }
    
    interval_pos_df <- data.frame(data.table::rbindlist(interval_positions)) %>% 
      dplyr::mutate(interval_size = endPOS - startPOS)
    
    Processed <- suppressWarnings(dplyr::left_join(correlation_df, 
                                                   interval_pos_df, by = c("trait", "CHROM", "POS"), 
                                                   copy = TRUE))
  } else {
    Processed <- mapping_df %>% 
      dplyr::mutate(strain = NA, value = NA, allele = NA, var.exp = NA, 
                    startPOS = NA, peakPOS = NA, endPOS = NA, 
                    peak_id = NA, interval_size = NA)
  }
  
  return(Processed)
}
system("echo begin mapping")

# run mapping
raw_mapping <- gwa_mapping(data = phenotype_data,
                        snpset = genotype_matrix,
                        kin_matrix = kinship_matrix)

# save mapping data set
readr::write_tsv(raw_mapping, 
                 path = glue::glue("{trait_name}_raw_mapping.tsv"),
                 col_names = T)

# process mapping data, define QTL
processed_mapping <- process_mapping_df(raw_mapping, 
                                        phenotype_data, 
                                        CI_size = as.numeric(args[8]), 
                                        snp_grouping = as.numeric(args[7]), 
                                        BF = QTL_cutoff,
                                        geno = genotype_matrix)

# save processed mapping data
readr::write_tsv(processed_mapping, 
                 path = glue::glue("{trait_name}_processed_mapping.tsv"),
                 col_names = T)

# generate manhattan plot
manhattan_plot <- manplot_edit(processed_mapping)

ggsave(manhattan_plot[[1]], 
       filename = glue::glue("{trait_name}_manplot.pdf"),
       height = 4, 
       width = 12)

# generate phenotype by genotype plot
pxg_df <- na.omit(processed_mapping) 

if( nrow(pxg_df) > 0 ){
  
  pxg_df <- pxg_df %>%
    dplyr::mutate(facet_marker = paste0(CHROM, ":", peakPOS))
  
  pxg_df%>%
    dplyr::group_by(allele, facet_marker)%>%
    dplyr::mutate(mean_pheno = mean(as.numeric(value), na.rm = T))%>%
    ggplot()+
    aes(x = factor(allele, levels = c(-1,1), labels = c("REF","ALT")))+
    geom_beeswarm(cex=2,priority='density',
                  aes(y = as.numeric(value)),
                  shape = 21, 
                  fill = "gray50",
                  size = 2)+
    geom_point(aes(y = mean_pheno), 
               fill = "red", 
               size = 3, 
               shape = 25)+
    theme_bw(15)+
    facet_grid(.~facet_marker)+
    labs(y = trait_name,
         x = "Genotype") +
    theme(axis.text.x = ggplot2::element_text(size = 14),
          axis.text.y = ggplot2::element_text(size = 14),
          axis.title.x = ggplot2::element_text(size = 16, face = "bold", color = "black", vjust = -0.3),
          axis.title.y = ggplot2::element_text(size = 16, face = "bold", color = "black", vjust = -0.3),
          strip.text.x = element_text(size = 14, face = "bold"))+
    theme(legend.position = "none")
  
  plot_width_scale <- length(unique(pxg_df$facet_marker))
  
  ggsave(filename = glue::glue("{trait_name}_pxgplot.pdf"),
         height = 4, 
         width = 4*plot_width_scale, limitsize = FALSE)
}


