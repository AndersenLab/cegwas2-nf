#!/usr/bin/env Rscript
library(tidyverse)
# library(cegwas2)

# input arguments
# 1 = LD file
# 2 = phenotype file
# 3 = gene file
# 4 = annotation file

args <- commandArgs(trailingOnly = TRUE)

pr_trait_ld <- data.table::fread(args[1])
phenotypes <- readr::read_tsv(args[2])

load(file = args[3])

analysis_trait <- colnames(phenotypes)[2]

colnames(phenotypes) <- c("strain", "Phenotype_Value")

query_regions <- pr_trait_ld %>%
    dplyr::select(CHROM, start_pos, end_pos)%>%
    dplyr::distinct()

query_regions

# update 20210330 KSE: use impute vcf for genotypes and hard vcf annotation file for annotations
# this could be bcsq or snpeff
annotations <- data.table::fread(args[4])

# annotation type -snpeff or bcsq
if("CONSEQUENCE" %in% names(annotations)) {
    ann_type <- "bcsq"
} else {
    ann_type <- "snpeff"
}

# do this for each QTL separately, then combine
annotation_out <- list()
for(r in 1:nrow(query_regions)){
    cq <- query_regions$CHROM[r]
    sq <- query_regions$start_pos[r]
    eq <- query_regions$end_pos[r]
    
    # pull variants from finemap impute -- don't need this?
    # impute_vcf <- data.table::fread(args[4]) %>%
    #     dplyr::select(marker:POS, REF:strains)
    
    # filter annotations to include variants within region
    # this excludes variants that are not annotated -- variants not within a gene
    annotation_out[[r]] <- annotations %>%
        dplyr::filter(CHROM == cq,
                      POS >= sq,
                      POS <= eq) %>% {
                          if(ann_type == "bcsq") dplyr::select(., -Strains) else .
                      } %>%
        # tidyr::separate_rows(Strains) %>% # is this going to crash R?
        # dplyr::filter(Strains %in% strains) %>%
        # dplyr::group_by_at(vars(-Strains)) %>%
        # dplyr::summarize(Strains = paste(Strains, collapse = ",")) %>%
        # dplyr::rename(strains = Strains) %>%
        tidyr::unite(marker, CHROM, POS, sep = "_")
    
    # add variants with no annotations???
    
}

# combine annotations for regions
annotation_df <- dplyr::bind_rows(annotation_out) %>%
    dplyr::left_join(pr_trait_ld, ., by = c("marker", "REF", "ALT")) %>% {
        if(ann_type == "bcsq") dplyr::rename(., gene_id = WORMBASE_ID) else .
    }

genes_in_region <- gene_ref_flat %>%
    dplyr::filter(wbgene %in% annotation_df$gene_id) %>%
    dplyr::select(gene_id = wbgene, strand, txstart, txend, feature_id = gene) %>%
    dplyr::arrange(txstart, feature_id)%>%
    dplyr::distinct(gene_id, feature_id, .keep_all = TRUE)

ugly_genes_in_region <- genes_in_region %>%
    dplyr::left_join(annotation_df, ., by = "gene_id") %>%
    dplyr::distinct(marker, CHROM, POS, log10p, peak_marker, strain, .keep_all = T) %>% {
        if(ann_type == "snpeff") dplyr::rename(., VARIANT_IMPACT = impact) else .
    }

tidy_genes_in_region <- ugly_genes_in_region %>% {
    if(ann_type == "bcsq") {
        # no gene_name, feature_type
        dplyr::select(., MARKER = marker, CHROM, POS, REF, ALT, MAF_variant = maf_marker_b,
                      WBGeneID = gene_id, GENE,
                      WBFeature_ID = feature_id, TRANSCRIPT_BIOTYPE = BIOTYPE, CONSEQUENCE, VARIANT_IMPACT,
                      NUCLEOTIDE_CHANGE = DNA_CHANGE, AMINO_ACID_CHANGE, BLOSUM, Grantham, Percent_Protein,
                      STRAND, TRANSCRIPTION_START_POS = txstart, TRANSCRIPTION_END_POS = txend,
                      PEAK_MARKER = peak_marker, PEAK_MAF = peak_maf, TRAIT = trait,
                      QTL_INTERVAL_START = start_pos, QTL_INTERVAL_END = end_pos,
                      VARIANT_LD_WITH_PEAK_MARKER = ld_r2, VARIANT_LOG10p = log10p,
                      STRAIN = strain, STRAIN_GENOTYPE = allele, DIVERGENT)
    } else {
        dplyr::select(., MARKER = marker, CHROM, POS, REF, ALT, MAF_variant = maf_marker_b,
                      GENE_NAME = gene_name, WBGeneID = gene_id, WBFeature_TYPE = feature_type,
                      WBFeature_ID = feature_id.x, VARIANT_IMPACT,
                      NUCLEOTIDE_CHANGE = DNA_CHANGE, AMINO_ACID_CHANGE,
                      STRAND = strand, TRANSCRIPTION_START_POS = txstart, TRANSCRIPTION_END_POS = txend,
                      PEAK_MARKER = peak_marker, PEAK_MAF = peak_maf, TRAIT = trait,
                      QTL_INTERVAL_START = start_pos, QTL_INTERVAL_END = end_pos,
                      VARIANT_LD_WITH_PEAK_MARKER = ld_r2, VARIANT_LOG10p = log10p,
                      STRAIN = strain, STRAIN_GENOTYPE = allele)
    }
}

write_tsv(tidy_genes_in_region,
          path = glue::glue("{analysis_trait}_{ann_type}_genes.tsv"))

for(r in 1:length(unique(ugly_genes_in_region$start_pos))){
    
    gene_df <- ugly_genes_in_region %>%
        dplyr::filter(start_pos == unique(ugly_genes_in_region$start_pos)[r]) %>%
        dplyr::distinct(gene_id, peak_marker, CHROM, strand, txstart, txend, start_pos, end_pos, log10p) %>%
        dplyr::group_by(gene_id) %>%
        dplyr::mutate(log10p = max(log10p, na.rm = T)) %>%
        dplyr::distinct()
    
    peak_variant <- as.numeric(strsplit(unique(gene_df$peak_marker), split = ":")[[1]][2])
    
    variant_df <- ugly_genes_in_region %>%
        dplyr::filter(start_pos == unique(ugly_genes_in_region$start_pos)[r]) %>%
        dplyr::distinct(CHROM, POS, log10p, VARIANT_IMPACT)
    
    variant_df$VARIANT_IMPACT[is.na(variant_df$VARIANT_IMPACT)] <- "Intergenic"
    
    xs <- unique(gene_df$start_pos)
    xe <- unique(gene_df$end_pos)
    gc <- unique(gene_df$CHROM)
    
    max_logp <- unique(max(variant_df$log10p, na.rm = T))/150
    
    gene_plot <- ggplot(gene_df) +
        geom_vline(aes(xintercept = peak_variant/1e6),
                   linetype=3, color = "cyan")+
        geom_segment(aes(x = ifelse(strand == "+", txstart/1e6, txend/1e6),
                         xend = ifelse(strand == "+", txend/1e6, txstart/1e6),
                         y = log10p,
                         yend = log10p),
                     arrow = arrow(length = unit(5, "points")), size = 1) +
        geom_segment(aes(x = POS/1e6,
                         xend = POS/1e6,
                         y = log10p+max_logp,
                         yend = log10p-max_logp,
                         color = VARIANT_IMPACT), data = variant_df) +
        scale_color_manual(values = c("MODIFIER" = "gray50",
                                      "LOW" = "gray30",
                                      "MODERATE" = "orange",
                                      "HIGH" = "red",
                                      "Intergenic" = "gray80"),
                           breaks = c("HIGH", "MODERATE", "LOW", "MODIFIER", "Intergenic"),
                           name = "EFFECT")+
        labs(x = "Genomic Position (Mb)",
             y = expression(-log[10](italic(p))))+
        theme_bw(18)+
        xlim(c(xs/1e6, xe/1e6)) +
        theme(legend.position = "top",
              panel.grid = element_blank())
    
    ggsave(gene_plot,
           filename = glue::glue("{analysis_trait}_{gc}:{xs}-{xe}_gene_plot_{ann_type}.pdf"),
           height=10, width = 14)
}
