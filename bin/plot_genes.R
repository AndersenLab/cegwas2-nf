#!/usr/bin/env Rscript
library(tidyverse)
library(cegwas2)

# input arguments
# 1 = LD file
# 2 = phenotype file
# 3 = gene file
# 4 = vcf


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

# verify that provided VCF has annotation field
test_out <- system(glue::glue("bcftools view {args[4]} | head -10000 | grep ANN"), intern=T)

if(length(test_out) > 0) {
    q_vcf <- args[4]
} else {
    system("echo 'Provided VCF does not have an ANN column, using CeNDR default'")
    if(grepl("20180527", args[5])){
        q_vcf <- glue::glue("http://storage.googleapis.com/elegansvariation.org/releases/{args[5]}/variation/WI.{args[5]}.soft-filter.vcf.gz")
    } else {
        q_vcf <- glue::glue("http://storage.googleapis.com/elegansvariation.org/releases/{args[5]}/variation/WI.{args[5]}.soft-filtered.vcf.gz")
    }
}

# update 20210312 KSE: avoid crashing by using one strain to get the snpeff annotations

snpeff_out <- list()
for(r in 1:nrow(query_regions)){
  cq <- query_regions$CHROM[r]
  sq <- query_regions$start_pos[r]
  eq <- query_regions$end_pos[r]

  snpeff_out[[r]] <- cegwas2::query_vcf(glue::glue("{cq}:{sq}-{eq}"),
                                        impact = "ALL",
                                        #samples = unique(phenotypes$strain),
                                        samples = "CB4856",
                                        vcf = q_vcf) %>%
    dplyr::select(CHROM:ALT, effect:transcript_biotype, nt_change:aa_change)%>%
    dplyr::distinct(CHROM, POS, .keep_all = T) %>%
    tidyr::unite(marker, CHROM, POS, sep = "_")
}

snpeff_df <- dplyr::bind_rows(snpeff_out) %>%
  dplyr::left_join(pr_trait_ld, ., by = c("marker", "REF", "ALT"))

genes_in_region <- gene_ref_flat %>%
  dplyr::filter(wbgene %in% snpeff_df$gene_id) %>%
  dplyr::select(gene_id = wbgene, strand, txstart, txend, feature_id = gene) %>%
  dplyr::arrange(txstart, feature_id)%>%
  dplyr::distinct(gene_id, feature_id, .keep_all = TRUE)

ugly_genes_in_region <- genes_in_region%>%
  dplyr::left_join(snpeff_df, ., by = "gene_id") %>%
  dplyr::distinct(marker, CHROM, POS, log10p, peak_marker, strain, impact, .keep_all = T)
  #dplyr::left_join(., phenotypes, by = "strain")

tidy_genes_in_region <- genes_in_region%>%
  dplyr::left_join(snpeff_df, ., by = "gene_id") %>%
  dplyr::distinct(marker, CHROM, POS, log10p, peak_marker, strain, impact, .keep_all = T) %>%
  dplyr::left_join(., phenotypes, by = "strain") %>%
  dplyr::select(MARKER = marker, CHROM, POS, REF, ALT, MAF_variant = maf_marker_b,
                GENE_NAME = gene_name, WBGeneID = gene_id, WBFeature_TYPE = feature_type,
                WBFeature_ID = feature_id.x, TRANSCRIPT_BIOTYPE = transcript_biotype, VARIANT_IMPACT = impact,
                NUCLEOTIDE_CHANGE = nt_change, AMINO_ACID_CHANGE = aa_change,
                STRAND = strand, TRANSCRIPTION_START_POS = txstart, TRANSCRIPTION_END_POS = txend,
                PEAK_MARKER = peak_marker, PEAK_MAF = peak_maf, TRAIT = trait,
                QTL_INTERVAL_START = start_pos, QTL_INTERVAL_END = end_pos,
                VARIANT_LD_WITH_PEAK_MARKER = ld_r2, VARIANT_LOG10p = log10p,
                STRAIN = strain, STRAIN_GENOTYPE = allele)

write_tsv(tidy_genes_in_region,
          path = glue::glue("{analysis_trait}_snpeff_genes.tsv"))

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
  dplyr::distinct(CHROM, POS, log10p, impact)

variant_df$impact[is.na(variant_df$impact)] <- "Intergenic"

xs <- unique(gene_df$start_pos)
xe <- unique(gene_df$end_pos)
gc <- unique(gene_df$CHROM)

max_logp <- unique(max(variant_df$log10p))/150

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
                   color = impact), data = variant_df) +
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
       filename = glue::glue("{analysis_trait}_{gc}:{xs}-{xe}_gene_plot.pdf"),
       height=10, width = 14)
}

