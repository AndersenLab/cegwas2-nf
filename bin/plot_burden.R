#!/usr/bin/env Rscript
library(tidyverse)

# args:
# 1 - trait name
# 2 - SKAT output
# 3 - VTprice output
args = commandArgs(trailingOnly=TRUE)

trait_name <- args[1]

# PLOT SKAT
skat <- data.table::fread(args[2])%>%
  dplyr::rowwise()%>%
  dplyr::mutate(CHROM = strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][1],
                POS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][1]),
                endPOS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][2]))%>%
  dplyr::mutate(size = abs(POS-endPOS)) %>%
  dplyr::filter(CHROM!="MtDNA")%>%
  dplyr::ungroup()%>%
  dplyr::filter(NumVar > 1)%>%
  dplyr::mutate(significant = ifelse(Pvalue < .05/n(), TRUE,FALSE ))%>%
  dplyr::arrange(Pvalue)

skat_plot <- skat%>%
    ggplot()+
    aes(x = POS/1e6, y = -log10(Pvalue), alpha = 0.5, color = significant)+
    geom_point()+
    scale_color_manual(values=c("black","red"))+
    facet_grid(.~CHROM, scales = "free", space = "free")+
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
                   strip.background = element_blank())+
    labs(x = "Genomic Position (Mb)", 
         y = expression(-log[10](italic(p))))

ggsave(filename = glue::glue("{trait_name}_SKAT.pdf"), 
       plot = skat_plot,
       height = 4, 
       width = 12)

# PLOT VTprice
benreg_burden <- read.table(args[3], header = T)
benreg_burden$Gene <- as.character(benreg_burden$Gene)
benreg_burden$RANGE <- as.character(benreg_burden$RANGE)

benreg_burden<-benreg_burden%>%
  dplyr::rowwise()%>%
  dplyr::mutate(CHROM = strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][1],
                POS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][1]),
                endPOS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][2]))%>%
  dplyr::mutate(size = abs(POS-endPOS)) %>%
  ungroup()%>%
  dplyr::mutate(significant = ifelse(PermPvalue < .05/n(), TRUE,FALSE ))%>%
  dplyr::filter(CHROM!="MtDNA",NumVar>1,size >500)

burden_manplot <- benreg_burden%>%
  ggplot()+
  aes(x = POS/1e6, y = Stat, size = NumVar, alpha = 0.5, color = significant)+
  geom_point(size = 0.25)+
  scale_color_manual(values=c("black","red"))+
  facet_grid(.~CHROM, scales = "free", space = "free")+
  theme_bw()+
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
                 strip.background = element_blank())+
  labs(x = "Genomic Position (Mb)", y = "Test Statisitic")

ggsave(filename = glue::glue("{trait_name}_VTprice.pdf"), 
       plot = burden_manplot,
       height = 4, 
       width = 12)