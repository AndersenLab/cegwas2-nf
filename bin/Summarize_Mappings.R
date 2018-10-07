#!/usr/bin/env Rscript
library(tidyverse)

list.files()

all_maps <- list()
for(i in 1:length(grep("mapping.tsv",list.files(), value = T))){
  
  t_pr_map <- readr::read_tsv(grep("mapping.tsv",list.files(), value = T)[i],
                              col_types = cols(
                                marker = col_character(),
                                CHROM = col_character(),
                                POS = col_integer(),
                                log10p = col_double(),
                                trait = col_character(),
                                BF = col_character(),
                                aboveBF = col_integer(),
                                strain = col_character(),
                                value = col_double(),
                                allele = col_integer(),
                                var.exp = col_double(),
                                startPOS = col_integer(),
                                peakPOS = col_integer(),
                                endPOS = col_integer(),
                                peak_id = col_character(),
                                interval_size = col_integer()
                              ))
  
  all_maps[[i]] <- t_pr_map
  
}

pr_map_df <- dplyr::bind_rows(all_maps)


significant <- pr_map_df %>%
  na.omit()

uniquemap <- significant %>%
  distinct(trait, log10p, .keep_all = TRUE)

# uniquemap <- uniquemap[-1,]

goodtraits <- uniquemap %>%
  ungroup() %>%
  dplyr::mutate(condition = "Arsenic") %>%
  dplyr::group_by(CHROM) %>%
  dplyr::arrange(CHROM, startPOS)

goodtraits$POS <- as.numeric(goodtraits$POS)
goodtraits$startPOS <- as.numeric(goodtraits$startPOS)
goodtraits$endPOS <- as.numeric(goodtraits$endPOS)
goodtraits$peakPOS <- as.numeric(goodtraits$peakPOS)
goodtraits$log10p <- as.numeric(goodtraits$log10p)


#Set chromosome boundaries
newrows <- goodtraits[1,]
newrows$POS <- as.numeric(newrows$POS)
newrows$startPOS <- as.numeric(newrows$startPOS)
newrows$endPOS <- as.numeric(newrows$endPOS)
newrows$peakPOS <- as.numeric(newrows$peakPOS)
newrows$log10p <- as.numeric(newrows$log10p)

newrows[1,] = c(NA,"I",1,"pc1",1,NA,NA,NA,NA,NA,NA,0,NA,14972282,NA, NA, "Arsenic", NA)
newrows[2,] = c(NA,"II",1,"pc1",NA,NA,NA,NA,NA,NA,NA,0,NA,15173999,NA, NA, "Arsenic", NA)
newrows[3,] = c(NA,"III",1,"pc1",NA,NA,NA,NA,NA,NA,NA,0,NA,13829314,NA, NA, "Arsenic", NA)
newrows[4,] = c(NA,"IV",1,"pc1",NA,NA,NA,NA,NA,NA,NA,0,NA,17450860,NA, NA, "Arsenic", NA)
newrows[5,] = c(NA,"V",1,"pc1",NA,NA,NA,NA,NA,NA,NA,0,NA,20914693,NA, NA, "Arsenic", NA)
newrows[6,] = c(NA,"X",1,"pc1",NA,NA,NA,NA,NA,NA,NA,0,NA,17748731,NA, NA, "Arsenic", NA)
newrows$POS <- as.numeric(newrows$POS)
newrows$startPOS <- as.numeric(newrows$startPOS)
newrows$endPOS <- as.numeric(newrows$endPOS)
newrows$peakPOS <- as.numeric(newrows$peakPOS)
newrows$log10p <- as.numeric(newrows$log10p)
newrows$trait <- ""
newrows <- dplyr::ungroup(newrows)

goodtraits$POS <- as.numeric(goodtraits$POS)
goodtraits$startPOS <- as.numeric(goodtraits$startPOS)
goodtraits$endPOS <- as.numeric(goodtraits$endPOS)
goodtraits$peakPOS <- as.numeric(goodtraits$peakPOS)
goodtraits$log10p <- as.numeric(goodtraits$log10p)

goodtraits%>%
  dplyr::ungroup()%>%
  ggplot()+
  aes(x=POS/1E6, y=trait)+
  theme_bw() +
  scale_fill_gradient(high = "#D7263D", low = "#0072B2",
                      name = expression(-log[10](italic(p))))+
  scale_color_gradient(high = "#D7263D", low = "#0072B2",
                       name = expression(-log[10](italic(p))))+
  geom_segment(aes(x = startPOS/1e6, y = trait, xend = endPOS/1e6, yend = trait, color = log10p), size = 2, alpha = 1) +
  geom_segment(data=newrows,aes(x = 0, y = trait, xend = endPOS/1e6, yend = trait), size = 2.5, alpha = 0) +
  geom_point(aes(fill=log10p),colour = "black",size = 2, alpha = 1, shape = 25)+
  xlab("Genomic Position (Mb)") + ylab("") +
  theme(strip.background = element_rect(colour = "black", fill = "white",
                                        size = 0.75, linetype = "solid")) +
  theme_bw(15) + 
  theme(strip.background = element_blank(),
        panel.background = element_blank(),
        strip.text.y = element_blank(),
        panel.spacing = unit(0.5, "lines"))+
  facet_grid(. ~ CHROM, scales = "free", space = "free")+
  ggplot2::labs(x = "Genomic Position (Mb)")

ggsave("Summarized_mappings.pdf", height = 12, width = 12)
