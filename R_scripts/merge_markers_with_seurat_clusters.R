#!/usr/bin/env Rscript

# Find kratom stem (MIT_AD_AE) markers for cell type annotation using mitsp_v1 genome and annotation
# MIT_AD_AE_clean_integrated_all_markers_minpct0.25_minlfc0.25_with_info_30jan25.txt is the marker list from seurat with gene info
# Use blastp file of top 4 kratom hits to poplar marker genes (bestHits_pop_mitsp_v1_markers.34828282.tsv)
# Used A. thaliana genes with peptide sequences downloaded from TAIR database.
# gene_descr is the description of Ath genes from TAIR
# pop_markers_tissues is the poplar tissue types that correspond to the marker. Obtained from Josh and Anita for leaf and stem tissues

library(tidyverse)
library(Seurat)
library(readxl)
library(RColorBrewer)
library(rcartocolor)
library(viridis)
library(patchwork)

setwd("/Users/anne_collart/Desktop/pop_marker_genes")

gene_descr <- read_delim("./TAIR_marker_descriptions.txt", delim = "\t", col_names = T)  %>%
  select(c(1, 3, 5)) #, col_names = c("ath", "type", "description"))
head(gene_descr)

tissues <- read_delim("./pop_markers_tissues.txt", delim = "\t", col_names = T)
head(tissues)

pop_markers <- read_delim("./bestHits_pop_mitsp_v1_markers.34841444.tsv", delim = "\t", col_names = c("Locus", "gene_ID", "pident",	"length",	"mismatch",	"gapopen",	"qstart",	"qend",	"sstart",	"send",	"evalue",	"bitscore"))
pop_markers$gene_ID <- gsub(pattern = "_", replacement = "-", pop_markers$gene_ID)
pop_markers$gene_ID <- gsub("\\.\\d+$", "", pop_markers$gene_ID)
pop_markers$Locus <- gsub("\\.\\d+$", "", pop_markers$Locus)

pop_markers <- pop_markers %>% # Replace the IPAP kratom genes from the blastp arabidopsis hits with the blastp hits from C. roseus
  mutate(gene_ID = ifelse(Locus == "AT4G15560", "Mitsp.v1.05-1G015600", gene_ID)) %>%
  mutate(gene_ID = ifelse(Locus == "AT5G62790","Mitsp.v1.15-1G014960", gene_ID)) %>%
  mutate(gene_ID = ifelse(Locus == "AT2G02500", "Mitsp.v1.07-1G017160", gene_ID)) %>%
  mutate(gene_ID = ifelse(Locus == "AT2G26930", "Mitsp.v1.14-2G007240", gene_ID))

head(pop_markers)

cluster_markers <- read_delim("./MIT_AD_AE_clean_integrated_all_markers_minpct0.25_minlfc0.25_with_info_30jan25.txt", delim = "\t", col_names = T) %>% rename("gene_ID" = "ID")
cluster_markers$gene_ID <- gsub(pattern = "_", replacement = "-", cluster_markers$gene_ID)
head(cluster_markers)

cluster_markers2 <- inner_join(cluster_markers, pop_markers, by = "gene_ID", relationship = "many-to-many") %>%
  inner_join(tissues, by = "Locus", relationship = "many-to-many") %>%
  inner_join(gene_descr, by = "Locus", relationship = "many-to-many") %>%
  select(1, 23, 22, 27, 10, 2, 8:9, 11, 21, 5:7, 25, 26) %>%
  group_by(cluster) %>%
  distinct(gene_ID, .keep_all = TRUE)
head(cluster_markers2)

write.table(cluster_markers2, file = "./MIT_AD_AE_seurat_markers_with_pop_markers_14feb25.txt", quote = F, col.names = T)
write.csv(cluster_markers2, file = "./MIT_AD_AE_seurat_markers_with_pop_markers_14feb25.csv")
