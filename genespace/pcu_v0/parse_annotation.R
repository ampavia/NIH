#!/usr/bin/env Rscript

devtools::install_github("jtlovell/GENESPACE")

library(GENESPACE)

#change the working directory path to match your own
wd = "/scratch/ac05869/gs_pcu_cro/wd"
path2mcscanx = '/home/ac05869/miniconda/envs/genespace/bin'


parsedPaths <- parse_annotations(
  rawGenomeRepo = "/scratch/ac05869/gs_pcu_cro/genomeRepo", 
  genomeDirs = c("pcu_v0_h1", "pcu_v0_h2"),
  genomeIDs = c("pcu_v0_h1", "pcu_v0_h2"),
  gffString = "gff3",
  faString = "fa",
  headerEntryIndex = 1, 
  gffIdColumn = "ID",
  genespaceWd = "/scratch/ac05869/gs_pcu_cro/wd")
