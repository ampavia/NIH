#!/usr/bin/env Rscript

devtools::install_github("jtlovell/GENESPACE")

library(GENESPACE)

#change the working directory path to match your own
wd = "/scratch/ac05869/gs_pcu_cro/wd"
path2mcscanx = '/home/ac05869/miniconda/envs/genespace/bin'
gpar <- init_genespace( wd = wd, path2mcscanx = path2mcscanx, nCores = 32,
                        genomeIDs = c( "pcu_v0_h1" , "pcu_v0_h2", "cro_v3"  ),
                        ploidy = c( 1, 1, 1) )
out <- run_genespace(gpar, overwrite = T)
