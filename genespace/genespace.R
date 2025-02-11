#!/usr/bin/env Rscript

devtools::install_github("jtlovell/GENESPACE")

library(GENESPACE)

#change the working directory path to match your own
wd = "/scratch/ac05869/kratom_sc_gs/wd"
path2mcscanx = '/home/ac05869/miniconda/envs/genespace/bin'
gpar <- init_genespace( wd = wd, path2mcscanx = path2mcscanx, nCores = 32,
                        genomeIDs = c( "mitr_v1" , "cro_v3"  ),
                        ploidy = c( 4, 1) )
out <- run_genespace(gpar, overwrite = T)
