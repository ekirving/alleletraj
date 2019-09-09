#!/usr/bin/env Rscript
library(coda, quietly=T)
library(data.table, quietly=T)

source('rscript/mcmc_utils.R')

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
burn_perc <- as.numeric(args[1])
thin <- as.numeric(args[2])
param_files <- args[3:length(args)]

# TODO remove when done testing
# prefix <- 'horse-DOM2-modsnp15000674-n100000000-s100-h0.5'
# burn_perc <- 0.5
# thin <- 100
# param_files <- c(paste0('data/selection/', prefix, '-chain1.param.gz'),
#                  paste0('data/selection/', prefix, '-chain2.param.gz'))

# load all the chains
chains.all = load_chains(param_files, burn_perc, thin, verbose=FALSE)

mpsrf <- tryCatch({
    # values substantially above 1.2 indicate lack of convergence
    gelman.diag(chains.all, multivariate=TRUE, autoburnin=FALSE)$mpsrf
}, error = function(e) {
    # output a dummy value to indicate that no valid MPSRF could be calculated
    100
})

cat(mpsrf)
