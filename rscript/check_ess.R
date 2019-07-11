#!/usr/bin/env Rscript
library(coda, quietly=T)
library(fitR, quietly=T)
library(data.table, quietly=T)
library(rjson, quietly=T)
library(stringr, quietly=T)

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
param_file <- args[1]
burn_perc <- as.numeric(args[2])

# TODO remove when done testing
# prefix <- 'horse-DOM2-modsnp9948933-n50000000-s1000-h0.5-chain2'
# param_file <- paste0('data/selection/', prefix, '.param')
# burn_perc <- 0.5

# load the chain and convert to an MCMC object
mcmc.chain <- mcmc(fread(param_file, header = T, sep = '\t', drop=c('gen')))
chain.length <- nrow(mcmc.chain)

# fix infinite errors
mcmc.chain[,'lnL'][mcmc.chain[,'lnL'] == '-Inf'] <- -.Machine$integer.max

# convert burn % to number of records
burnin = burn_perc * chain.length

# burn in the chain (thinning is already done)
mcmc.burn <- mcmc(
    burnAndThin(mcmc.chain, burn = burnin),
    start = burnin)

# calculate the ESS for all params
ess <- effectiveSize(mcmc.burn)
cat(toJSON(ess, indent = 2))
