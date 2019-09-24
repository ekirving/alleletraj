#!/usr/bin/env Rscript
library(coda, quietly=T)
library(fitR, quietly=T)
library(data.table, quietly=T)
library(rjson, quietly=T)
library(stringr, quietly=T)

source('rscript/mcmc_utils.R')

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
burn_perc <- as.numeric(args[1])
thin <- as.numeric(args[2])
ess_file  <- args[3]
psrf_file <- args[4]
trace_png <- args[5]
gelman_png <- args[6]
param_files <- args[7:length(args)]

# TODO remove when done testing
# prefix <- 'horse-DOM2-modsnp15000674-n100000000-s100-h0.5'
# burn_perc <- 0.5
# thin <- 100
# ess_file  <- paste0('data/selection/', prefix, '-chainAll.ess')
# psrf_file <- paste0('data/selection/', prefix, '-chainAll.psrf')
# trace_png <- paste0('data/pdf/selection/', prefix, '-chainAll-trace-pt1.png')
# gelman_png <- paste0('data/pdf/selection/', prefix, '-chainAll-gelman-pt1.png')
# param_files <- c(paste0('data/selection/', prefix, '-chain1.param.gz'),
#                  paste0('data/selection/', prefix, '-chain2.param.gz'))

cat("Analysing combined chains.", "\n\n")

# load all the chains
chains.all = load_chains(param_files, burn_perc, thin, drop=c('gen', 'first_nonzero'))

# print the summary stats
print(summary(chains.all))

# calculate the ESS for all params across all the replicates
ess <- effectiveSize(chains.all)
cat(toJSON(ess, indent = 2), file=ess_file)

cat("Effective Sample Size.\n")
print(ess)
cat("\n")

# plot the combined traces
cat("Plotting combined traces.", "\n\n")
png(file=str_replace(trace_png, 'pt1', 'pt%d'), width=7, height=7, units='in', res=300)
plot(chains.all)
off <- dev.off()

gelman <- tryCatch({
    # NB. values substantially above 1 indicate lack of convergence.
    gelman.diag(chains.all, multivariate=TRUE, autoburnin=FALSE)
}, error = function(e) {
    # not all models can compute a value multivariate PSRF
    gelman.diag(chains.all, multivariate=FALSE, autoburnin=FALSE)
})

gelman.list <- gelman$psrf[,1]
gelman.list['mpsrf'] <- if (is.null(gelman$mpsrf)) NA else gelman$mpsrf;
cat(toJSON(gelman.list, indent = 2), file=psrf_file)

cat("Gelman and Rubin's convergence diagnostic.", "\n")
print(gelman)

if (gelman$psrf['lnL', 1] > 1.1) {
    cat(paste0("WARNING: PSRF of likelihood above threshold = ", round(gelman$psrf['lnL', 1], 3)))
}

# do any params have infinite variance (very bad!)
inf.params <- names(gelman$psrf[,1][gelman$psrf[,1] == 'Inf'])

if (length(inf.params) > 0) {
    # drop the offending columns so we can print the Gelman plot
    chains.all <- chains.all[,!(colnames(chains.all[[1]]) %in% inf.params)]
}

cat("Plotting the Gelman and Rubin's convergence diagnostic.", "\n\n")
png(file=str_replace(gelman_png, 'pt1', 'pt%d'), width=7, height=7, units='in', res=300)
gelman.plot(chains.all)
off <- dev.off()
