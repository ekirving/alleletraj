#!/usr/bin/env Rscript
library(coda, quietly=T)
library(fitR, quietly=T)
library(data.table, quietly=T)
library(rjson, quietly=T)

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
burn_perc <- as.numeric(args[1])
thin <- as.numeric(args[2])
ess_file  <- args[3]
psrf_file <- args[4]
trace_pdf <- args[5]
gelman_pdf <- args[6]
param_files <- args[7:length(args)]

# TODO remove when done testing
# prefix <- 'horse-DOM2-modsnp10654592-n50000000-s1000-h0.5'
# burn_perc <- 0.2
# thin <- 1000
# ess_file  <- paste0('data/selection/', prefix, '-chainAll.ess')
# psrf_file <- paste0('data/selection/', prefix, '-chainAll.psrf')
# trace_pdf <- paste0('data/pdf/selection/', prefix, '-chainAll-trace.pdf')
# # trace_png <- paste0('data/pdf/selection/', prefix, '-chainAll-trace-pt%d.png')
# gelman_pdf <- paste0('data/pdf/selection/', prefix, '-chainAll-gelman.pdf')
# param_files <- c(paste0('data/selection/', prefix, '-chain1.param.gz'),
#                  paste0('data/selection/', prefix, '-chain2.param.gz'),
#                  paste0('data/selection/', prefix, '-chain3.param.gz'),
#                  paste0('data/selection/', prefix, '-chain4.param.gz'))

cat("Analysing combined chains.", "\n\n")

chains <- c()

for (i in param_files) {
    # load the chain
    mcmc.chain <- mcmc(fread(i, header = T, sep = '\t', drop = c('gen', 'first_nonzero')))

    # fix infinite errors
    mcmc.chain[,'lnL'][mcmc.chain[,'lnL'] == '-Inf'] <- -.Machine$integer.max

    # convert burn % to number of records
    burnin = burn_perc * nrow(mcmc.chain)

    # burn in the chain (thinning is already done)
    chains[[i]] <- mcmc(
        burnAndThin(mcmc.chain, burn = burnin),
        start = burnin * thin,
        thin = thin)
}

chains.all = mcmc.list(chains)

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
pdf(file=trace_pdf)
# png(file=trace_png, width=7, height=7, units='in', res=300)
plot(chains.all)
off <- dev.off()

cat("Plotting the Gelman and Rubin's convergence diagnostic.", "\n\n")
pdf(file=gelman_pdf)
gelman.plot(chains.all)
off <- dev.off()

# NB. values substantially above 1 indicate lack of convergence.
gelman <- gelman.diag(chains.all, multivariate=TRUE, autoburnin=FALSE)

gelman.list <- gelman$psrf[,1]
gelman.list['mpsrf'] <- gelman$mpsrf
cat(toJSON(gelman.list, indent = 2), file=psrf_file)

cat("Gelman and Rubin's convergence diagnostic.", "\n")
print(gelman)

if (gelman$psrf['lnL', 1] > 1.1) {
    cat(paste0("WARNING: PSRF of likelihood above threshold = ", round(gelman$psrf['lnL', 1], 3)))
}
