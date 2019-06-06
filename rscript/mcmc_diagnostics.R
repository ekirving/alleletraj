#!/usr/bin/env Rscript
library(coda, quietly=T)
library(fitR, quietly=T)
library(data.table, quietly=T)

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
burn.perc<- args[2]

# TODO remove when done testing
prefix <- 'horse-DOM2-modsnp9016431-chain1-n50000-s100-h0.5'
burn.perc<- 0.2

cat("Loading MCMC...\n")
cat("Prefix: ", prefix, "\n")

param_file <- paste0('data/selection/', prefix, '.param')

# load the chain and convert to an MCMC object
mcmc.chain <- mcmc(fread(param_file, header = T, sep = '\t', drop=c('gen')))

chain.length <- nrow(mcmc.chain)
cat("Chain length: ", chain.length, "\n")

burn.num <- burn.perc* chain.length
cat("Burn in: ", burn.num, "\n")

# check the acceptance rate (ideal is 0.234)
cat("Acceptance Rate (ideal is 0.234) = ", 1 - rejectionRate(mcmc.chain)[1], "\n\n")

# burn in the chain
mcmc.burn <- mcmc(burnAndThin(mcmc.chain, burn=burn.num), start=burn.num)  # TODO chain is already thinned ('s' flag)

# NB. we do not thin the chain because there is no need
# see https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.2041-210X.2011.00131.x

# print the summary stats
print(summary(mcmc.burn))

# calculate the ESS for all params
ess <- effectiveSize(mcmc.burn)
write.table(data.frame(as.list(ess)), file=paste0('data/selection/', prefix, '.ess'), sep='\t')

cat("Effective Sample Size.\n")
print(ess)
cat("\n")

if (min(ess) < 100) {
    cat(paste0("WARNING: ESS below threshold. min(ess) = ", min(ess), '; ', param_file))
}

# plot ESS vs. burn-in
cat("Plotting ESS vs. burn-in.", "\n\n")
pdf(file=paste0('data/pdf/selection/', prefix, '-ess-burn.pdf'), width=21, height=14)
plotESSBurn(mcmc.chain, step.size=round(burn.num/2, 0))
off <- dev.off()

cat("Plotting autocorrelation.\n\n")
pdf(file=paste0('data/pdf/selection/', prefix, '-burn-thin-autocorr.pdf'))
autocorr.plot(mcmc.burn, lag.max=50)
off <- dev.off()

cat("Plotting the trace.\n\n")
pdf(file=paste0('data/pdf/selection/', prefix, '-burn-trace.pdf'))
# png(file=paste0('data/pdf/selection/', prefix, '-burn-trace-pt%d.png'), width=7, height=7, units='in', res=300)
plot(mcmc.burn)
off <- dev.off()
