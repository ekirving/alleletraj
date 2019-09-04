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
thin <- as.numeric(args[3])
diff_units <- as.numeric(args[4])
pop_size <- as.numeric(args[5])
ess_file <- args[6]
map_file <- args[7]
ess_pdf <- args[8]
auto_pdf <- args[9]
trace_png <- args[10]

# TODO remove when done testing
# prefix <- 'horse-DOM2-modsnp3020195-n100000000-s100-h0.5-chain1'
# param_file <- paste0('data/selection/', prefix, '.param.gz')
# burn_perc <- 0.5
# thin <- 1000
# diff_units <- 274400
# pop_size <- 17150
# ess_file <- paste0('data/selection/', prefix, '.ess')
# map_file <- paste0('data/selection/', prefix, '.map')
# ess_pdf <- paste0('data/pdf/selection/', prefix, '-ess-burn.pdf')
# auto_pdf <- paste0('data/pdf/selection/', prefix, '-burn-thin-autocorr.pdf')
# trace_png <- paste0('data/pdf/selection/', prefix, '-burn-trace-pt1.png')

cat("Loading MCMC...\n")
cat("Param: ", param_file, "\n")

# load the chain and convert to an MCMC object (and cast infinite values to NA)
chain <- fread(param_file, header = T, sep = '\t', drop=c('gen'), na.strings=c('inf', '-inf'))

# drop NAs
chain <- na.omit(chain)

# convert to MCMC object
mcmc.chain <- mcmc(chain)
chain.length <- nrow(mcmc.chain)

# convert burn % to number of records
burnin = burn_perc * chain.length

cat("Chain length: ", chain.length * thin, "\n")
cat("Burn in: ", burnin * thin, "\n")
cat("Thin: ", thin, "\n")

# check the acceptance rate / except this doesn't work for pre-thinnned chains
# cat("Acceptance Rate (ideal is 0.234) = ", 1 - rejectionRate(mcmc.chain)[1], "\n\n")

# drop all the sample_time and first_nonzero params as they dominate the ESS space
mcmc.notime <- mcmc(chain[,c('lnL','pathlnL','alpha1','alpha2','F','age','end_freq')])

# plot ESS vs. burn-in
cat("Plotting ESS vs. burn-in.", "\n\n")
pdf(file=ess_pdf)  # , width=21, height=14
plotESSBurn(mcmc.notime, longest.burn.in=chain.length * 0.5, step.size=chain.length * 0.1)
off <- dev.off()

# burn in the chain (thinning is already done)
mcmc.burn <- mcmc(
    burnAndThin(mcmc.chain, burn = burnin),
    start = burnin * thin,
    thin = thin)

# print the summary stats
print(summary(mcmc.burn))

# calculate the ESS for all params
ess <- effectiveSize(mcmc.burn)
cat(toJSON(ess, indent = 2), file=ess_file)

cat("Effective Sample Size.\n")
print(ess)
cat("\n")

if (min(ess) < 100) {
    cat(paste0("WARNING: ESS below threshold. min(ess) = ", min(ess), "\n\n"))
}

# calculate the maximum a posteriori of all the parameters
map <- c()
for (param in varnames(mcmc.burn)) {
    d = density(na.omit(mcmc.burn[,param]))
    map[[param]] <- d$x[which.max(d$y)]

    if (any(str_detect(param, c('age', 'time', 'first')))) {
        # convert diffusion units into calendar years
        map[[param]] <- -round(map[[param]] * diff_units, 0)

    } else if (param =='end_freq') {
        # convert end_freq back into an actual frequency
        map[[param]] <- (1 - cos(map[[param]])) / 2

    } else if (str_detect(param, 'alpha')) {
        # convert alpha params into s values
        map[[param]] <- map[[param]] / (2 * pop_size)
    }
}
cat(toJSON(map, indent = 2), file=map_file)

cat("Maximum a posteriori.\n")
print(map)
cat("\n")

cat("Plotting autocorrelation.\n\n")
pdf(file=auto_pdf)
autocorr.plot(mcmc.burn, lag.max=50)
off <- dev.off()

cat("Plotting the trace.\n\n")
png(file=str_replace(trace_png, 'pt1', 'pt%d'), width=7, height=7, units='in', res=300)
plot(mcmc.burn)
off <- dev.off()

cat("FINISHED\n")
