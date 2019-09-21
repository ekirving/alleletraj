#!/usr/bin/env Rscript
require(R.utils, quietly=T)
source('rscript/mcmc_path_utilities.r')

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_prefix <- args[2]
units <- as.numeric(args[3])
burn_perc <- as.numeric(args[4])
pdf_file <- args[5]

# TODO remove when done testing
# input_file <- 'data/selection/horse-DOM2-modsnp7146976.input'
# output_prefix <- 'data/selection/horse-DOM2-modsnp7146976-n100000-s100-h0.5-chain1'
# units <- 274400
# burn_perc <- 0.5
# pdf_file <- 'data/pdf/selection/horse-DOM2-modsnp7146976-n100000-s100-h0.5-chain1-traj.pdf'

# load the samples input file
samples <- read.table(input_file, col.names=c('derived_count', 'sample_size', 'bin_high', 'bin_low'))

# get the data in the right format
samples$freq <- samples$derived_count/samples$sample_size
samples$time <- rowMeans(samples[c('bin_high', 'bin_low')])
samples$size <- round(log(samples$sample_size +1))/2

# get the chain length
chain.length <- countLines(paste0(output_prefix, '.traj.gz'))[1]

# convert burn % to number of records
burnin <- burn_perc * chain.length

# load the MCMC run (WARNING: very CPU and memory intensive!)
paths <- read.path(output_prefix, chain.length, burnin)

# handle error condition where time and traj are not of equal length
if (length(paths$traj) != length(paths$time)) {
    warning("Time and trajectory files are not of equal length")
    trunc <- min(c(length(paths$traj), length(paths$time)))
    paths$traj <- paths$traj[1:trunc]
    paths$time <- paths$time[1:trunc]
}

# handle truncated final line
if (length(tail(paths$traj, n=1)[[1]]) != length(tail(paths$time, n=1)[[1]])) {
    warning("Final line of time and trajectory files are not of equal length")
    paths$traj <- paths$traj[1:(length(paths$traj)-1)]
    paths$time <- paths$time[1:(length(paths$time)-1)]
}

# plot the trajectoy
pdf(file=pdf_file, width=12, height=7.5)
plot.posterior.paths(paths, samples$freq, samples$time, samples$size, units, xlim=c(-25000/units,0))
dev.off()
