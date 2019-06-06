#!/usr/bin/env Rscript
source("rscript/path_utilities.r")

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_prefix <- args[2]
units <- as.numeric(args[3])
burnin <- as.numeric(args[4])
thin <- as.numeric(args[5])
pdf_file <- args[6]

# TODO remove when done testing
# input_file <- 'data/selection/horse-DOM2-modsnp9016431.input'
# output_prefix <- 'data/selection/horse-DOM2-modsnp9016431-n50000-s100-h0.5-chain1'
# units <- 2 * 54787 * 8
# burnin <- 5
# thin <- 5
# pdf_file <- 'data/pdf/selection/horse-DOM2-modsnp9016431-n50000-s100-h0.5-chain1-traj.pdf'

# load the samples input file
samples <- read.table(input_file, col.names=c('derived_count', 'sample_size', 'bin_high', 'bin_low'))

# get the data in the right format
samples$freq <- samples$derived_count/samples$sample_size
samples$time <- rowMeans(samples[c('bin_high', 'bin_low')])

# load the MCMC run (WARNING: very CPU and memory intensive!)
paths <- read.path(output_prefix)

# TODO rescale time by units

# plot the trajectoy
pdf(file=pdf_file, width = 8, height = 6)
plot.posterior.paths(paths, samples$freq, samples$time, burnin=burnin/thin)
dev.off()

