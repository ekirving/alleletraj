#!/usr/bin/env Rscript
suppressWarnings(source("rscript/path_utilities.r"))

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_prefix <- args[2]
pop_size <- strtoi(args[3])
burn_in <- strtoi(args[3])

# TODO remove when done testing
# setwd('/Users/Evan/Dropbox/Code/alleletraj')
# input_file <- 'selection/pig-EUD-modsnp_25822357.input'
# output_prefix <- 'selection/pig-EUD-modsnp_25822357'
# pop_size <- strtoi('20563')
# burn_in <- 500

# load the samples input file
samples <- read.table(input_file, col.names=c('derived_count', 'sample_size', 'bin_high', 'bin_low'))

# get the data in the right format
samples$freq <- samples$derived_count/samples$sample_size
samples$time <- rowMeans(samples[c('bin_high', 'bin_low')], na.rm=TRUE) / pop_size

# load the MCMC run (very slow and memory costly)
paths <- read.path(output_prefix)

pdf(file=paste0(output_prefix, ".pdf"), width = 8, height = 6)

# plot the trajectory
plot.posterior.paths(paths, samples$freq, samples$time, burnin=burn_in)

dev.off()
