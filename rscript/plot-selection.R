#!/usr/bin/env Rscript
suppressMessages(source("rscript/path_utilities.r"))

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_prefix <- args[2]
gen_time <- strtoi(args[3])
pop_size <- strtoi(args[4])
burn_in <- strtoi(args[5])

# TODO remove when done testing
# setwd('/Users/Evan/Dropbox/Code/alleletraj')
# input_file <- 'selection/pig-EUD-modsnp_71891.input'
# output_prefix <- 'selection/pig-EUD-modsnp_71891'
# gen_time <- strtoi('5')
# pop_size <- strtoi('20563')
# burn_in <- 5

# load the samples input file
samples <- read.table(input_file, col.names=c('derived_count', 'sample_size', 'bin_high', 'bin_low'))

# get the data in the right format
samples$freq <- samples$derived_count/samples$sample_size
samples$time <- rowMeans(samples[c('bin_high', 'bin_low')]) / (2 * pop_size * gen_time)

# load the MCMC run (WARNING: very CPU and memory costly!)
paths <- read.path(output_prefix)

# plot the trajectory
pdf(file=paste0('pdf/', basename(output_prefix), "-noburnin.pdf"), width = 8, height = 6)
plot.posterior.paths(paths, samples$freq, samples$time)
dev.off()

# plot the trajectoy
pdf(file=paste0('pdf/', basename(output_prefix), ".pdf"), width = 8, height = 6)
plot.posterior.paths(paths, samples$freq, samples$time, burnin=burn_in)
dev.off()
