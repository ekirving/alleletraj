#!/usr/bin/env Rscript
library(ggplot2)
library(reshape2)

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
stub=args[1]
label=args[2]

# TODO remove when done testing
# setwd('/Users/Evan/Dropbox/Code/alleletraj')
# stub = 'all-snps-ages'
# label = 'Derived Allele vs. Sample Age'

# read the data file
dat <- read.table(paste0("tsv/", stub, ".tsv"), header=TRUE, sep ='\t') #, stringsAsFactors=F)

# drop any NA data
dat <- na.omit(dat)

# transversions only (i.e. no transitions [A <-> G] and [C <-> T]
# dat <- dat[!(dat$ancestral %in% c('A', 'G') & dat$derived %in% c('A', 'G')),]
# dat <- dat[!(dat$ancestral %in% c('C', 'T') & dat$derived %in% c('C', 'T')),]

pdf(file=paste0("pdf/", stub, ".pdf"), width = 8, height = 7)

ggplot(dat, aes(x = oldest_derived, y = oldest_sample, fill=trait)) +
    geom_abline(color='grey') +
    geom_point(size=2, shape=21) +
    theme(
        plot.title = element_text(hjust = 0.5),
        legend.key = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        panel.background = element_blank()
    ) +
    ggtitle(label) +
    scale_fill_discrete(name="QTL Trait Class") +
    scale_x_continuous(limits = c(0, 13000), breaks = seq(0, 13000, by = 1000)) +
    scale_y_continuous(limits = c(0, 13000), breaks = seq(0, 13000, by = 1000)) +
    xlab("Earliest derived allele (BP)") +
    ylab("Earliest sample (BP)")

dev.off()
