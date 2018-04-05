#!/usr/bin/env Rscript
library(ggplot2)
library(reshape2)

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
stub=args[1]
label=args[2]

# TODO remove when done testing
# setwd('/Users/Evan/Dropbox/Code/alleletraj')
# stub = 'sample_bin-05500-05001'
# label = '5001 - 5500 BP'

# read the data file
dat = read.table(paste0("tsv/", stub, ".tsv"), header=TRUE, sep ='\t') #, stringsAsFactors=F)

# add 'Unknown' to the factor levels
levels(dat$status) <- c(levels(dat$status), 'Unknown')

# set NA to 'Unknown'
dat$status[is.na(dat$status)] <- 'Unknown'

# order by map_prcnt DESC
dat <- dat[order(-rank(dat$map_prcnt)),]

# perserve the ordering by resetting the levels
dat$accession <- factor(dat$accession, levels = dat$accession)

# pdf(file=paste0("pdf/", stub, ".pdf"), width = length(dat$accession)/3, height = max(dat$map_prcnt)/3)
pdf(file=paste0("pdf/", stub, ".pdf"), width = 7, height = 5)

ggplot(dat, aes(x = accession, y = map_prcnt, fill=status)) +
    geom_bar(stat='identity', width=1) +
    theme(
        plot.title = element_text(hjust = 0.5),
          legend.key = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5),
          panel.background = element_blank()
          ) +
    ggtitle(label) +
    scale_y_continuous(breaks = seq(0, 100, by = 1)) +
    xlab("Accession") +
    ylab("Map %")

dev.off()
