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

# load the sample median ages
age <- read.table("tsv/all-snps-counts.tsv", header=TRUE, sep ='\t')

# count the number of SNPs covered in 100 year bins
# age.hist <- hist(age$median, breaks = seq(0, 13000, by=100), plot=FALSE)

# log10 transform the result, as high coverage samples have millions of SNPs
# age.log10 <- data.frame(age.hist$mids, log10(age.hist$counts)*500)
# colnames(age.log10) <- c('mid', 'log10cnt')

pdf(file=paste0("pdf/", stub, ".pdf"), width = 8, height = 7)

ggplot(dat) +
    geom_abline(color='grey') +
    geom_histogram(data=age,binwidth = 100, aes(x = median, weight = 1/200), alpha = 0.3) +
    # geom_bar(stat="identity", data=age.log10, aes(x = mid, y = log10cnt)) +
    geom_point(size=2, shape=21, aes(x = oldest_derived, y = oldest_sample, fill=trait)) +
    theme(
        plot.title = element_text(hjust = 0.5),
        legend.key = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        panel.background = element_blank()
    ) +
    ggtitle(label) +
    scale_fill_discrete(name="QTL Trait Class") +
    scale_x_continuous(limits = c(0, 13150), breaks = seq(0, 13000, by = 1000)) +
    scale_y_continuous(limits = c(0, 13150), breaks = seq(0, 13000, by = 1000),
                       sec.axis = sec_axis(trans=~.*200, name = "SNP count",
                                           breaks = seq(0, 3e6, by = 5e5))) +
    xlab("Earliest derived allele (BP)") +
    ylab("Earliest sample (BP)")

dev.off()
