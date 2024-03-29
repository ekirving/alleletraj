#!/usr/bin/env Rscript
library(RMySQL)
library(ggplot2)
library(ggridges)
library(reshape2)
library(tidyr)
library(RColorBrewer)
library(dplyr)

# get the command line arguments
# args <- commandArgs(trailingOnly = TRUE)
# modsnp_id <- args[1]
# bin_width <- args[2]

# TODO remove when done testing
setwd('/Users/Evan/Dropbox/Code/alleletraj')
modsnp_id <- '25822357'
bin_width <- 500

# connect to the remote server
mydb = dbConnect(MySQL(), user='root', password='', dbname='alleletraj_pig', host='localhost')

# fetch the details of the SNP
rs = dbSendQuery(mydb, paste0(
    "SELECT ms.*
       FROM modern_snps ms
      WHERE ms.id = ", modsnp_id
))

modsnp <- fetch(rs, n=-1)

# get the data and all the reads for this SNP
rs = dbSendQuery(mydb, paste0(
     "SELECT s.accession,
             CASE
                 WHEN COALESCE(s.gmm_status, s.status) LIKE '%wild%'     THEN 'Wild'
                 WHEN COALESCE(s.gmm_status, s.status) LIKE '%domestic%' THEN 'Domestic'
                 WHEN SUBSTR(`group`, 3, 1) = 'W' THEN 'Wild'
                 WHEN SUBSTR(`group`, 3, 1) = 'D' THEN 'Domestic'
                 WHEN haplogroup = 'Y1' THEN 'Domestic'
                 ELSE 'NA'
             END AS `status`,
             sb.bin,
             COALESCE(c14.confident, sd.confident) confident,
             COALESCE(c14.lower, sd.lower, sd.median + 100) lower,
             COALESCE(c14.upper, sd.upper, sd.median - 100) upper,
             GROUP_CONCAT(sr.base ORDER BY sr.base = ms.ancestral DESC SEPARATOR '/') genotype
        FROM modern_snps ms
        JOIN sample_reads sr
          ON sr.chrom = ms.chrom
         AND sr.site = ms.site
         AND sr.called = 1
        JOIN samples s
          ON s.id = sr.sample_id
   LEFT JOIN sample_bins sb
          ON sb.sample_id = s.id
   LEFT JOIN sample_dates sd
          ON s.age = sd.age
   LEFT JOIN sample_dates_c14 c14
          ON c14.accession = s.accession
       WHERE ms.id = ", modsnp_id, "
         AND s.valid = 1
    GROUP BY s.id"
))

# fetch the resultset from the DB
ages <- fetch(rs, n=-1)

# drop any NA data
ages <- na.omit(ages)

# calculate the mean and standard deviation (assume range = μ ± 2σ, i.e. 95%)
ages$mean <- (ages$lower + ages$upper) / 2
ages$sd <- (ages$lower - ages$upper) / 4

# add the genotype the label name
ages$label <- paste0(ages$accession, " [", substr(paste0(ages$genotype, "/."), 1, 3), "]")

# calculate the median age of each bin
ages$median <- as.integer(gsub( " .*$", "", ages$bin)) - bin_width/2

# find the max value of the oldest bin
max_lower <- ceiling(max(ages$lower)/bin_width)*bin_width

# sort the data by lower bound
ages <- ages[order(-ages$mean,-ages$sd),]

# set the ordering for things
statuses <- c('Wild','Domestic','NA')
alleles <- c(modsnp[['ancestral']], modsnp[['derived']])

# simulate 100 variable for each distribution
sim.data <- data.frame(apply(ages[,c('mean','sd')], 1, function(x) rnorm(10000, mean=x[1], sd=x[2])))
colnames(sim.data) <- ages$label

# melt the data
sim.melt <- melt(sim.data)
colnames(sim.melt) <- c('label', 'value')

# merge the status back in
sim.melt <- merge(sim.melt, ages[,c('label', 'status')], by = "label")

# convert diploids into haploids
ages.grp <- separate_rows(ages[c('status', 'median', 'genotype')], c('genotype'), sep='/')

# get the count of alleles in each group
ages.grp <- ages.grp %>%
    group_by(status, median, genotype) %>%
    summarise(count = n())

# setup the title of the plot
title <- paste0("SNP chr", modsnp['chrom'], ":", modsnp['site'],
                " [", modsnp['ancestral'], "/", modsnp['derived'], "]",
                " DAF=", round(modsnp['daf'], 2),
                " (#", modsnp_id, ")")

# show all available colour palettes
# display.brewer.all()

# bind the fill elements to colours
colour_map <- setNames(brewer.pal(5,"Set1")[c(2,3,1,4,5)],
                       c(statuses, alleles))

# TODO fixe me!
colour_map <- c(colour_map, c("D-A" = "#b2df8a",
                              "D-T" = "#33a02c",
                              "N-T" = "#e31a1c",
                              "W-A" = "#a6cee3",
                              "W-T" = "#1f78b4"))

# TODO fixe me!
alleles_long <- c("D-T", "D-A", "W-T", "W-A")

ages.grp$statgeno <- paste0(substr(ages.grp$status, 1, 1), '-', ages.grp$genotype)

# pdf(file=paste('rscript/', basename, '.pdf', sep=''), width = 10, height = length(ages$accession)/7)

ggplot() +

    # display the Wild alelle counts as a stacked column chart
    geom_col(data=ages.grp[ages.grp$status == 'Wild',],
             aes(x=median - bin_width/4, y=count, fill=statgeno, colour=statgeno),
             width=bin_width/2) +

    # display the Domestic alelle counts as a stacked column chart
    geom_col(data=ages.grp[ages.grp$status == 'Domestic',],
             aes(x=median + bin_width/4, y=count, fill=statgeno, colour=statgeno),
             width=bin_width/2) +

    # display the sample dates as a rigline plot
    geom_density_ridges(data=sim.melt,
                        aes(x = value, y = label, fill=status, color=status),
                        scale = 3) +

    # set a manual scale that only shows the status elements
    scale_fill_manual(name="Status",
                      breaks=statuses,
                      values=colour_map,
                      guide=guide_legend(override.aes = list(color=colour_map[statuses]))) +

    # setup a dummy scale to show the genotypes (we need to override the aesthetic)
    scale_color_manual(name="Genotype",
                       breaks=alleles_long,
                       values=rep("lightgrey", 7),
                       labels=c('Ancestral (T)', 'Derived (A)','Ancestral (T)', 'Derived (A)'),
                       guide=guide_legend(override.aes = list(fill=colour_map[alleles_long],
                                                              color=colour_map[alleles_long]))) +

    # set the breaks for the x-axis
    scale_x_continuous(breaks = seq(0, max_lower, by = bin_width),
                       labels = seq(0, max_lower, by = bin_width)/1000,
                       minor_breaks = NULL,
                       expand = c(0.01, 0)) +

    # add a little padding to the y-axis
    scale_y_discrete(expand = c(0.01, 1)) +

    # label the plot and the axes
    ggtitle(title) +
    xlab("kyr BP") +
    ylab("Samples") +

    # tweak the theme
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_text(family="Courier", vjust = 0),
          panel.grid.major.x=element_line(colour="lightgrey"))

# dev.off()
