#!/usr/bin/env Rscript
library(RMySQL)
library(ggplot2)
library(ggridges)
library(reshape2)
library(tidyr)
library(RColorBrewer)

# get the command line arguments
# args <- commandArgs(trailingOnly = TRUE)
# modsnp_id <- args[1]
# bin_width <- args[2]

# TODO remove when done testing
setwd('/Users/Evan/Dropbox/Code/alleletraj')
modsnp_id <- '8130095'
bin_width <- 500

# connect to the remote server
mydb = dbConnect(MySQL(), user='root', password='', dbname='allele_trajectory', host='localhost')

# execute the query (which caches the resultset)
rs = dbSendQuery(mydb, paste0(
    "SELECT ms.*
       FROM modern_snps ms
      WHERE ms.id = ", modsnp_id
))

modsnp <- fetch(rs, n=-1)

# execute the query (which caches the resultset)
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
         AND sr.snp = 1
        JOIN samples s
          ON s.id = sr.sample_id
         AND s.species = 'pig'
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

# sort the data by lower bound
ages <- ages[order(-ages$mean,-ages$sd),]

# set the order of the status factors
ages$status <- factor(ages$status, c('Wild','Domestic','NA'))

# simulate 100 variable for each distribution
sim.data <- data.frame(apply(ages[,c('mean','sd')], 1, function(x) rnorm(10000, mean=x[1], sd=x[2])))
colnames(sim.data) <- ages$label

# melt the data
sim.melt <- melt(sim.data)
colnames(sim.melt) <- c('label', 'value')

# merge the status back in
sim.melt <- merge(sim.melt, ages[,c('label', 'status')], by = "label")

# convert diploids into haploids
alleles <- separate_rows(ages[c('status', 'bin', 'genotype')], c('genotype'), sep='/')

# set the order of the genotype factors
alleles$genotype <- factor(alleles$genotype, c('T','A'))  # TODO fix this

# calculate the median age of each bin
alleles$median <- as.integer(gsub( " .*$", "", alleles$bin)) - bin_width/2

# get the count of alleles in each group
alleles.grp <- alleles %>%
    group_by(status, median, genotype) %>%
    summarise(count = n())

# show all available colour palettes
# display.brewer.all()

# setup the title of the plot
title <- paste0("SNP chr", modsnp['chrom'], ":", modsnp['site'],
                " [", modsnp['ancestral'], "/", modsnp['derived'], "]",
                " maf=", round(modsnp['maf'], 2),
                " (#", modsnp_id, ")")

# bind the fill elements to colours
colour_map <- setNames(brewer.pal(5,"Set1")[c(2,3,1,4,5)],
                       c(levels(sim.melt$status), levels(alleles$genotype)))

# pdf(file=paste('rscript/', basename, '.pdf', sep=''), width = 10, height = length(ages$accession)/7)

ggplot() +

    # display the Domestic alelle counts as a stacked column chart
    geom_col(data=alleles.grp[alleles.grp$status == 'Domestic',],
         aes(x=median - bin_width/4, y=count, fill=genotype, colour=genotype),
         width=bin_width/2) +

    # display the Wild alelle counts as a stacked column chart
    geom_col(data=alleles.grp[alleles.grp$status == 'Wild',],
             aes(x=median + bin_width/4, y=count, fill=genotype, colour=genotype),
             width=bin_width/2) +

    # display the sample dates as a rigline plot
    geom_density_ridges(data=sim.melt,
                        aes(x = value, y = label, fill=status, color=status),
                        scale = 3) +

    # set a manual scale that only shows the status elements
    scale_fill_manual(name="Status",
                      breaks=levels(sim.melt$status),
                      values=colour_map,
                      guide=guide_legend(override.aes = list(color=colour_map[levels(sim.melt$status)]))) +

    # setup a dummy scale to show the genotypes (we need to override the aesthetic)
    scale_color_manual(name="Genotype",
                       breaks=levels(alleles$genotype),
                       values=rep("lightgrey", 5),
                       guide=guide_legend(override.aes = list(fill=colour_map[levels(alleles$genotype)],
                                                              color=colour_map[levels(alleles$genotype)]))) +

    # set the breaks for the x-axis
    scale_x_continuous(breaks = seq(0, 11500, by = bin_width),
                       labels = seq(0, 11500, by = bin_width)/1000,
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
          # set the colour of the grid lines
          panel.grid.major.x=element_line(colour="lightgrey"))

# dev.off()
