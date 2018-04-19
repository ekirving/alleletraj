#!/usr/bin/env Rscript
library(ggplot2)
library(RMySQL)
library(reshape2)

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
modsnp_id <- args[1]

# TODO remove when done testing
setwd('/Users/Evan/Dropbox/Code/alleletraj')
modsnp_id <- '8130095'  # 51 reads

# connect to the remote server
mydb = dbConnect(MySQL(), user='root', password='', dbname='allele_trajectory', host='localhost')

# execute the query (which caches the resultset)
rs = dbSendQuery(mydb,
     "SELECT s.id, s.accession,
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
             GROUP_CONCAT(sr.base ORDER BY sr.base) alleles
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
       WHERE ms.id = 8130095
         AND s.valid = 1
    GROUP BY s.id
      HAVING confident IS NOT NULL"
)

# fetch the resultset from the DB
dat = fetch(rs, n=-1)

# import the sample ages
ages <- read.table(paste('rscript/', basename, '.tsv', sep=''), header = T)

# drop missing dates
ages <- ages[!is.na(ages$confident),]

# calculate the mean and standard deviation (assume range = μ ± 2σ, i.e. 95%)
ages$mean <- (ages$lower + ages$upper) / 2
ages$sd <- (ages$lower - ages$upper) / 4

# add the number of SNPs to the label name
ages$label <- paste(ages$accession, " (#", ages$num_snp, " SNPs)", sep="")

# sort the data by lower bound
ages <- ages[order(-ages$mean,-ages$sd),]

# simulate 100 variable for each distribution
sim.data <- data.frame(apply(ages[,c('mean','sd')], 1, function(x) rnorm(10000, mean=x[1], sd=x[2])))
colnames(sim.data) <- ages$label

# melt the data
sim.melt <- melt(sim.data)

pdf(file=paste('rscript/', basename, '.pdf', sep=''), width = 10, height = length(ages$accession)/7)

# display it as a rigline plot
ggplot(sim.melt, aes(x = value, y = variable, fill = variable, color=variable)) +
    geom_density_ridges(scale = 3) +
    theme_minimal(base_size = 10) + theme(axis.text.y = element_text(vjust = 0)) +
    scale_x_continuous(expand = c(0.01, 0)) +
    scale_y_discrete(expand = c(0.01, 0)) +
    scale_fill_cyclical(values = c("#1f78b4", "#33a02c")) +
    scale_color_cyclical(values = c("#D3D3D3"))

dev.off()
