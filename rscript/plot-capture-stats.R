#!/usr/bin/env Rscript
library(ggplot2)
library(RMySQL)
library(reshape2)

setwd('/Users/Evan/Dropbox/Code/alleletraj')

# connect to the remote server
mydb = dbConnect(MySQL(), user='root', password='', dbname='allele_trajectory', host='localhost')

# execute the query (which caches the resultset)
rs = dbSendQuery(mydb,
    "SELECT sq.accession,
            NULLIF(sq.mapped_q30, 0) AS cap_mapq30,
            NULLIF(sq.duplicates, 0) AS cap_duplic,
            AVG(NULLIF(sq2.mapped_q30, 0)) AS avg_mapq30,
            AVG(NULLIF(sq2.duplicates, 0)) AS avg_duplic
       FROM sample_quality sq
       JOIN sample_quality sq2
         ON sq.accession = sq2.accession
        AND sq2.pool NOT IN ('OXF47', 'OXF57')
      WHERE sq.pool IN ('OXF47', 'OXF57')
   GROUP BY sq.accession"
)

# fetch the resultset from the DB
dat = fetch(rs, n=-1)

# drop any NA data
dat <- na.omit(dat)

# compute the unique endogenous fraction
dat$cap_endo <- dat$cap_mapq30 * (1 - dat$cap_duplic)
dat$avg_endo <- dat$avg_mapq30 * (1 - dat$avg_duplic)

scale_max <- 5
scale_step <- 1

ggplot(dat, aes(x = cap_endo, y = avg_endo)) +
    geom_abline(color='grey') +
    geom_point(size=2, shape=21) +
    geom_smooth(method='lm') +
    ggtitle('Effect of Whole Genome Capture') +
    scale_x_continuous(limits = c(0, scale_max), breaks = seq(0, scale_max, by = scale_step)) +
    scale_y_continuous(limits = c(0, scale_max), breaks = seq(0, scale_max, by = scale_step)) +
    xlab("Unique Endogenous (Capture)") +
    ylab("Unique Endogenous (Pooled)")

dat.melt <- melt(dat[,c('accession', 'cap_endo', 'avg_endo')], id.vars = c("accession"))

ggplot(dat.melt) +
    geom_violin(aes(factor(variable), value)) +
    xlab("Pooled vs. Capture") +
    ylab("Unique Endogenous (log10)") +
    scale_y_continuous(trans='log10', limits=c(1e-3, 110),
                       breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
                       labels = c("0.001%", "0.01%", "0.1%", "1%", "10%", "100%")) +
    scale_x_discrete(limit = c("avg_endo", "cap_endo"), labels = c("Pooled", "Capture"))

# dev.off()
