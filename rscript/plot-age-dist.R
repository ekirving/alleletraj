# install.packages("ggridges")

library(ggplot2)
library(ggridges)
library(reshape2)

setwd('/Users/Evan/Dropbox/Code/alleletraj')

# basename <- 'dates-all-samples'
# basename <- 'dates-modsnp-547444'
# basename <- 'dates-modsnp-4499449'
# basename <- 'dates-modsnp-8130030'
# basename <- 'dates-modsnp-8130095'
# basename <- 'dates-modsnp-8130326'
basename <- 'dates-modsnp-5242077'

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
