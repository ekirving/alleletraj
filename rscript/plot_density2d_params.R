#!/usr/bin/env Rscript

library(readr, quietly = T)
library(dplyr, quietly = T)
library(ggplot2, quietly = T)
library(viridis, quietly = T)
library(stringr, quietly = T)
library(coda, quietly = T)

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
species <- args[1]
population <- args[2]
burnin <- strtoi(args[3])
pop_size <- strtoi(args[4])
gen_time <- strtoi(args[5])

# TODO remove when done testing
# setwd('/Users/Evan/Dropbox/Code/alleletraj')
# species <- 'horse'
# population <- 'DOM2'
# burnin <- 1000 # i.e. 20%
# pop_size <- 16000
# gen_time <- 8

# make a helper function to load the param files
read_tsv_param <- function(filename) {

    params <- suppressMessages(read_tsv(filename, col_names = T))
    params <- params[(burnin+1):nrow(params),]

    if (nrow(params) < 4000) {
        return()
    }

    # enforce a min ESS
    if (min(effectiveSize(params[,!names(params) %in% c("gen")])) < 100) {
        return()
    }

    # extract the modsnp id
    params$id <- as.integer(str_extract(filename, "(?<=-)[0-9]+"))

    params
}

# load all the param files
files <- list.files(path = "selection", pattern = paste(species, population, "(.+).param", sep = '-'), full.names = T)
mcmc.params <- lapply(files, read_tsv_param) %>% bind_rows()

# convert diffusion units into calendar years
mcmc.params$ageyrs <- mcmc.params$age * 2 * pop_size * gen_time

# function for plotting the MCMC params
plot_density2d <- function(param, ylab) {

    # set default limits, breaks and labels
    min_x <- -50000
    max_x <- 0
    brk_w <- 5000
    lim_x <- 150000

    # setup shared layers for the density2d plot
    gglayers <- list(

        # set the colour scheme for the density plot
        scale_fill_viridis(name = "Density", direction = -1),

        # set the breaks for the x-axis
        scale_x_continuous(limits = c(min_x - lim_x, max_x),
                           breaks = seq(min_x, max_x, by = brk_w),
                           labels = seq(min_x, max_x, by = brk_w) / 1000,
                           minor_breaks = NULL),

        # using limits() drops all data points that are not within the specified range,
        # causing discontinuity of the density plot, while coord_cartesian() zooms
        # without losing the data points.
        coord_cartesian(xlim = c(min_x, max_x)),

        # label the plot and the axes
        xlab("kyr BP"),

        # tweak the theme
        theme_minimal(base_size = 10),

        # remove the vertical grid lines
        theme(panel.grid.major.x = element_blank())
    )

    pdf(file = paste0('rscript/pdf/', species, '-', population ,'-density2d-age-', param, '.pdf'),
        width = 16, height = 9)

    g <- ggplot(mcmc.params) +

        # plot density2d
        stat_density_2d(
            aes_string(x = 'ageyrs', y = param, fill = '..level..'),
            geom = "polygon", show.legend = F, bins = 20
        ) +

        # label the y-axis
        ylab(ylab) +

        # add the common layers
        gglayers

    # draw the finished plot
    print(g)

    dev.off()
}

# ------------------------------------------------------------------------------
# plot the End frequency
# ------------------------------------------------------------------------------

# convert end_freq back into an actual frequency
mcmc.params$freq <- (1 - cos(mcmc.params$end_freq)) / 2

plot_density2d(param = 'freq', ylab = 'End Frequency')

# ------------------------------------------------------------------------------
# plot the Selection coefficients
# ------------------------------------------------------------------------------

# convert alpha params into s values
mcmc.params$s1 <- mcmc.params$alpha1 / (2 * pop_size)
mcmc.params$s2 <- mcmc.params$alpha2 / (2 * pop_size)

plot_density2d(param = 's1', ylab = expression(paste("s"[1])))
plot_density2d(param = 's2', ylab = expression(paste("s"[2])))
