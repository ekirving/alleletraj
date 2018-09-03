library(readr)
library(dplyr)
library(ggplot2)
library(viridis)
library(coda)

setwd('/Users/Evan/Dropbox/Code/alleletraj')

# TODO command line params
pop_size <- 16000
gen_time <- 8
burnin <- 1000 # i.e. 20%

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
files = list.files(path="selection", pattern="*.param", full.names = T)
mcmc.params = suppressMessages(lapply(files, read_tsv_param)) %>% bind_rows()

# convert diffusion units into calendar years
mcmc.params$ageyrs <- mcmc.params$age * 2 * pop_size * gen_time
mcmc.params$firstyrs <- mcmc.params$first_nonzero * 2 * pop_size * gen_time

# convert end_freq back into an actual frequency
mcmc.params$freq <- (1-cos(mcmc.params$end_freq))/2

max_age <- -50000
brk_width <- 5000

# ------------------------------------------------------------------------------
# --                           Age vs Alpha                                   --
# ------------------------------------------------------------------------------

pdf(file=paste('rscript/mcmc-2d-age-alpha.pdf', sep=''), width = 8, height = 4.5)

ggplot(mcmc.params) +

    # plot 2D density of Age vs Alpha
    stat_density_2d(aes(x=ageyrs, y=alpha1, fill=..level..), geom="polygon", show.legend=F) +

    scale_fill_viridis(name = "Density", direction = -1) +

    scale_x_continuous(limits = c(max_age - 150000, 0),
                       breaks = seq(max_age, 0, by = brk_width),
                       labels = seq(max_age, 0, by = brk_width)/1000,
                       minor_breaks = NULL) +

    # using limits() drops all data points that are not within the specified range,
    # causing discontinuity of the density plot, while coord_cartesian() zooms
    # without losing the data points.
    coord_cartesian(xlim = c(max_age, 0)) +

    # label the plot and the axes
    xlab("kyr BP") +
    ylab("Alpha") +

    # tweak the theme
    theme_minimal(base_size = 10)

dev.off()

# ------------------------------------------------------------------------------
# --                           Age vs End Freq                                --
# ------------------------------------------------------------------------------

pdf(file=paste('rscript/mcmc-2d-age-end_freq.pdf', sep=''), width = 8, height = 4.5)

ggplot(mcmc.params) +

    # plot 2D density of plot Age vs Freq
    stat_density_2d(aes(x=ageyrs, y=freq, fill=..level..), geom="polygon", show.legend=F) +

    scale_fill_viridis(name = "Density", direction = -1) +

    scale_x_continuous(limits = c(max_age - 150000, 0),
                       breaks = seq(max_age, 0, by = brk_width),
                       labels = seq(max_age, 0, by = brk_width)/1000,
                       minor_breaks = NULL) +

    # using limits() drops all data points that are not within the specified range,
    # causing discontinuity of the density plot, while coord_cartesian() zooms
    # without losing the data points.
    coord_cartesian(xlim = c(max_age, 0)) +

    # label the plot and the axes
    xlab("kyr BP") +
    ylab("Frequency") +

    # tweak the theme
    theme_minimal(base_size = 10)

dev.off()

# ------------------------------------------------------------------------------
# --                           Age vs End Freq                                --
# ------------------------------------------------------------------------------

pdf(file=paste('rscript/mcmc-2d-alpha-end_freq.pdf', sep=''), width = 8, height = 4.5)

ggplot(mcmc.params) +

    # plot 2D density of plot Age vs Freq
    stat_density_2d(aes(x=alpha1, y=freq, fill=..level..), geom="polygon", show.legend=F) +

    scale_fill_viridis(name = "Density", direction = -1) +

    # label the plot and the axes
    xlab("Alpha1") +
    ylab("Frequency") +

    # tweak the theme
    theme_minimal(base_size = 10)

dev.off()
