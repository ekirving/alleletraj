library(readr)
library(dplyr)
library(ggplot2)

setwd('/Users/Evan/Dropbox/Code/alleletraj')


# TODO command line params
pop_size <- 16000
gen_time <- 8
burnin <- 1000 # i.e. 20%

# load all the param files
files = list.files(path="selection", pattern="*.param", full.names = T)
mcmc.params = suppressMessages(lapply(files, function(x) read_tsv(x, col_names = F, skip = burnin + 1))) %>% bind_rows()
names(mcmc.params) <- suppressMessages(names(read_tsv(files[1], col_names = T, n_max = 0)))

# convert diffusion units into calendar years
mcmc.params$ageyrs <- mcmc.params$age * 2 * pop_size * gen_time
mcmc.params$firstyrs <- mcmc.params$first_nonzero * 2 * pop_size * gen_time

# convert end_freq back into an actual frequency
mcmc.params$freq <- (1-cos(mcmc.params$end_freq))/2

max_age <- -50000
brk_width <- 5000

ggplot(mcmc.params) +

    # plot 2D density of Age vs Alpha
    stat_density_2d(aes(x=ageyrs, y=alpha1, fill=..level..), geom="polygon",
                    colour="white", show.legend=F) +

    # # plot 2D density of plot Age vs Freq
    # stat_density_2d(aes(x=ageyrs, y=freq, fill=..level..), geom="polygon",
    #                 colour="white", show.legend=F) +

    scale_x_continuous(limits = c(max_age - 5000, 0),
                       breaks = seq(max_age, 0, by = brk_width),
                       labels = seq(max_age, 0, by = brk_width)/1000,
                       minor_breaks = NULL) +

    # using limits() drops all data points that are not within the specified range,
    # causing discontinuity of the density plot, while coord_cartesian() zooms
    # without losing the data points.
    coord_cartesian(xlim = c(max_age, 0)) +

    # label the plot and the axes
    xlab("kyr BP") +
    # ylab("Frequency") +
    ylab("Alpha") +

    # tweak the theme
    theme_minimal(base_size = 10)

