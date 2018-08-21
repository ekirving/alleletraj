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

min_x <- -50
max_x <- 250
brk_width <- 50

ggplot(mcmc.params, aes(alpha1)) +

    # show a density plot of the alpha1 selection coefficient
    geom_density(fill="#48a4e2", colour="#0d4869") +

    # set the breaks for the x-axis
    scale_x_continuous(limits = c(min_x-50, max_x+50),
                       breaks = seq(min_x, max_x, by = brk_width),
                       labels = seq(min_x, max_x, by = brk_width),
                       minor_breaks = NULL,
                       expand = c(0.01, 0)) +

    # using limits() drops all data points that are not within the specified range,
    # causing discontinuity of the density plot, while coord_cartesian() zooms
    # without losing the data points.
    coord_cartesian(xlim = c(min_x, max_x)) +

    # label the plot and the axes
    xlab("Alpha1") +
    ylab("Density") +

    # tweak the theme
    theme_minimal(base_size = 10) +
    theme(
        # remove the vertical grid lines
        panel.grid.major.x = element_blank()
    )
