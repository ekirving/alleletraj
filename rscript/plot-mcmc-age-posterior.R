library(readr)
library(dplyr)
library(ggplot2)

setwd('/Users/Evan/Dropbox/Code/alleletraj')

# load all the param files
files = list.files(path="selection", pattern="*.param", full.names = T)
mcmc.params = suppressMessages(lapply(files, read_tsv) %>% bind_rows())

# convert diffusion units into calendar years
mcmc.params$ageyrs <- mcmc.params$age * 2 * pop_size * gen_time

# TODO command line params
pop_size <- 16000
gen_time <- 8

dom1.age <- -5500
dom2.age <- -4000

max_age <- -100000
brk_width <- 10000

ggplot(mcmc.params, aes(ageyrs)) +

    # show a density plot of the calander ages
    geom_density(adjust=0.5, fill="darkgrey") +

    # plot the ages of the main domestication events
    geom_vline(xintercept=c(dom1.age, dom2.age), linetype = "dashed", colour = 'red') +

    # set the breaks for the x-axis
    scale_x_continuous(limits = c(max_age, 0),
                       breaks = seq(max_age, 0, by = brk_width),
                       labels = seq(max_age, 0, by = brk_width)/1000,
                       minor_breaks = NULL,
                       expand = c(0.01, 0)) +

    # label the plot and the axes
    xlab("kyr BP") +
    ylab("Density") +

    # tweak the theme
    theme_minimal(base_size = 10)
