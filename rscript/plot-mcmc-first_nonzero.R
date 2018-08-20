library(readr)
library(dplyr)
library(ggplot2)

setwd('/Users/Evan/Dropbox/Code/alleletraj')

burnin <- 1000 # i.e. 20%

# load all the param files
files = list.files(path="selection", pattern="*.param", full.names = T)
mcmc.params = suppressMessages(lapply(files, function(x) read_tsv(x, col_names = F, skip = burnin + 1))) %>% bind_rows()
names(mcmc.params) <- suppressMessages(names(read_tsv(files[1], col_names = T, n_max = 0)))

# convert diffusion units into calendar years
mcmc.params$ageyrs <- mcmc.params$age * 2 * pop_size * gen_time
mcmc.params$firstyrs <- mcmc.params$first_nonzero * 2 * pop_size * gen_time

# TODO command line params
pop_size <- 16000
gen_time <- 8

dom1.age <- -5500
dom2.age <- -4000

max_age <- -50000
brk_width <- 5000

ggplot(mcmc.params, aes(firstyrs)) +

    # show a histogram
    geom_histogram(fill="darkgrey") +

    # plot the ages of the main domestication events
    geom_vline(xintercept=c(dom1.age, dom2.age), linetype = "dashed", colour = 'red') +

    # label the plot and the axes
    xlab("years BP") +
    ylab("Density") +

    # tweak the theme
    theme_minimal(base_size = 10)
