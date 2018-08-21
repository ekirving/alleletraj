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

# convert end_freq back into an actual frequency
mcmc.params$freq <- (1-cos(mcmc.params$end_freq))/2

# TODO command line params
pop_size <- 16000
gen_time <- 8

max_age <- -50000
brk_width <- 5000

ggplot(mcmc.params) +

    # show a histogram
    # geom_density_2d(aes(ageyrs, alpha1)) +
    geom_density_2d(aes(ageyrs, freq)) +

    scale_x_continuous(limits = c(max_age - 5000, 0),
                       breaks = seq(max_age, 0, by = brk_width),
                       labels = seq(max_age, 0, by = brk_width)/1000,
                       minor_breaks = NULL,
                       expand = c(0.01, 0))


ggplot(mcmc.params) +
    geom_density_2d(aes(alpha1, freq))


ggplot(mcmc.params) +
    geom_contour(aes(alpha1, freq), stat="density2d", bins = 30)
