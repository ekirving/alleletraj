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
mcmc.params$firstyrs <- mcmc.params$first_nonzero * 2 * pop_size * gen_time

dom1.age <- -5500
dom2.age <- -4000

ggplot(mcmc.params) +

    # show a histogram of the first non-zero derived allele
    geom_histogram(aes(firstyrs), fill="#48a4e2", colour="#0d4869", binwidth = 500) +

    # plot the ages of the main domestication events
    geom_vline(xintercept=c(dom1.age, dom2.age), linetype = "dashed", colour = 'red') +

    # label the plot and the axes
    xlab("years BP") +
    ylab("Count") +

    # tweak the theme
    theme_minimal(base_size = 10) +
    theme(
        # remove the vertical grid lines
        panel.grid.major.x = element_blank()
    )
