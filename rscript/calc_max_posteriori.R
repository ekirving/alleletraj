library(stringr)
library(coda)

setwd('/Users/Evan/Dropbox/Code/alleletraj')

# TODO command line params
pop_size <- 16000
gen_time <- 8
burnin <- 1000 # i.e. 20%

a_posteriori <- function(filename) {
    # load the data (and apply burnin)
    params <- suppressMessages(read_tsv(filename, col_names = T))
    params <- params[(burnin+1):nrow(params),]

    if (nrow(params) < 4000) {
        return()
    }

    if (min(effectiveSize(params[,!names(params) %in% c("gen")])) < 100) {
        return()
    }

    # convert diffusion units into calendar years
    params$ageyrs <- params$age * 2 * pop_size * gen_time

    # convert end_freq back into an actual frequency
    params$freq <- (1-cos(params$end_freq))/2

    # extract the modsnp id
    id <- as.integer(str_extract(filename, "(?<=-)[0-9]+"))

    # get the max a posteriori values...

    # alpha1
    d = density(na.omit(params$alpha1))
    alpha1 <- d$x[which.max(d$y)]

    # age
    d = density(na.omit(params$ageyrs))
    age <- d$x[which.max(d$y)]

    # end freq
    d = density(na.omit(params$freq))
    end_freq <- d$x[which.max(d$y)]

    data.frame(id, alpha1, age, end_freq)
}

# load all the param files
files = list.files(path="selection", pattern="*.param", full.names = T)
mcmc.post = lapply(files, a_posteriori) %>% bind_rows()
