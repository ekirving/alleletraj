library(readr)
library(dplyr)
library(ggplot2)
library(ggridges)
library(RMySQL)
library(stringr)
library(viridis)
library(stringr)

setwd('/Users/Evan/Dropbox/Code/alleletraj')

# TODO command line params
pop_size <- 16000
gen_time <- 8
burnin <- 1000 # i.e. 20%

# make a helper function to load the param files
read_tsv_param <- function(filename) {
    ret <- suppressMessages(read_tsv(filename, col_names = F, skip = burnin + 1))
    if (nrow(ret) > 1) {
        # extract the modsnp id
        ret$source <- as.integer(str_extract(filename, "(?<=-)[0-9]+"))
    }
    ret
}

# load all the param files
files = list.files(path="selection", pattern="*.param", full.names = T)
mcmc.params = lapply(files, read_tsv_param) %>% bind_rows()
names(mcmc.params) <- suppressMessages(c(names(read_tsv(files[1], col_names = T, n_max = 0)), 'id'))

# convert diffusion units into calendar years
mcmc.params$ageyrs <- mcmc.params$age * 2 * pop_size * gen_time

# convert end_freq back into an actual frequency
mcmc.params$freq <- (1-cos(mcmc.params$end_freq))/2

# connect to the remote server
mydb = dbConnect(MySQL(), user='root', password='', dbname='alleletraj_horse', host='localhost')

# fetch the details of the SNP
rs = dbSendQuery(mydb,
    "SELECT DISTINCT ms.id, t.class AS trait
       FROM qtls q
       JOIN traits t
         ON t.id = q.trait_id
       JOIN ensembl_variants ev
         ON ev.rsnumber = q.peak
       JOIN modern_snps ms
         ON ms.variant_id = ev.id
       JOIN qtl_snps qs
         ON qs.qtl_id = q.id
        AND qs.modsnp_id = ms.id
      WHERE q.associationType = 'Association'
        AND q.valid = 1")

modsnp <- fetch(rs, n=-1)

# join the traits onto the param files
mcmc.params <- left_join(mcmc.params, modsnp, by = 'id')

dom1.age <- -5500
dom2.age <- -4000

max_age <- -50000
brk_width <- 5000

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

pdf(file=paste('rscript/mcmc-ridgeline-age.pdf', sep=''), width = 8, height = 4.5)

ggplot() +

    # display the sample dates as a rigline plot
    stat_density_ridges(data=mcmc.params,
                        aes(x=ageyrs, y=trait, fill=0.5 - abs(0.5-..ecdf..)),
                        scale = 3,
                        bandwidth=1000,
                        geom = "density_ridges_gradient", calc_ecdf = TRUE) +

    scale_fill_viridis(name = "Posterior", direction = -1) +

    # plot the ages of the main domestication events
    geom_vline(xintercept=c(dom1.age, dom2.age), linetype = "dashed", colour = '#c94904') +

    # set the breaks for the x-axis
    scale_x_continuous(limits = c(max_age - 150000, 0),
                       breaks = seq(max_age, 0, by = brk_width),
                       labels = seq(max_age, 0, by = brk_width)/1000,
                       minor_breaks = NULL,
                       expand = c(0.01, 0)) +

    # using limits() drops all data points that are not within the specified range,
    # causing discontinuity of the density plot, while coord_cartesian() zooms
    # without losing the data points.
    coord_cartesian(xlim = c(max_age, 0)) +

    # label the plot and the axes
    xlab("kyr BP") +
    ylab("Trait") +

    # tweak the theme
    theme_minimal(base_size = 10) +
    theme(
        # remove the vertical grid lines
        panel.grid.major.x = element_blank()
    )

dev.off()

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------


pdf(file=paste('rscript/mcmc-ridgeline-end_freq.pdf', sep=''), width = 8, height = 4.5)

ggplot() +

    # display the end_freq as a rigline plot
    stat_density_ridges(data=mcmc.params,
                        aes(x=freq, y=trait, fill=0.5 - abs(0.5-..ecdf..)),
                        scale = 3,
                        geom = "density_ridges_gradient", calc_ecdf = TRUE) +

    scale_fill_viridis(name = "Posterior", direction = -1) +

    # label the plot and the axes
    xlab("End Frequency") +
    ylab("Trait") +

    # tweak the theme
    theme_minimal(base_size = 10) +
    theme(
        # remove the vertical grid lines
        panel.grid.major.x = element_blank()
    )

dev.off()

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

pdf(file=paste('rscript/mcmc-ridgeline-alpha.pdf', sep=''), width = 8, height = 4.5)

min_x <- -50
max_x <- 250
brk_width <- 50

ggplot() +

    # display the end_freq as a rigline plot
    stat_density_ridges(data=mcmc.params,
                        aes(x=alpha1, y=trait, fill=0.5 - abs(0.5-..ecdf..)),
                        # scale = 3,
                        geom = "density_ridges_gradient", calc_ecdf = TRUE) +

    scale_fill_viridis(name = "Posterior", direction = -1) +

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
    ylab("Trait") +

    # tweak the theme
    theme_minimal(base_size = 10) +
    theme(
        # remove the vertical grid lines
        panel.grid.major.x = element_blank()
    )

dev.off()

