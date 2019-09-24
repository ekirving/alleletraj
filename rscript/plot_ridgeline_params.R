#!/usr/bin/env Rscript
library(readr, quietly = T)
library(dplyr, quietly = T)
library(ggplot2, quietly = T)
library(ggridges, quietly = T)
library(RMySQL, quietly = T)
library(stringr, quietly = T)
library(viridis, quietly = T)
library(stringr, quietly = T)
library(coda, quietly = T)
library(forcats, quietly = T)
library(grid, quietly = T)
library(gtable, quietly = T)

source('rscript/mcmc_utils.R')

# TODO should density ridges be multiplied by the number of SNPs?

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
species <- args[1]
population <- args[2]
burn_perc <- as.numeric(args[3])
pop_size <- as.numeric(args[4])
gen_time <- as.numeric(args[5])

# TODO remove when done testing
# setwd('/Users/Evan/Dropbox/Code/alleletraj')
# species <- 'horse'
# population <- 'DOM2'
# burn_perc <- 0.5
# pop_size <- 17150
# gen_time <- 8

# do not use scientific notation
options(scipen=999)

# TODO find a better way to do this
db_name <- if (species == 'horse') 'alleletraj_horse_equcab2_rel38' else 'alleletraj_cattle_umd3.1_rel38'

# connect to the remote server
mydb <- dbConnect(MySQL(), user = 'root', password = '', dbname = db_name, host = 'localhost')

# fetch all the good GWAS models
rs <- dbSendQuery(mydb, paste0(
    "SELECT DISTINCT s.modsnp_id, s.no_modern, s.mispolar, s.const_pop,
            s.no_age, s.length, s.thin, s.model, s.frac, sp.chains
       FROM selection s
       JOIN selection_psrf sp
         ON sp.selection_id = s.id
       JOIN modern_snps ms
         ON ms.id = s.modsnp_id
       JOIN qtls q
         ON q.chrom = ms.chrom
        AND q.site= ms.site
      WHERE s.population = '", population, "'
        AND s.gwas = 1
        AND sp.valid = 1
        AND q.valid = 1"
))

models <- fetch(rs, n = -1)

# load all the param files
mcmc.params <- apply(models, 1, load_models) %>% bind_rows()

# fetch the details of every GWAS association
rs <- dbSendQuery(mydb,
    "SELECT DISTINCT ms.id, t.class, t.name AS trait
       FROM qtls q
       JOIN traits t
         ON t.id = q.trait_id
       JOIN modern_snps ms
         ON ms.chrom = q.chrom
        AND ms.site = q.site
      WHERE q.valid = 1"
)

traits <- fetch(rs, n = -1)

# standardise sentance casing
traits$trait <- str_to_sentence(traits$trait)

# join the traits onto the param files
mcmc.params <- inner_join(mcmc.params, traits, by = 'id')

# add a dummy variable so we can aggregate all SNPs together
mcmc.params$allsnp <- 'GWAS hits'

# count how many SNPs each class and trait has
mcmc.params <- mcmc.params %>%
    group_by(class) %>% mutate(class.n = paste0(class, ' (n=', n_distinct(id, trait), ')')) %>%
    group_by(trait) %>% mutate(trait.n = paste0(trait, ' (n=', n_distinct(id, trait), ')')) %>%
    group_by(allsnp) %>% mutate(allsnp.n = paste0(allsnp, ' (n=', n_distinct(id, trait), ')'))

# function for plotting the MCMC params
plot_ridgeline <- function(param, xlab, min_x, max_x, brk_w, lim_x = NULL, x_breaks = NULL,
                        x_labels = NULL, bandwidth = NULL, vline = NULL) {

    # set default limits, breaks and labels
    if (is.null(lim_x)) {
        lim_x <- brk_w
    }

    if (is.null(x_breaks)) {
        x_breaks <- seq(min_x, max_x, by = brk_w)
    }

    if (is.null(x_labels)) {
        x_labels <- x_breaks
    }

    # TODO consider -ve log10 transform
    # # https://stackoverflow.com/questions/37446064/i-need-ggplot-scale-x-log10-to-give-me-both-negative-and-positive-numbers-as-o
    # signed_log <- scales::trans_new("signed_log", transform=function(x) sign(x)*log(abs(x)), inverse=function(x) sign(x)*exp(abs(x)))
    # 
    # # https://stackoverflow.com/questions/14504869/histogram-with-negative-logarithmic-scale-in-r?noredirect=1&lq=1
    # trans_asinh <- scales::trans_new("asinh", transform = function(x) asinh(x),  inverse = function(x) sinh(x))

    # setup shared layers for the selection plots
    gglayers <- list(

        # add optional vertial line(s) to the plot
        vline,

        # set the colour scheme for the density plot
        scale_fill_viridis(name = "Posterior", direction = -1),

        # set the breaks for the x-axis
        scale_x_continuous(limits = c(min_x - lim_x, max_x),
                           breaks = x_breaks,
                           labels = x_labels,
                           minor_breaks = NULL,
                           expand = c(0.01, 0)),

        # using limits() drops all data points that are not within the specified range,
        # causing discontinuity of the density plot, while coord_cartesian() zooms
        # without losing the data points.
        coord_cartesian(xlim = c(min_x, max_x)),

        # label the x axis
        xlab(xlab),

        # tweak the theme
        theme_minimal(base_size = 10),

        # remove the vertical grid lines
        theme(panel.grid.major.x = element_blank())
    )

    # --------------------------------------------------------------------------
    # joint density plot
    # --------------------------------------------------------------------------

    pdf(file = paste0('data/pdf/plots/', species, '-', population ,'-ridgeline-', param, '-all.pdf'),
        width = 16, height = 4.5)

    # built the plot, but don't display it yet
    g <- ggplot() +

        # add the density plot
        stat_density_ridges(
            data = mcmc.params,
            aes_string(x = param, y = 'allsnp.n', fill = '0.5 - abs(0.5-..ecdf..)'),
            geom = "density_ridges_gradient", calc_ecdf = TRUE,
            bandwidth = bandwidth,  # controls smoothing in density plot
            rel_min_height = 0.005, # trim the trailing lines
            scale = 1000000         # set massive scale so single ridgeline uses all space
        ) +

        # label the y-axis
        ylab("Density") +

        # add the common layers
        gglayers

    # draw the finished plot
    print(g)

    dev.off()

    # --------------------------------------------------------------------------
    # classes as ridgelines
    # --------------------------------------------------------------------------

    pdf(file = paste0('data/pdf/plots/', species, '-', population ,'-ridgeline-', param, '-classes.pdf'),
        width = 16, height = 9)

    # built the plot, but don't display it yet
    g <- ggplot() +

        # add the trait classes
        stat_density_ridges(
            data = mcmc.params,
            aes_string(x = param, y = 'fct_rev(class.n)', fill = '0.5 - abs(0.5-..ecdf..)'),
            geom = "density_ridges_gradient", calc_ecdf = TRUE,
            bandwidth = bandwidth,  # controls smoothing in density plot
            rel_min_height = 0.01,  # trim the trailing lines
            scale = 2               # controls vertival overlap
        ) +

        # label the y-axis
        ylab("Class") +

        # add the common layers
        gglayers

    # draw the finished plot
    print(g)

    dev.off()

    # --------------------------------------------------------------------------
    # traits as ridgelines, facted by class
    # --------------------------------------------------------------------------

    pdf(file = paste0('data/pdf/plots/', species, '-', population ,'-ridgeline-', param, '-traits.pdf'),
        width = 16, height = 9)

    # built the plot, but don't display it yet
    g <- ggplot() +

        # add the traits
        stat_density_ridges(
            data = mcmc.params,
            aes_string(x = param, y = 'fct_rev(trait.n)', fill = '0.5 - abs(0.5-..ecdf..)'),
            geom = "density_ridges_gradient", calc_ecdf = TRUE,
            bandwidth = bandwidth,  # controls smoothing in density plot
            rel_min_height = 0.01,  # trim the trailing lines
            scale = 2               # controls vertival overlap
        ) +

        # use facets to split the traits by class
        facet_wrap(facets = vars(class.n), ncol = 1, scales = "free_y") +

        # label the y-axis
        ylab("Trait") +

        # add the common layers
        gglayers

    # make a table object to control the facet layout
    gt = ggplot_gtable(ggplot_build(g))

    # count how many traits each class has
    ntraits <- mcmc.params %>%
        group_by(class) %>%
        summarise(ntraits = n_distinct(trait))

    # set the heights of each facet relative to the number of traits
    for (i in 1:nrow(ntraits)) {
        row <- i * 5 + 3
        height <- ntraits[i,2]/5
        gt$heights[row] = unit(height, units = "null")
    }

    # helper funtion which shows the margins of the layout
    # gtable_show_layout(gt)

    # draw the finished plot
    grid.draw(gt)

    dev.off()
}

# ------------------------------------------------------------------------------
# plot the Age of the allele
# ------------------------------------------------------------------------------

# convert diffusion units into calendar years
mcmc.params$ageyrs <- mcmc.params$age * 2 * pop_size * gen_time

# define the limits of the plot
min_x <- -50000
max_x <- 0
brk_w <- 5000
lim_x <- abs(min_x) * 10
bandw <- 1000

# show x-axis label in units of kya
x_breaks <- seq(min_x, max_x, by = brk_w)
x_labels <- x_breaks / 1000

# add a dashed vertical line for the two main domestication timepoints
vline <- geom_vline(xintercept = c(-5500, -4000), linetype = "dashed", 
                    colour = '#c94904')

plot_ridgeline(param = 'ageyrs', xlab = 'kyr BP', min_x = min_x, max_x = max_x, 
               brk_w = brk_w, lim_x = lim_x, x_breaks = x_breaks, 
               x_labels = x_labels, bandwidth = bandw, vline = vline)

# ------------------------------------------------------------------------------
# plot the End frequency
# ------------------------------------------------------------------------------

# convert end_freq back into an actual frequency
mcmc.params$freq <- (1 - cos(mcmc.params$end_freq)) / 2

# define the limits of the plot
min_x <- -0.03
max_x <- 1
brk_w <- 0.25
lim_x <- 0

# show x-axis label in units of 25%
x_breaks <- seq(0, 1, by = 0.25)
x_labels <- seq(0, 1, by = 0.25)

# TODO fix issue with negative values in the density estimage due to many near zero value
plot_ridgeline(param = 'freq', xlab = 'End Frequency', min_x = min_x, 
               max_x = max_x, brk_w = brk_w, lim_x = lim_x, 
               x_breaks = x_breaks, x_labels = x_labels)

# ------------------------------------------------------------------------------
# plot the Selection coefficients
# ------------------------------------------------------------------------------

# convert alpha params into s values
mcmc.params$s1 <- mcmc.params$alpha1 / (2 * pop_size)
mcmc.params$s2 <- mcmc.params$alpha2 / (2 * pop_size)

# define the limits of the plot
min_x <- -0.005
max_x <- 0.015
brk_w <- 0.005

# add a dashed vertical line at 0
vline <- geom_vline(xintercept = 0, linetype = "dashed", colour = '#c94904')

# plot s1
plot_ridgeline(param = 's1', xlab = expression(paste("s"[1])), min_x = min_x, 
               max_x = max_x, brk_w = brk_w, vline = vline)

# plot s2
plot_ridgeline(param = 's2', xlab = expression(paste("s"[2])), min_x = min_x, 
               max_x = max_x, brk_w = brk_w, vline = vline)
