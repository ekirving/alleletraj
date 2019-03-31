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

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
species <- args[1]
population <- args[2]
burnin <- strtoi(args[3])
pop_size <- strtoi(args[4])
gen_time <- strtoi(args[5])

# TODO remove when done testing
# setwd('/Users/Evan/Dropbox/Code/alleletraj')
# species <- 'horse'
# population <- 'DOM2'
# burnin <- 1000 # i.e. 20%
# pop_size <- 16000
# gen_time <- 8

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
files <- list.files(path = "selection", pattern = paste(species, population, "(.+).param", sep = '-'), full.names = T)
mcmc.params <- lapply(files, read_tsv_param) %>% bind_rows()

# add a dummy variable to make ggridges behave like a normal density plot
mcmc.params$density <- ''

# connect to the remote server
mydb <- dbConnect(MySQL(), user = 'root', password = '', dbname = 'alleletraj_horse', host = 'localhost')

# TODO this is explicitly dropping all non-GWAS hits, which is fine for now but needs updating later
# fetch the details of the SNP
rs <- dbSendQuery(mydb,
    "SELECT DISTINCT ms.id, t.class, t.name AS trait
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

traits <- fetch(rs, n = -1)

# standardise sentance casing
traits$trait <- str_to_sentence(traits$trait)

# join the traits onto the param files
mcmc.params <- inner_join(mcmc.params, traits, by = 'id')

# function for plotting the MCMC params
plot_params <- function(param, xlab, min_x, max_x, brk_w, lim_x = NULL, x_breaks = NULL,
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

    pdf(file = paste0('rscript/pdf/', species, '-', population ,'-ridgeline-', param, '-all.pdf'), width = 16, height = 4.5)

    # built the plot, but don't display it yet
    g <- ggplot() +

        # add the density plot
        stat_density_ridges(
            data = mcmc.params,
            aes_string(x = param, y = 'density', fill = '0.5 - abs(0.5-..ecdf..)'),
            geom = "density_ridges_gradient", calc_ecdf = TRUE,
            bandwidth = bandwidth,  # controls smoothing in density plot
            rel_min_height = 0.005, # trim the trailing lines
            scale = 10000           # large vertial overlap for single ridgeline
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

    # scale height relative to the number of classes
    pdf.height <- length(unique(mcmc.params$class)) * 0.75

    pdf(file = paste0('rscript/pdf/', species, '-', population ,'-ridgeline-', param, '-classes.pdf'), width = 16, height = pdf.height)

    # built the plot, but don't display it yet
    g <- ggplot() +

        # add the trait classes
        stat_density_ridges(
            data = mcmc.params,
            aes_string(x = param, y = 'fct_rev(class)', fill = '0.5 - abs(0.5-..ecdf..)'),
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

    # scale height relative to the number of traits
    pdf.height <- length(unique(mcmc.params$trait)) * 0.225

    pdf(file = paste0('rscript/pdf/', species, '-', population ,'-ridgeline-', param, '-traits.pdf'), width = 16, height = pdf.height)

    # built the plot, but don't display it yet
    g <- ggplot() +

        # add the traits
        stat_density_ridges(
            data = mcmc.params,
            aes_string(x = param, y = 'fct_rev(trait)', fill = '0.5 - abs(0.5-..ecdf..)'),
            geom = "density_ridges_gradient", calc_ecdf = TRUE,
            bandwidth = bandwidth,  # controls smoothing in density plot
            rel_min_height = 0.01,  # trim the trailing lines
            scale = 2               # controls vertival overlap
        ) +

        # use facets to split the traits by class
        facet_wrap(facets = vars(class), ncol = 1, scales = "free_y") +

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

# show x-axis label in units of kya
x_breaks <- seq(-50000, 0, by = 5000)
x_labels <- x_breaks/1000

# add a dashed vertical line for the two main domestication timepoints
vline <- geom_vline(xintercept = c(-5500, -4000), linetype = "dashed", colour = '#c94904')

plot_params(param = 'ageyrs', xlab = 'kyr BP', min_x = -50000, max_x = 0, brk_w = 5000, lim_x = 150000, x_breaks = x_breaks, x_labels = x_labels, bandwidth = 1000, vline = vline)

# ------------------------------------------------------------------------------
# plot the End frequency
# ------------------------------------------------------------------------------

# convert end_freq back into an actual frequency
mcmc.params$freq <- (1-cos(mcmc.params$end_freq))/2

x_breaks <- seq(0, 1, by = 0.25)
x_labels <- seq(0, 1, by = 0.25)

# TODO fix issue with negative values in the density estimage due to many near zero value
plot_params(param = 'freq', xlab = 'End Frequency', min_x = -0.03, max_x = 1, brk_w = 0.25, lim_x = 0, x_breaks = x_breaks, x_labels = x_labels)

# ------------------------------------------------------------------------------
# plot the Selection coefficients
# ------------------------------------------------------------------------------

# convert alpha params into s values
mcmc.params$s1 <- mcmc.params$alpha1 / (2 * pop_size)
mcmc.params$s2 <- mcmc.params$alpha2 / (2 * pop_size)

# add a dashed vertical line at 0
vline <- geom_vline(xintercept = 0, linetype = "dashed", colour = '#c94904')

plot_params(param = 's1', xlab = expression(paste("s"[1])), min_x = -0.005, max_x =  0.015, brk_w = 0.005, vline = vline)
plot_params(param = 's2', xlab = expression(paste("s"[2])), min_x = -0.005, max_x =  0.015, brk_w = 0.005, vline = vline)
