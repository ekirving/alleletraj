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

setwd('/Users/Evan/Dropbox/Code/alleletraj')

# TODO command line params
pop_size <- 16000
gen_time <- 8
burnin <- 1000 # i.e. 20%
species <- 'horse'
population <- 'DOM2'


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
files <- list.files(path="selection", pattern=paste(species, population, "(.+).param", sep = '-'), full.names = T)
mcmc.params <- lapply(files, read_tsv_param) %>% bind_rows()

# convert diffusion units into calendar years
mcmc.params$ageyrs <- mcmc.params$age * 2 * pop_size * gen_time

# convert end_freq back into an actual frequency
mcmc.params$freq <- (1-cos(mcmc.params$end_freq))/2

# convert alpha params into s
mcmc.params$s1 <- mcmc.params$alpha1 / (2 * pop_size )
mcmc.params$s2 <- mcmc.params$alpha2 / (2 * pop_size )

# add a dummy variable to allow single line density plots
mcmc.params$density <- ''

# connect to the remote server
mydb <- dbConnect(MySQL(), user='root', password='', dbname='alleletraj_horse', host='localhost')

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

traits <- fetch(rs, n=-1)

# standardise sentance casing
traits$trait <- str_to_sentence(traits$trait)

# join the traits onto the param files
mcmc.params <- inner_join(mcmc.params, traits, by = 'id')

# ------------------------------------------------------------------------------
# --                           Age                                            --
# ------------------------------------------------------------------------------

dom1.age <- -5500
dom2.age <- -4000

max_age <- -50000
brk_width <- 5000

# all traits...

pdf(file=paste('rscript/pdf/mcmc-age.pdf', sep=''), width = 8, height = 4.5)

ggplot() +

    # display the sample dates as a rigline plot
    stat_density_ridges(data=mcmc.params,
                        aes(x=ageyrs, y=density, fill=0.5 - abs(0.5-..ecdf..)),
                        scale = 100000,
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
    ylab("Density") +

    # tweak the theme
    theme_minimal(base_size = 10) +
    theme(
        # remove the vertical grid lines
        panel.grid.major.x = element_blank()
    )

dev.off()


# trait classes as ridgelines...

pdf(file=paste('rscript/pdf/mcmc-age-ridgeline.pdf', sep=''), width = 8, height = 4.5)

ggplot() +

    # display the sample dates as a rigline plot
    stat_density_ridges(data=mcmc.params,
                        aes(x=ageyrs, y=class, fill=0.5 - abs(0.5-..ecdf..)),
                        scale = 2,
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


# traits from each class as ridgelines...

classes <- unique(mcmc.params$class)

for (cls in classes) {

    # count the number of traits
    num_traits <- length(unique(mcmc.params[mcmc.params$class == cls,]$trait))

    pdf(file=paste0('rscript/pdf/mcmc-age-', str_to_lower(cls), '-ridgeline.pdf'), width = 8, height = max(num_traits * 0.8, 1.6))

    print(ggplot() +

        # display the data as a rigline plot
        stat_density_ridges(data=mcmc.params[mcmc.params$class == cls,],
                            aes(x=ageyrs, y=trait, fill=0.5 - abs(0.5-..ecdf..)),
                            scale = ifelse(num_traits > 1, 1.5, 100000),
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
        ylab(cls) +

        # tweak the theme
        theme_minimal(base_size = 10) +
        theme(
            # remove the vertical grid lines
            panel.grid.major.x = element_blank()
        ))

    dev.off()
}

# ------------------------------------------------------------------------------
# --                            End Freq                                      --
# ------------------------------------------------------------------------------

# all traits...

pdf(file=paste('rscript/pdf/mcmc-end_freq.pdf', sep=''), width = 8, height = 4.5)

ggplot() +

    # display the end_freq as a rigline plot
    stat_density_ridges(data=mcmc.params,
                        aes(x=freq, y=density, fill=0.5 - abs(0.5-..ecdf..)),
                        scale = 2,
                        geom = "density_ridges_gradient", calc_ecdf = TRUE) +

    scale_fill_viridis(name = "Posterior", direction = -1) +

    # label the plot and the axes
    xlab("End Frequency") +
    ylab("Density") +

    # tweak the theme
    theme_minimal(base_size = 10) +
    theme(
        # remove the vertical grid lines
        panel.grid.major.x = element_blank()
    )

dev.off()


# traits as ridgelines...

pdf(file=paste('rscript/pdf/mcmc-end_freq-ridgeline.pdf', sep=''), width = 8, height = 4.5)

ggplot() +

    # display the end_freq as a rigline plot
    stat_density_ridges(data=mcmc.params,
                        aes(x=freq, y=class, fill=0.5 - abs(0.5-..ecdf..)),
                        scale = 2,
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
# --                                S1                                        --
# ------------------------------------------------------------------------------

min_x <- -0.005
max_x <-  0.015
brk_width <- 0.005

# the common settings for the selection plots
gglayers <- list(

    # add a dashed vertical line at 0
    geom_vline(xintercept=0, linetype = "dashed", colour = '#c94904'),

    # set the colour scheme for the density plot
    scale_fill_viridis(name = "Posterior", direction = -1),

    # set the breaks for the x-axis
    scale_x_continuous(limits = c(min_x-brk_width, max_x+brk_width),
                       breaks = seq(min_x, max_x, by = brk_width),
                       labels = seq(min_x, max_x, by = brk_width),
                       minor_breaks = NULL,
                       expand = c(0.01, 0)),

    # using limits() drops all data points that are not within the specified range,
    # causing discontinuity of the density plot, while coord_cartesian() zooms
    # without losing the data points.
    coord_cartesian(xlim = c(min_x, max_x)),

    # label the x axis (using subscript notation)
    xlab(ifelse(s=='s1',
                expression(paste("s"[1])),
                expression(paste("s"[2])))),

    # tweak the theme
    theme_minimal(base_size = 10),

    # remove the vertical grid lines
    theme(panel.grid.major.x = element_blank())
)

for (s in c('s1','s2')) {

    # --------------------------------------------------------------------------
    # all data as a density plot
    # --------------------------------------------------------------------------

    pdf(file=paste0('rscript/pdf/', species, '-', population ,'-ridgeline-', s, '-all.pdf'), width = 16, height = 4.5)

    # built the plot, but don't display it yet
    g <- ggplot() +

        # add the density plot
        stat_density_ridges(data=mcmc.params,
                            aes_string(x=s, y='density', fill='0.5 - abs(0.5-..ecdf..)'),
                            geom = "density_ridges_gradient", calc_ecdf = TRUE,
                            rel_min_height = 0.005, # trim the trailing lines
                            scale = 2               # controls vertival overlap
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

    pdf(file=paste0('rscript/pdf/', species, '-', population ,'-ridgeline-', s, '-classes.pdf'), width = 16, height = pdf.height)

    # built the plot, but don't display it yet
    g <- ggplot() +

        # add the trait classes
        stat_density_ridges(data=mcmc.params,
                            aes_string(x=s, y='fct_rev(class)', fill='0.5 - abs(0.5-..ecdf..)'),
                            geom = "density_ridges_gradient", calc_ecdf = TRUE,
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

    pdf(file=paste0('rscript/pdf/', species, '-', population ,'-ridgeline-', s, '-traits.pdf'), width = 16, height = pdf.height)

    # built the plot, but don't display it yet
    g <- ggplot() +

        # add the traits
        stat_density_ridges(data=mcmc.params,
                            aes_string(x=s, y='fct_rev(trait)', fill='0.5 - abs(0.5-..ecdf..)'),
                            geom = "density_ridges_gradient", calc_ecdf = TRUE,
                            rel_min_height = 0.01,  # trim the trailing lines
                            scale = 2               # controls vertival overlap
        ) +

        # use facets to split the traits by class
        facet_wrap(facets=vars(class), ncol=1, scales="free_y") +

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
