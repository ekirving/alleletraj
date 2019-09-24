#!/usr/bin/env Rscript
library(data.table, quietly=T)

load_chains <- function(param_files, burn_perc, thin=100, select=NULL, drop=NULL, verbose=TRUE) {
    chains <- c()
    start <- NA

    for (i in param_files) {
        # load the chain (and cast infinite values to NA)
        chain <- fread(i, header = T, sep = '\t', select=select, drop=drop, na.strings=c('inf', '-inf'))

        # drop NAs
        chain <- na.omit(chain)

        chain.length <- nrow(chain)

        # convert burn % to number of records
        burnin = round(burn_perc * chain.length)

        if (is.na(start)) {
            # make sure we're consistent about the start point
            start <- burnin * thin
        }

        # burn in the chain (thinning is already done)
        chains[[i]] <- mcmc(
            chain[(burnin + 1):chain.length,],
            start = start,
            thin = thin)
    }

    # check if chains are all equal length
    min_iter <- min(unlist(lapply(chains, niter)))
    max_iter <- max(unlist(lapply(chains, niter)))

    # normalise chain length (necessary if some models did not run to completion)
    if (min_iter != max_iter) {
        if (verbose) {
            cat(paste0("WARNING: Chains are not of equal length. min=", min_iter, " max=", max_iter, " (", round(min_iter/max_iter*100, 0), "%)", "\n"))
        }
        chains <- lapply(chains, function(x) {
            num_iter <- niter(x)
            offset <- num_iter - min_iter + 1
            mcmc(x[offset:num_iter,],
                 start = round(burn_perc/(1-burn_perc) * min_iter),
                 thin = thin)
        })
    }

    mcmc.list(chains)
}

# compose the path to the param files from the hyperparameters of the model
param_file_paths <- function(model) {
    path <- paste0('data/selection/', species, '-', population, '-modsnp', trimws(model[['modsnp_id']]))
    for (flag in c('no_modern', 'mispolar', 'const_pop', 'no_age')) {
        if (model[flag] == 1) {
            path <- paste0(path, '-', flag)
        }
    }
    path <- paste0(path, '-n', model['length'], '-s', model['thin'], '-h', model['model'], '-F', model['frac'])
    paste0(path, '-chain', unlist(strsplit(model[['chains']], ",")), '.param.gz')
}

load_models <- function(model) {
    param_files <- param_file_paths(model)

    # ignore all the time params
    select <- c('lnL', 'pathlnL', 'alpha1', 'alpha2', 'F', 'age','end_freq', 'first_nonzero')

    # load all the chains
    chains.all <- load_chains(param_files, burn_perc, select=select, verbose=F)

    # merge the replicates
    chains.all <- as.data.frame(rbind(chains.all[[1]], chains.all[[2]]))

    # add the modsnp ID
    chains.all$id <- as.numeric(model[['modsnp_id']])

    chains.all
}
