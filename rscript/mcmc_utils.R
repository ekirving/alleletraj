#!/usr/bin/env Rscript

load_chains <- function(param_files, burn_perc, thin, verbose=TRUE) {
    chains <- c()

    start <- NA

    for (i in param_files) {
        # load the chain (and cast infinite values to NA)
        chain <- fread(i, header = T, sep = '\t', drop=c('gen', 'first_nonzero'), na.strings=c('inf', '-inf'))

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
            cat(paste0("WARNING: Chains are not of equal length. min=", min_iter, " max=", max_iter, " (", round(min_iter/max_iter*100, 0), "%)"))
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
