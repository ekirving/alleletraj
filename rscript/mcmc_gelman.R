
cat("\n\n", "--------------", "\n\n")
cat("Analysing combined chains.", "\n\n")

chains.all = mcmc.list(chains)

# print the summary stats
print(summary(chains.all))

cat("Effective Sample Size.\n")
print(effectiveSize(chains.all))
cat("\n")

# plot the combined traces
cat("Plotting combined traces.", "\n\n")
# pdf(file=paste0('bayes/', prefix, "-", graph_code, '-burn-trace-0.pdf'))
png(file=paste0('bayes/', prefix, "-", graph_code, '-burn-trace-0-pt%d.png'), width=7, height=7, units='in', res=300)
plot(chains.all)
off <- dev.off()

cat("Plotting the Gelman and Rubin's convergence diagnostic.", "\n\n")
pdf(file=paste0('bayes/', prefix, "-", graph_code, '-burn-gelman.pdf'))
gelman.plot(chains.all)
off <- dev.off()

# NB. values substantially above 1 indicate lack of convergence.
gelman <- gelman.diag(chains.all, multivariate=FALSE, autoburnin=FALSE)

cat("Gelman and Rubin's convergence diagnostic.", "\n")
print(gelman)

if (gelman$psrf['likelihood', 1] > 1.1) {
    cat(paste0("WARNING: PSRF of likelihood above threshold = ", round(gelman$psrf['likelihood', 1], 3),
               ' ./bayes/', prefix, '-', graph_code))
}
