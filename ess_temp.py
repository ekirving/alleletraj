#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas
from py_coda import mcmc, read_pycoda

chain_file = "selection.old/horse-DOM2-1362521.param"
dat_prefix = "selection.old/horse-DOM2-1362521"
pdf_prefix = "pdf.old/horse-DOM2-1362521-"
burn_in = 1000

df = pandas.read_csv(chain_file, header=0, sep="\t")

mcmcobj = mcmc(df.columns[1:], df.loc[burn_in:, 'lnL':].values, thin=1)

mcmcobj.plot_traces(savefig=pdf_prefix)
mcmcobj.plot_autocorr(savefig=pdf_prefix)

# see https://pymc-devs.github.io/pymc/modelchecking.html#formal-methods
with open(dat_prefix + ".geweke", "w") as fout:
    mcmcobj.geweke(savefig=pdf_prefix, fp=fout)

with open(dat_prefix + ".stats", "w") as fout:
    mcmcobj.get_stats(fp=fout)

with open(dat_prefix + ".heidel", "w") as fout:
    mcmcobj.heidelberger_welch(fp=fout)
