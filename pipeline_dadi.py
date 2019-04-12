#!/usr/bin/env python
# -*- coding: utf-8 -*-

import dadi
import luigi
import math
import pickle
import random
import scipy.stats as st

# import my custom modules
from pipeline_consts import *
from pipeline_snp_call import PolarizeVCF, SubsetSNPsVCF
from pipeline_utils import PipelineTask, run_cmd

# number of sequential epochs to test
DADI_EPOCHS = 5

# how many independent runs should we do to find global maximum of params
DADI_REPLICATES = 50

# number of points to use in the grid
DADI_GRID_PTS = 100


def dadi_epoch_model(params, ns, pts):
    """
    Sequential epoch model for dadi.

    :param params: Times and population sizes of the epochs (e.g. t1, t2, ..., n1, n2)
    :param ns: The number of samples
    :param pts: The number of points in the grid
    """

    # each epoch has a time and a size
    num_epoch = len(params) / 2

    # split the two types of params
    time = params[num_epoch:]
    size = params[:num_epoch]

    # make the grid
    grid = dadi.Numerics.default_grid(pts)

    # one-dimensional phi for a constant-sized population
    phi = dadi.PhiManip.phi_1D(grid)

    for i in range(num_epoch):
        # integrate a 1-dimensional phi forward
        phi = dadi.Integration.one_pop(phi, grid, time[i], size[i])

    # compute sample Spectrum from population frequency distribution phi
    fs = dadi.Spectrum.from_phi(phi, ns, [grid])

    return fs


class EasySFS(PipelineTask):
    """
    Calculate the Site Frequency Spectrum.

    :type species: str
    :type population: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    folded = luigi.BoolParameter()

    def requires(self):
        return SubsetSNPsVCF(self.species, self.population)

    def output(self):
        return luigi.LocalTarget('sfs/{}/dadi/{}.sfs'.format(self.basename, self.population))

    def run(self):

        # get all the samples to use
        samples = [s for s in SAMPLES[self.species][self.population] if s not in SFS_EXCLUSIONS[self.species]]

        # make a sample/population file
        pop_file = 'sfs/{}.pops'.format(self.basename)
        with open(pop_file, 'w') as fout:
            for sample in samples:
                fout.write('{}\t{}\n'.format(sample, self.population))

        params = {
            'vcf':  self.input().path,
            'pops': pop_file,
            'out':  self.basename,
            'proj': len(samples) * 2,  # don't project down
            'fold': '--unfolded' if not self.folded else ''
        }

        # pipe 'yes' into easySFS to get past the interactive prompt which complains about excluded samples
        cmd = "echo 'yes' | easySFS.py -a -i {vcf} -p {pops} -o sfs/{out} --proj {proj} {fold}".format(**params)

        run_cmd([cmd], shell=True)


class DadiEpochOptimizeParams(PipelineTask):
    """
    Optimise the log likelihood of the model parameters for the given SFS.

    :type species: str
    :type population: str
    :type folded: bool
    :type epoch: int
    :type n: int
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    folded = luigi.BoolParameter()
    epoch = luigi.IntParameter()
    n = luigi.IntParameter()

    def requires(self):
        return EasySFS(self.species, self.population, self.folded)

    def output(self):
        return [luigi.LocalTarget("sfs/{}.{}".format(self.basename, ext)) for ext in ['pickle', 'log']]

    def run(self):

        # unpack the inputs/outputs
        sfs_file = self.input()
        pkl_file, log_file = self.output()

        # load the frequency spectrum
        fs = dadi.Spectrum.from_file(sfs_file.path)

        # set the upper and lower parameter bounds
        lower = [.01] * self.epoch + [0] * self.epoch
        upper = [100] * self.epoch + [5] * self.epoch

        # TODO pick random starting values (bounded by lower/upper)
        # start = [random.uniform(lower[i], upper[i]) for i in range(0, self.epoch * 2)]

        # pick random starting values (bounded by 0-2)
        start = st.uniform.rvs(scale=2, size=self.epoch * 2)

        # optimize log(params) to fit model to data using Nelder-Mead algorithm
        p_opt = dadi.Inference.optimize_log_fmin(start, fs, dadi_epoch_model, lower_bound=lower, upper_bound=upper,
                                                 pts=DADI_GRID_PTS, verbose=50, output_file=log_file.path,
                                                 full_output=True)

        # fit the optimised model
        model = dadi_epoch_model(p_opt[0], fs.sample_sizes, DADI_GRID_PTS)

        # calculate theta given the model
        theta = dadi.Inference.optimal_sfs_scaling(model, fs)

        # save the relevant information
        best = {'epoch': self.epoch, 'n': self.n, 'lnL': -p_opt[1], 'params': p_opt[0], 'theta': theta}

        # save the results by pickling them in a file
        with pkl_file.open('w') as fout:
            pickle.dump(best, fout)


class DadiEpochMaximumLikelihood(PipelineTask):
    """
    Find the model run with the maximum log likelihood out of all the replicates.

    Because dadi gets stuck easily in local maxima we run multiple replicates.

    :type species: str
    :type population: str
    :type folded: bool
    :type epoch: int
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    folded = luigi.BoolParameter()
    epoch = luigi.IntParameter()

    def requires(self):
        for n in range(1, DADI_REPLICATES + 1):
            yield DadiEpochOptimizeParams(self.species, self.population, self.folded, self.epoch, n)

    def output(self):
        return luigi.LocalTarget("sfs/{}-maxlnL.pickle".format(self.basename))

    def run(self):

        params = []

        # load all the pickled param values from the replicate runs
        for pkl_file, log_file in self.input():
            with pkl_file.open('r') as fin:
                params.append(pickle.load(fin))

        # find the params that produced the highest maximum likelihood
        max_lnl = max(params, key=lambda x: x['lnL'])

        # save the results by pickling them in a file
        with self.output().open('w') as fout:
            pickle.dump(max_lnl, fout)


class DadiEpochBestModel(PipelineTask):
    """
    Find the best fitting model across all epochs, using the Akaike information criterion (AIC).

    :type species: str
    :type population: str
    :type folded: bool
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    folded = luigi.BoolParameter()

    def requires(self):
        for epoch in range(1, DADI_EPOCHS + 1):
            yield DadiEpochMaximumLikelihood(self.species, self.population, self.folded, epoch)

    def output(self):
        return luigi.LocalTarget("sfs/{}.pickle".format(self.basename))

    def run(self):
        epochs = []

        # load the pickled max lnL params from each of the epoch models
        for pkl_file in self.input():
            with pkl_file.open('r') as fin:
                epochs.append(pickle.load(fin))

        # calculate the AIC - i.e. AIC = 2 * num_params - 2 * lnL
        for epoch in epochs:
            epoch['aic'] = (2 * len(epoch['params'])) - (2 * epoch['lnL'])

        # get the min AIC
        min_aic = min(e['aic'] for e in epochs)

        # compute the relative likelihood of each model - i.e. exp((AICmin − AICi)/2)
        for epoch in epochs:
            epoch['relL'] = math.exp((min_aic - epoch['aic']) / 2)

        # reverse sort by relative likelihood
        epochs.sort(key=lambda x: x['relL'], reverse=True)

        if epochs[1]['relL'] > 0.05:
            raise Exception('ERROR: Cannot reject second best model based on relative likelihood of AIC')

        # save the results by pickling them in a file
        with self.output().open('w') as fout:
            pickle.dump(epochs, fout)


class CountCallableSites(PipelineTask):
    """
    Count the number of callable sites, as dadi needs this number to estimate the ancestral population size from theta.

    :type species: str
    :type population: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()

    def requires(self):
        return PolarizeVCF(self.species, self.population)

    def output(self):
        return luigi.LocalTarget('sfs/{}.L'.format(self.basename))

    def run(self):

        # count all unique sites
        size = run_cmd(["bcftools query -f '%CHROM %POS\\n' {} | uniq | wc -l".format(self.input().path)], shell=True)

        with self.output().open('w') as fout:
            fout.write(size)


class DadiModelDemography(luigi.WrapperTask):
    """
    Find the best fitting of 5 sequential epoch models (i.e. 1 epoch, 2 epoch, etc.).
    """

    def requires(self):
        # run multiple epochs, with multiple independent replicates
        return DadiEpochBestModel('horse', 'DOM2', False)


if __name__ == '__main__':
    luigi.run()
