#!/usr/bin/env python
# -*- coding: utf-8 -*-

import luigi
import math
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import dadi
import pickle
import random

# import my custom modules
from alleletraj.consts import MUTATION_RATE
from snp_call import PolarizeVCF, WholeGenomeSNPsVCF
from alleletraj import utils

# number of sequential epochs to test
DADI_EPOCHS = 5

# how many independent replicates should we run to find the global maximum params (dadi can get stuck in local maxima)
DADI_REPLICATES = 500

# number of points to use in the grid
DADI_GRID_PTS = 100

# maximum relative log likelihood to not reject the second best model
DADI_MAX_RELATIVE_LL = 0.10


def dadi_n_epoch(params, ns, pts):
    """
    Sequential epoch model for dadi.

    :param params: Population sizes and times of the epochs (e.g. n1, n2... t1, t2)
    :param ns: Number of samples in resulting Spectrum
    :param pts: Number of grid points to use in integration
    """

    # how many epochs does this model have (each epoch has 2 params)
    epochs = len(params) / 2

    # nu: Ratio of contemporary to ancient population size
    # t:  Time in the past at which size change happened (in units of 2*Na generations)
    nu, t = params[:epochs], params[epochs:]

    # make the grid
    grid = dadi.Numerics.default_grid(pts)

    # one-dimensional phi for a constant-sized population
    phi = dadi.PhiManip.phi_1D(grid)

    for i in range(epochs):
        # integrate a 1-dimensional phi forward
        phi = dadi.Integration.one_pop(phi, grid, t[i], nu[i])

    # compute sample Spectrum from population frequency distribution phi
    fs = dadi.Spectrum.from_phi(phi, ns, [grid])

    return fs


class EasySFS(utils.PipelineTask):
    """
    Calculate the Site Frequency Spectrum.

    :type species: str
    :type population: str
    :type folded: bool
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    folded = luigi.BoolParameter()

    resources = {'cpu-cores': 1, 'ram-gb': 96}

    def requires(self):
        return WholeGenomeSNPsVCF(self.species)

    def output(self):
        return [luigi.LocalTarget('data/sfs/{}/dadi/{}.{}'.format(self.basename, self.population, ext))
                for ext in ['sfs', 'log']]

    def run(self):
        # unpack the outputs
        sfs_file, log_file = self.output()

        # get all the samples to use (excluding those not explicitly flagged for the SFS calculation)
        samples = [sample for sample in self.samples if self.samples[sample]['sfs'] == '1']

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

        # NOTE easySFS expects the REF allele to be ancestral, rather than using the INFO/AA field
        # pipe 'yes' into easySFS to get past the interactive prompt which complains about excluded samples
        cmd = "echo 'yes' | easySFS.py -a -f -i {vcf} -p {pops} -o sfs/{out} --proj {proj} {fold}".format(**params)

        log = utils.run_cmd([cmd], shell=True)

        # save the output
        with log_file.open('w') as fout:
            fout.write(log)


class DadiEpochOptimizeParams(utils.PipelineTask):
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
        # trim the n value from the folder name
        folder = self.basename.rpartition('-')[0]
        return [luigi.LocalTarget('data/dadi/{}/{}.{}'.format(folder, self.basename, ext)) for ext in ['pkl', 'log']]

    def run(self):
        # unpack the inputs/outputs
        sfs_file, _ = self.input()
        pkl_file, log_file = self.output()

        # load the frequency spectrum
        fs = dadi.Spectrum.from_file(sfs_file.path)

        # set the upper and lower parameter bounds (0.01 < nu < 100 | 0 < T < 5)
        lower = [.01] * self.epoch + [0] * self.epoch
        upper = [100] * self.epoch + [5] * self.epoch

        # pick random starting values, bounded by the upper and lower parameter limits
        start = [random.uniform(lower[i], upper[i]) for i in range(0, self.epoch * 2)]

        # make the extrapolating version of our demographic model function.
        dadi_n_epoch_extrap = dadi.Numerics.make_extrap_log_func(dadi_n_epoch)

        # make sure the output folder exists
        log_file.makedirs()

        # optimize log(params) to fit model to data using Nelder-Mead algorithm
        p_opt = dadi.Inference.optimize_log_fmin(start, fs, dadi_n_epoch_extrap, lower_bound=lower, upper_bound=upper,
                                                 pts=DADI_GRID_PTS, verbose=50, output_file=log_file.path,
                                                 full_output=True)

        # fit the optimised model
        model = dadi_n_epoch_extrap(p_opt[0], fs.sample_sizes, DADI_GRID_PTS)

        # calculate theta given the model
        theta = dadi.Inference.optimal_sfs_scaling(model, fs)

        # save the relevant information
        best = {'epoch': self.epoch, 'n': self.n, 'lnL': -p_opt[1], 'params': p_opt[0], 'theta': theta}

        # save the results by pickling them in a file
        with pkl_file.open('w') as fout:
            pickle.dump(best, fout)


class DadiEpochMaximumLikelihood(utils.PipelineTask):
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
        for n in range(1, int(DADI_REPLICATES) + 1):
            yield DadiEpochOptimizeParams(self.species, self.population, self.folded, self.epoch, n)

    def output(self):
        return luigi.LocalTarget('data/dadi/{}/{}-maxlnL.pkl'.format(self.basename, self.basename))

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


class DadiEpochBestModel(utils.PipelineTask):
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
        yield EasySFS(self.species, self.population, self.folded)

        for epoch in range(1, DADI_EPOCHS + 1):
            yield DadiEpochMaximumLikelihood(self.species, self.population, self.folded, epoch)

    def output(self):
        return [luigi.LocalTarget('data/dadi/{}.{}'.format(self.basename, ext)) for ext in ['pkl', 'pdf']]

    def run(self):
        # unpack the inputs/outputs
        (sfs_file, _), epoch_pkls = self.input()[0], self.input()[1:]
        pkl_out, pdf_out = self.output()

        epochs = []

        # load the pickled max lnL params from each of the epoch models
        for pkl_file in epoch_pkls:
            with pkl_file.open('r') as fin:
                epochs.append(pickle.load(fin))

        # calculate the AIC = 2 * num_params - 2 * lnL
        for epoch in epochs:
            epoch['aic'] = (2 * len(epoch['params'])) - (2 * epoch['lnL'])

        # get the min AIC
        min_aic = min(e['aic'] for e in epochs)

        # compute the relative likelihood of each model = exp((AICmin − AICi)/2)
        for epoch in epochs:
            epoch['relL'] = math.exp((min_aic - epoch['aic']) / 2)

        # reverse sort by relative likelihood
        epochs.sort(key=lambda x: x['relL'], reverse=True)

        # reject modelling if 2nd best model has a relative likelihood greater than acceptable
        if epochs[1]['relL'] > DADI_MAX_RELATIVE_LL:
            raise Exception('ERROR: Cannot reject second best model based on relative likelihood.\n' + str(epochs[:2]))

        # load the frequency spectrum
        fs = dadi.Spectrum.from_file(sfs_file.path)

        # make the extrapolating version of our demographic model function.
        dadi_n_epoch_extrap = dadi.Numerics.make_extrap_log_func(dadi_n_epoch)

        # fit the optimised model
        model = dadi_n_epoch_extrap(epochs[0]['params'], fs.sample_sizes, DADI_GRID_PTS)

        # plot the figure
        fig = plt.figure(1)
        dadi.Plotting.plot_1d_comp_multinom(model, fs, fig_num=1, plot_masked=True)
        fig.savefig(pdf_out.path)
        plt.close(fig)

        # save the results by pickling them in a file
        with pkl_out.open('w') as fout:
            pickle.dump(epochs, fout)


class CountCallableSites(utils.PipelineTask):
    """
    Count the number of callable sites, as dadi needs this number to estimate the ancestral population size from theta.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        for chrom in self.chromosomes:
            yield PolarizeVCF(self.species, chrom)

    def output(self):
        return luigi.LocalTarget('data/dadi/{}.L'.format(self.basename))

    def run(self):
        total = 0

        # count all unique sites
        for vcf_file in self.input():
            size = utils.run_cmd(["bcftools query -f '%CHROM %POS\\n' {} | uniq | wc -l"
                                 .format(vcf_file.path)], shell=True)
            total += int(size)

        with self.output().open('w') as fout:
            fout.write(str(total))


class DadiDemography(utils.PipelineTask):
    """
    Convert the best fitting dadi model into a demography file for `selection`

    :type species: str
    :type population: str
    :type folded: bool
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    folded = luigi.BoolParameter(default=False)

    def requires(self):
        yield DadiEpochBestModel(self.species, self.population, self.folded)
        yield CountCallableSites(self.species)

    def output(self):
        return [luigi.LocalTarget('data/dadi/{}.{}'.format(self.basename, ext)) for ext in ['pop', 'nref']]

    def run(self):
        # unpack the inputs/outputs
        (pkl_file, _), size_file = self.input()
        pop_file, nfef_file = self.output()

        # load the best epoch model
        with pkl_file.open('r') as fin:
            best = pickle.load(fin)[0]

        # unpack the values
        theta, params, epoch = best['theta'], list(best['params']), best['epoch']

        # get the mutation rate
        rate = MUTATION_RATE[self.species]

        # get the count of all callable sites
        length = int(size_file.open().read())

        # dadi scales population size by 2*Nref, where Nref is the size of the most ancient population
        # in dadi, θ = 4*Nref*µ, where µ = number of mutations per generation, so to solve for Nref
        nref = theta / (4 * rate * length)

        # save the Nref, so we can interpret the modelling results
        with nfef_file.open('w') as fout:
            fout.write(str(int(nref)))

        # unpack the dadi model params
        sizes, times = params[:best['epoch']], params[epoch:]

        # reverse the param order (newest first, oldest last)
        sizes.reverse()
        times.reverse()

        # save the demography file
        with pop_file.open('w') as fout:
            age = 0.0
            for i in range(epoch):
                # dadi times are durations for each epoch, so we need to make them cumulative
                age += times[i]

                # save the data
                fout.write('{:.4f}\t0.0\t-{}\n'.format(sizes[i], age))

            # add the infinite time-point
            fout.write('1.0\t0\t-Inf\n')


class DadiPipeline(utils.PipelineWrapperTask):
    """
    Find the best fitting of 5 sequential epoch ∂a∂i models (i.e. 1 epoch, 2 epochs, etc.).

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        for pop in self.populations:
            yield DadiDemography(self.species, pop)


if __name__ == '__main__':
    luigi.run()
