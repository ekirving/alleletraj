#!/usr/bin/env python
# -*- coding: utf-8 -*-

import dadi
import luigi
import pickle
import random
import scipy.stats as st

# import my custom modules
from pipeline_consts import *
from pipeline_snp_call import PolarizeVCF, SubsetSNPsVCF
from pipeline_utils import PipelineTask, run_cmd, dump

# number of sequential epochs to test
DADI_MAX_EPOCHS = 5

# how many independent runs should we do to find global maximum of params
DADI_MAX_REPS = 50  # TODO put back to 50
# DADI_MAX_ITER = 2


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
        return [luigi.LocalTarget("sfs/{}-params.{}".format(self.basename, ext)) for ext in ['pickle', 'log']]

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

        print("start = {}".format(start))

        # optimize log(params) to fit model to data using Nelder-Mead algorithm
        p_opt = dadi.Inference.optimize_log_fmin(start, fs, dadi_epoch_model, lower_bound=lower, upper_bound=upper,
                                                 pts=100, verbose=50, output_file=log_file.path, full_output=True)

        # save the params by pickling them in a file
        with pkl_file.open('w') as fout:
            pickle.dump(p_opt, fout)

        # # calculate the best-fit model AFS
        # model = dadi_epoch_model(p_opt, fs.sample_sizes, self.grid_size)
        #
        # # likelihood of the data given the model AFS
        # ll_model = dadi.Inference.ll_multinom(model, fs)
        #
        # # the optimal value of theta given the model
        # theta = dadi.Inference.optimal_sfs_scaling(model, fs)
        #
        # # record the best fitting params for this run
        # best = [ll_model, theta] + list(p_opt)
        #
        # # save the params by pickling them in a file
        # with pkl_file.open('w') as fout:
        #     pickle.dump(best, fout)


class DadiEpochMaximumLikelihood(PipelineTask):
    """
    Find the model run with the maximum likelihood out of all the replicates.

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
        for n in range(1, DADI_MAX_REPS + 1):
            yield DadiEpochOptimizeParams(self.species, self.population, self.folded, self.epoch, n)

    def output(self):
        return luigi.LocalTarget("sfs/{}-best.{}".format(self.basename, 'pickle'))

    def run(self):

        params = []

        # load all the pickled param values from the replicate runs
        for pkl_file, log_file in self.input():
            with pkl_file.open('r') as fin:
                params.append(pickle.load(fin))

        # TODO find the best one
        print('\n\n\n----------------------\n\n\n')
        dump(params)
        print('\n\n\n----------------------\n\n\n')


class DadiModelDemography(luigi.WrapperTask):
    """
    Find the best fitting of 5 sequential epoch models (i.e. 1 epoch, 2 epoch, etc.).
    """

    def requires(self):
        # run multiple epochs, with multiple independent replicates
        for epoch in range(1, DADI_MAX_EPOCHS + 1):
            yield DadiEpochMaximumLikelihood('horse', 'DOM2', False, epoch)


if __name__ == '__main__':
    luigi.run()
