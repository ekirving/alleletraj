#!/usr/bin/env python
# -*- coding: utf-8 -*-

import dadi
import luigi
import pickle
import random

# import my custom modules
from pipeline_consts import *
from pipeline_snp_call import PolarizeVCF, SubsetSNPsVCF
from pipeline_utils import PipelineTask, run_cmd

# number of sequential epochs to test
DADI_MAX_EPOCHS = 5

# how many iterations to use to optimise params
# DADI_MAX_ITER = 50  # TODO put back to 50
DADI_MAX_ITER = 2


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


class DadiOptimizeParams(PipelineTask):
    """
    Optimise the log likelihood of the model parameters for the given SFS.

    :type species: str
    :type population: str
    :type folded: bool
    :type epochs: int
    :type n: int
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    folded = luigi.BoolParameter()
    epochs = luigi.IntParameter()
    n = luigi.IntParameter()

    def requires(self):
        return EasySFS(self.species, self.population, self.folded)

    def output(self):
        return [luigi.LocalTarget("sfs/{}-optima.{}".format(self.basename, ext)) for ext in ['pickle', 'log']]

    def run(self):

        # unpack the inputs/outputs
        sfs_file = self.input()
        pkl_file, log_file = self.output()

        # load the frequency spectrum
        spec = dadi.Spectrum.from_file(sfs_file.path)

        # set the upper and lower parameter bounding
        lower = [.01] * self.epochs + [0] * self.epochs
        upper = [100] * self.epochs + [5] * self.epochs

        # pick random starting values (bounded by lower/upper)
        start = [random.uniform(lower[i], upper[i]) for i in range(0, self.epochs * 2)]

        print("start = {}".format(start))

        # optimize log(params) to fit model to data using Nelder-Mead algorithm
        best = dadi.Inference.optimize_log_fmin(start, spec, dadi_epoch_model, lower_bound=lower, upper_bound=upper,
                                                pts=100, verbose=50, output_file=log_file.path, full_output=True)

        # save the optimal params by pickling them in a file
        with pkl_file.open('w') as fout:
            pickle.dump(best, fout)


class DadiModelDemography(luigi.WrapperTask):
    """
    Find the best fitting of 5 sequential epoch models (i.e. 1 epoch, 2 epoch, etc.).
    """

    def requires(self):

        # run multiple sequential epochs
        for epochs in range(1, DADI_MAX_EPOCHS + 1):

            # each each epoch, run multiple replicates
            for n in range(1, DADI_MAX_ITER + 1):
                yield DadiOptimizeParams('horse', 'DOM2', False, epochs, n)


if __name__ == '__main__':
    luigi.run()
