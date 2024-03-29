#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
import csv
import math

# third party modules
import luigi
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import dadi
import pickle
import random

# local modules
from alleletraj.const import MUTATION_RATE
from vcf import PolarizeVCF, WholeAutosomeSNPsVCF
from alleletraj import utils

# number of sequential epochs to test
DADI_MAX_EPOCHS = 5  # TODO increase this, as 5 was the best model for cattle and it would be preferable to overshoot

# how many independent replicates should we run to find the global maximum params (dadi can get stuck in local maxima)
DADI_REPLICATES = 1000

# number of points to use in the grid
DADI_GRID_PTS = 100

# maximum relative log likelihood to not reject the second best model
DADI_MAX_RELATIVE_LL = 0.10

# is spectrum folded or polarised
DADI_FOLDED = True
DADI_UNFOLDED = False


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


class EasySFSPopFile(utils.DatabaseTask):
    """
    Make sample and population files for the just those samples used to calculate the SFS.

    :type species: str
    :type population: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()

    def output(self):
        return [luigi.LocalTarget('data/sfs/{}.{}'.format(self.basename, ext)) for ext in ['pops', 'spl']]

    def run(self):
        pop_file, spl_file = self.output()

        # get the list of samples to use for the SFS calculation
        samples = self.dbc.get_records('samples', {'population': self.population, 'ancient': 0, 'sfs': 1}, key='name')

        # make both the samples and populations files
        with pop_file.open('w') as pop_fout, spl_file.open('w') as spl_fout:
            for sample in samples:
                pop_fout.write('{}\t{}\n'.format(sample, self.population))
                spl_fout.write('{}\n'.format(sample))


class EasySFS(utils.DatabaseTask):
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
        yield EasySFSPopFile(self.species, self.population)
        yield WholeAutosomeSNPsVCF(self.species)

    def output(self):
        return [luigi.LocalTarget('data/sfs/{}/dadi/{}.{}'.format(self.basename, self.population, ext))
                for ext in ['sfs', 'log']]

    def run(self):
        # unpack the params
        (pop_file, _), vcf_file = self.input()
        sfs_file, log_file = self.output()

        # get the number of samples in the pop file
        count = utils.run_cmd(['wc', '-l', pop_file.path])
        num_samples = int(count.split()[0])

        params = {
            'vcf': vcf_file.path,
            'pops': pop_file.path,
            'out': self.basename,
            'proj': num_samples * 2,  # don't project down
            'fold': '--unfolded' if not self.folded else ''
        }

        # TODO include all populations so we only have to load the VCF once!!
        # NOTE easySFS expects the REF allele to be ancestral, rather than using the INFO/AA field
        # pipe 'yes' into easySFS to get past the interactive prompt which complains about excluded samples
        cmd = "echo 'yes' | easySFS.py -a -f -i {vcf} -p {pops} -o data/sfs/{out} --proj {proj} {fold}".format(**params)

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

        # make a deterministic random seed (helps keep everything easily reproducible)
        seed = int('{}{}'.format(self.epoch, self.n))
        random.seed(seed)

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
        yield EasySFS(self.species, self.population, self.folded)

        for n in range(1, DADI_REPLICATES + 1):
            yield DadiEpochOptimizeParams(self.species, self.population, self.folded, self.epoch, n)

    def output(self):
        yield luigi.LocalTarget('data/dadi/{}/{}-maxlnL.pkl'.format(self.basename, self.basename))
        yield luigi.LocalTarget('data/dadi/{}.params'.format(self.basename))
        yield luigi.LocalTarget('data/dadi/{}.pdf'.format(self.basename))

    def run(self):
        # unpack the params
        (sfs_file, _), pkl_files = self.input()[0], self.input()[1:]
        pkl_out, params_file, pdf_file = self.output()

        params = []

        # load all the pickled param values from the replicate runs
        for pkl_file, _ in pkl_files:
            with pkl_file.open('r') as fin:
                params.append(pickle.load(fin))

        # find the params that produced the highest maximum likelihood
        max_lnl = max(params, key=lambda x: x['lnL'])

        # load the frequency spectrum
        fs = dadi.Spectrum.from_file(sfs_file.path)

        # make the extrapolating version of our demographic model function.
        dadi_n_epoch_extrap = dadi.Numerics.make_extrap_log_func(dadi_n_epoch)

        # fit the optimised model
        model = dadi_n_epoch_extrap(max_lnl['params'], fs.sample_sizes, DADI_GRID_PTS)

        # plot the figure
        fig = plt.figure(1)
        dadi.Plotting.plot_1d_comp_multinom(model, fs, fig_num=1, plot_masked=True)
        fig.savefig(pdf_file.path)
        plt.close(fig)

        # save the pickle and the params
        with pkl_out.open('w') as fout:
            pickle.dump(max_lnl, fout)

        with params_file.open('w') as fout:
            fout.write('{}'.format(max_lnl))


class CountChromSites(utils.DatabaseTask):
    """
    Count the number of callable sites in a chromosome.

    :type species: str
    :type population: str
    :type chrom: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        yield EasySFSPopFile(self.species, self.population)
        yield PolarizeVCF(self.species, self.chrom)

    def output(self):
        return luigi.LocalTarget('data/sfs/{}.size'.format(self.basename))

    def run(self):
        # unpack the params
        (_, spl_file), vcf_file = self.input()
        size_file = self.output()

        # count the unique sites
        cmd = "bcftools view --samples-file {} --exclude-uncalled --exclude-types indels,mnps,bnd,other {} | " \
              "bcftools query --format '%CHROM %POS\\n' | uniq | wc -l".format(spl_file.path, vcf_file.path)

        size = utils.run_cmd([cmd], shell=True)

        with size_file.open('w') as fout:
            fout.write('{}'.format(size))


class CountCallableSites(utils.DatabaseTask):
    """
    Count the total number of callable sites, so dadi cant estimate the ancestral population size from theta.

    :type species: str
    :type population: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()

    def requires(self):
        for chrom in self.autosomes:
            yield CountChromSites(self.species, self.population, chrom)

    def output(self):
        return luigi.LocalTarget('data/dadi/{}.L'.format(self.basename))

    def run(self):
        # unpack the params
        chrom_files = self.input()
        size_file = self.output()

        total = 0

        # sum all the chromosome sizes
        for chrom_file in chrom_files:
            with chrom_file.open() as fin:
                total += int(fin.read())

        with size_file.open('w') as fout:
            fout.write('{}'.format(total))


class DadiEpochDemography(utils.PipelineTask):
    """
    Convert the best fitting dadi model for a given epoch into a demography file for `selection`.

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
        yield DadiEpochMaximumLikelihood(self.species, self.population, self.folded, self.epoch)
        yield CountCallableSites(self.species, self.population)

    def output(self):
        (pkl_file, _, _), _ = self.input()
        yield pkl_file  # pass on the pickle file
        yield luigi.LocalTarget('data/dadi/{}.pop'.format(self.basename))
        yield luigi.LocalTarget('data/dadi/{}.nref'.format(self.basename))

    def run(self):
        # unpack the inputs/outputs
        (pkl_file, _, _), size_file = self.input()
        _, pop_file, nref_file = self.output()

        # load the best params for the epoch
        with pkl_file.open('r') as fin:
            best = pickle.load(fin)

        # unpack the values
        theta, params, epoch = best['theta'], list(best['params']), best['epoch']

        # get the mutation rate
        rate = MUTATION_RATE[self.species]

        # get the count of all callable sites
        length = int(size_file.open().read())

        # in dadi, θ = 4*Nref*µ*L, where µ = mutation rate and L = size of the region used to estimate the SFS
        # so to solve for Nref...
        nref = theta / (4 * rate * length)

        # save the Nref, so we can interpret the modelling results
        with nref_file.open('w') as fout:
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


class DadiBestEpochModel(utils.PipelineTask):
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
        for epoch in range(1, DADI_MAX_EPOCHS + 1):
            yield DadiEpochDemography(self.species, self.population, self.folded, epoch)

    def output(self):
        return [luigi.LocalTarget('data/dadi/{}.{}'.format(self.basename, ext)) for ext in ['pop', 'nref', 'aic']]

    def run(self):
        # unpack the params
        epoch_pkls = self.input()
        pop_out, nref_out, aic_file = self.output()

        epochs = []
        files = {}

        # load the pickled max lnL params from each of the epoch models
        for pkl_file, pop_in, nref_in in epoch_pkls:
            with pkl_file.open('r') as fin:
                epoch = pickle.load(fin)

            epochs.append(epoch)
            files[epoch['epoch']] = (pop_in, nref_in)

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
            print('WARNING: Cannot reject second best model based on relative likelihood.\n' + str(epochs[:2]))

        # copy the files from the the best epoch
        pop_in, nref_in = files[epochs[0]['epoch']]
        pop_in.copy(pop_out.path)
        nref_in.copy(nref_out.path)

        # save the results
        with aic_file.open('w') as tsv_file:
            fields = ['relL', 'aic', 'lnL', 'theta', 'epoch', 'params', 'n']
            writer = csv.DictWriter(tsv_file, fieldnames=fields, delimiter='\t')
            writer.writeheader()

            for epoch in epochs:
                writer.writerow(epoch)


class DadiConstPopModel(utils.PipelineTask):
    """
    Make a constant population size from from the most recent epoch from the best fitting multi-epoch model.

    :type species: str
    :type population: str
    :type folded: bool
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    folded = luigi.BoolParameter()

    def requires(self):
        return DadiBestEpochModel(self.species, self.population, self.folded)

    def output(self):
        return [luigi.LocalTarget('data/dadi/{}-const_pop.{}'.format(self.basename, ext)) for ext in
                ['pop', 'nref', 'aic']]

    def run(self):
        pop_in, nref_in, _ = self.input()
        pop_out, nref_out, aic_out = self.output()

        # get the effective population size of the most recent epoch
        with pop_in.open('r') as pop_fin, nref_in.open('r') as nref_fin:
            nref = float(pop_fin.read().split()[0]) * int(nref_fin.read())

        # make the const pop model
        with pop_out.open('w') as pop_fout, nref_out.open('w') as nref_out:
            pop_fout.write('1.0\t0.0\t-Inf')
            nref_out.write(str(int(nref)))

        # make a blank aic file for symmetry
        with aic_out.open('w') as fout:
            fout.write('')


class DadiBestModel(utils.PipelineWrapperTask):
    """
    Wrapper task to return the best fitting model.

    :type species: str
    :type population: str
    :type folded: bool
    :type const_pop: bool
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    folded = luigi.BoolParameter()
    const_pop = luigi.BoolParameter(default=False)

    def requires(self):
        if self.const_pop:
            return DadiConstPopModel(self.species, self.population, self.folded)
        else:
            return DadiBestEpochModel(self.species, self.population, self.folded)

    def output(self):
        return self.input()


class DadiPipeline(utils.PipelineWrapperTask):
    """
    Find the best fitting of 5 sequential epoch ∂a∂i models (i.e. 1 epoch, 2 epochs, etc.).

    :type species: str
    :type population: str
    """
    species = luigi.Parameter()
    population = luigi.OptionalParameter(default=None)

    def requires(self):

        populations = [self.population] if self.population else self.list_populations(modern=True)

        for pop in populations:
            for folded in [DADI_FOLDED, DADI_UNFOLDED]:
                for const_pop in [True, False]:
                    yield DadiBestModel(self.species, pop, folded, const_pop)


if __name__ == '__main__':
    luigi.run()
