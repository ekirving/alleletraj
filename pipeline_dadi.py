#!/usr/bin/env python
# -*- coding: utf-8 -*-

import dadi
import logging
import luigi
import pickle

from datetime import datetime

# import my custom modules
from pipeline_consts import *
from pipeline_snp_call import PolarizeVCF, SubsetSNPsVCF
from pipeline_utils import PipelineTask, run_cmd, LogBuffer


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
    polarised = luigi.BoolParameter()

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
            'fold': '--unfolded' if self.polarised else ''
        }

        # pipe 'yes' into easySFS to get past the interactive prompt which complains about excluded samples
        cmd = "echo 'yes' | easySFS.py -a -i {vcf} -p {pops} -o sfs/{out} --proj {proj} {fold}".format(**params)

        run_cmd([cmd], shell=True)


class DadiOptimizeParams(PipelineTask):
    """
    Optimise the log likelihood of the model parameters for the given SFS.

    # TODO add param types
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    polarised = luigi.BoolParameter()

    model = luigi.Parameter()
    scenario = luigi.Parameter()
    grid_size = luigi.ListParameter()
    upper_bound = luigi.ListParameter()
    lower_bound = luigi.ListParameter()
    param_start = luigi.ListParameter(significant=False)

    n = luigi.IntParameter()

    def requires(self):
        return EasySFS(self.species, self.population, self.polarised)

    def output(self):
        return luigi.LocalTarget("sfs/{}.opt".format(self.basename))

    def run(self):

        # load the frequency spectrum
        fs = dadi.Spectrum.from_file(self.input().path)
        ns = fs.sample_sizes

        # get the demographic model to test
        func = getattr(dadi.Demographics2D, self.model)

        # Make the extrapolating version of our demographic model function.
        func_ex = dadi.Numerics.make_extrap_log_func(func)

        # keep a list of the optimal params
        p_best = []

        # buffer the logs created by dadi so we can inspect them
        log_buffer = LogBuffer()
        logger = logging.getLogger('Inference')
        logger.addHandler(logging.StreamHandler(log_buffer))

        # run the optimisation many times
        for i in range(0, DADI_MAX_ITER):

            # Perturb our parameters before optimization. This does so by taking each
            # parameter a up to a factor of two up or down.
            p_perturb = dadi.Misc.perturb_params(self.param_start,
                                                 fold=1,
                                                 upper_bound=self.upper_bound,
                                                 lower_bound=self.lower_bound)

            print('Started optimization: {:<12} | {:<8} | {:<8} | {:<9} | {:<15} | n={:>3} | i={:>3} |'.format(
                self.group, self.pop1, self.pop2, self.model, self.scenario, self.n, i))

            start = datetime.now()

            # do the optimization...
            p_opt = dadi.Inference.optimize_log(p_perturb, fs, func_ex, self.grid_size,
                                                lower_bound=self.lower_bound,
                                                upper_bound=self.upper_bound,
                                                fixed_params=self.fixed_params,
                                                verbose=20,
                                                maxiter=DADI_MAX_ITER)

            end = datetime.now()
            diff = (end - start).total_seconds() / 60

            print('Finshed optimization: {:<12} | {:<8} | {:<8} | {:<9} | {:<15} | n={:>3} | i={:>3} | t={:>3.1f} mins'.format(
                self.group, self.pop1, self.pop2, self.model, self.scenario, self.n, i, diff))

            # reset the log buffer
            log_buffer.log = []

            # Calculate the best-fit model AFS.
            model = func_ex(p_opt, ns, self.grid_size)

            # Likelihood of the data given the model AFS.
            ll_model = dadi.Inference.ll_multinom(model, fs)

            # The optimal value of theta given the model.
            theta = dadi.Inference.optimal_sfs_scaling(model, fs)

            # get the buffered warnings
            warnings = " ".join(log_buffer.log)

            # we only care about non-masked data
            if "Model is masked" not in warnings:

                print('Maximum log composite likelihood: {0}'.format(ll_model))
                print('Optimal value of theta: {0}'.format(theta))

                # record the best fitting params for this run
                p_best.append([ll_model, theta] + list(p_opt))

                # sort the params (largest ll_model first)
                p_best.sort(reverse=True)

            try:
                # run the optimisation again, starting from the best fit we've seen so far
                self.param_start = p_best[0][2:]
            except IndexError:
                # otherwise, generate a set of new random starting params (so we don't get stuck in bad param space)
                self.param_start = random_params(self.lower_bound, self.upper_bound, self.fixed_params)

        # if we've run the iteration DADI_MAX_ITER times and not found any non-masked params then we've failed
        if not p_best:
            raise Exception("{} | {} | {} | n={}: FAILED to find any non-masked params".format(self.group,
                                                                                               self.model,
                                                                                               self.scenario,
                                                                                               self.n))
        # save the list of optimal params by pickeling them in a file
        with self.output().open('w') as fout:
            pickle.dump(p_best, fout)


class DadiModelDemography(luigi.WrapperTask):
    """
    Find the best fitting of 5 sequential epoch models (i.e. 1 epoch, 2 epoch, etc.).
    """

    def requires(self):
        return EasySFS('horse', 'DOM2')


if __name__ == '__main__':
    luigi.run()
