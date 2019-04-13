#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import luigi
import os
import unicodecsv as csv

# import my custom modules
from pipeline_consts import *
from pipeline_utils import PipelineTask, db_conn, run_cmd, trim_ext


# the population history is either: constant, or a fully specified complex demography
MCMC_POP_CONST = 'const'
MCMC_POP_DEMOG = 'demog'

# number of MCMC cycles to run
MCMC_CYCLES = int(5e7)

# frequency of sampling from the posterior
MCMC_SAMPLE_FREQ = 10000

# fraction of MCMC cycles to discard as burn in
MCMC_BURN_IN = 0.2

# frequency of printing output to the screen
MCMC_PRINT = 1000

# fraction of the allele frequency to update during a trajectory update move
MCMC_FRACTION = 20

# derived allele is: 1=dominant, 0=recessive, 0.5=additive
MCMC_DOMINANCE = 0.5

# random number seed
MCMC_RANDOM_SEED = 234395

# minimum number of time bins
MCMC_MIN_BINS = 3


class SelectionInputFile(PipelineTask):
    """
    Generate the 4-column sample input file for `selection` (Schraiber et al., 2016)

    See https://github.com/Schraiber/selection

    The sample input is a 4-column, white-space-separated file. Each line corresponds to a single sample, which can be
    just one individual or many individuals from approximately the same time period pooled together.

    For each sample, each column is
    1. the number of derived alleles
    2. the sample size (in haploid genomes)
    3. the most ancient end of the possible age of the sample (i.e. the oldest it could be)
    4. the most recent end of the possible age of the sample (i.e. the youngest it could be)
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    modsnp_id = luigi.IntParameter()

    def output(self):
        return luigi.LocalTarget("selection/{}.input".format(self.basename))

    def run(self):

        dbc = db_conn(self.species)

        gen_time = GENERATION_TIME[self.species]
        pop_size = POPULATION_SIZE[self.species][self.population]

        # resolve pseudo-population DOM2WLD
        pop = self.population if self.population != 'DOM2WLD' else "DOM2', 'WILD"

        bins = dbc.get_records_sql("""
            # get the ancient frequencies in each bin
            SELECT SUM(sr.base = ms.derived) AS derived_count,
                   COUNT(sr.id) AS sample_size,
                   -(age - (age % 500) + 500) / (2 * {pop_size} * {gen_time}) AS bin_high,
                   -(age - (age % 500) + 1) / (2 * {pop_size} * {gen_time}) AS bin_low
              FROM modern_snps ms
              JOIN sample_reads sr
                ON sr.chrom = ms.chrom
               AND sr.site = ms.site
               AND sr.called = 1
              JOIN samples s
                ON s.id = sr.sample_id
             WHERE ms.id = {modsnp_id}
               AND s.age IS NOT NULL
               AND s.status IN ('{population}')
          GROUP BY bin_high

             UNION

            # add the modern frequency 
            SELECT derived_count, ancestral_count + derived_count, 0, 0
              FROM modern_snps ms
             WHERE ms.id = {modsnp_id}

          ORDER BY bin_high
               """.format(modsnp_id=self.modsnp_id, population=pop, gen_time=gen_time,
                          pop_size=pop_size), key=None)

        if len(bins) < MCMC_MIN_BINS:
            raise RuntimeError('ERROR: Insufficient time bins to run `selection` (n={})'.format(len(bins)))

        # write the sample input file
        with open(self.output().path, "wb") as tsv_file:

            fields = ['derived_count', 'sample_size', 'bin_high', 'bin_low']
            writer = csv.DictWriter(tsv_file, fieldnames=fields, delimiter='\t')

            # write the data to disk
            for b in bins:
                writer.writerow(b)


class SelectionRunMCMC(PipelineTask):
    """
    Run `selection` for the given SNP.
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    modsnp_id = luigi.IntParameter()
    pop_hist = luigi.Parameter()
    mcmc_cycles = luigi.IntParameter()
    mcmc_freq = luigi.IntParameter()

    priority = PRIORITY_HIGH

    def requires(self):
        return SelectionInputFile(self.species, self.population, self.modsnp_id)

    def output(self):
        return [luigi.LocalTarget("selection/{}.{}".format(self.basename, ext))
                for ext in ['log', 'param', 'time', 'traj']]

    def run(self):

        # compose the input and output file paths
        input_file = self.input().path
        log_file = self.output()[0].path
        output_prefix = trim_ext(log_file)
        
        # get path of the population history file
        pop_hist = 'data/selection/{}-{}-{}.pop'.format(self.species, self.population, self.pop_hist)

        try:
            # run `selection`
            run_cmd(['sr',
                     '-D', input_file,        # path to data input file
                     '-P', pop_hist,          # path to population size history file
                     '-o', output_prefix,     # output file prefix
                     '-a',                    # flag to infer allele age
                     '-h', MCMC_DOMINANCE,    # derived allele is: 1=dominant, 0=recessive, 0.5=additive
                     '-n', self.mcmc_cycles,  # number of MCMC cycles to run
                     '-s', self.mcmc_freq,    # frequency of sampling from the posterior
                     '-f', MCMC_PRINT,        # frequency of printing output to the screen
                     '-F', MCMC_FRACTION,     # fraction of the allele frequency to update during a trajectory move
                     '-e', MCMC_RANDOM_SEED,  # random number seed
                     ], stdout=log_file)

        except RuntimeError as e:
            # delete the unfinished *.time and *.traj files
            self.output()[2].remove()
            self.output()[3].remove()

            raise RuntimeError(e)

        # TODO measure ESS and enforce threshold
        # https://www.rdocumentation.org/packages/LaplacesDemon/versions/16.1.0/topics/ESS
        # https://cran.r-project.org/web/packages/coda/index.html


class SelectionPlot(PipelineTask):
    """
    Plot the allele trajectory.

    :type mcmc_cycles: int
    :type mcmc_freq: int
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    modsnp_id = luigi.IntParameter()
    pop_hist = luigi.Parameter(default=MCMC_POP_CONST)
    mcmc_cycles = luigi.IntParameter(default=MCMC_CYCLES)
    mcmc_freq = luigi.IntParameter(default=MCMC_SAMPLE_FREQ)

    priority = PRIORITY_MAX
    resources = {'cpu-cores': CPU_CORES_ONE, 'ram-gb': 64}

    def requires(self):
        return SelectionRunMCMC(*self.all_params())

    def output(self):
        return luigi.LocalTarget("pdf/{}.pdf".format(self.basename))

    def run(self):

        # compose the input and output file paths
        input_file = "selection/{}-{}-{}.input".format(self.species, self.population, self.modsnp_id)
        output_prefix = trim_ext(self.input()[0].path)

        gen_time = GENERATION_TIME[self.species]
        pop_size = POPULATION_SIZE[self.species][self.population]
        burn_in = self.mcmc_cycles / self.mcmc_freq * MCMC_BURN_IN

        try:
            # plot the allele trajectory
            run_cmd(['Rscript', 'rscript/plot-selection.R', input_file, output_prefix, gen_time, pop_size, burn_in])

            # TODO plot the strength of selection (a1 and a2)

        except RuntimeError as e:
            # delete the broken PDF
            if os.path.isfile(self.output().path):
                self.output().remove()

            raise RuntimeError(e)


class SelectionHorseGWAS(luigi.WrapperTask):
    """
    Run `selection` on all the direct GWAS hits for horses.
    """

    def requires(self):

        dbc = db_conn('horse')

        # run for DOM2 and DOM2 + WLD
        pops = ['DOM2', 'DOM2WLD']

        # TODO WTF? why is the qtl_snps join dropping 417 rows?
        # get the modsnp_id for every GWAS hit
        modsnps = dbc.get_records_sql("""
            SELECT DISTINCT ms.id
              FROM qtls q
              JOIN ensembl_variants ev 
                ON ev.rsnumber = q.peak
              JOIN modern_snps ms
                ON ms.variant_id = ev.id
              JOIN qtl_snps qs
                ON qs.qtl_id = q.id
               AND qs.modsnp_id = ms.id 
             WHERE q.associationType = 'Association'
               AND q.valid = 1""")

        # TODO why do these 5 fail?
        bad = [5860440, 4993470, 7833738, 14185088, 9410404]

        for pop in pops:
            for modsnp_id in modsnps:
                if modsnp_id not in bad:
                    yield SelectionPlot('horse', pop, modsnp_id, MCMC_POP_CONST)


class SelectionHorseGWASFlankingSNPs(luigi.WrapperTask):
    """
    Run `selection` on all the direct GWAS hits for horses.
    """

    def requires(self):

        dbc = db_conn('horse')

        # get the modsnp_id for every GWAS hit
        modsnps = dbc.get_records_sql("""
            SELECT DISTINCT qs.modsnp_id AS id
              FROM qtls q
              JOIN qtl_snps qs
                ON qs.qtl_id = q.id
               AND qs.best IS NOT NULL
             WHERE q.associationType = 'Association'
               AND q.valid = 1""")

        for modsnp_id in modsnps:
            yield SelectionPlot('horse', 'DOM2WLD', modsnp_id, MCMC_POP_CONST)


class SelectionHorseTest(luigi.WrapperTask):
    """
    Run `selection` on all a sub-set of SNPs to test MCMC params.
    """

    def requires(self):

        # run for DOM2 and DOM2 + WLD
        pops = ['DOM2', 'DOM2WLD']

        # hand-picked set of SNPs
        modsnps = [5049478, 5102390]

        # run for both constant pop and full demography
        histories = [MCMC_POP_CONST,
                     MCMC_POP_DEMOG]

        # try a few different options
        params = [(1e6, 1e2),
                  (1e6, 1e3),
                  (1e7, 1e3),
                  (1e7, 1e4),
                  (1e8, 1e4)]

        for pop in pops:
            for modsnp_id in modsnps:
                for pop_hist in histories:
                    for cycles, freq in params:
                        yield SelectionPlot('horse', pop, modsnp_id, pop_hist, int(cycles), int(freq))


if __name__ == '__main__':
    luigi.run()
