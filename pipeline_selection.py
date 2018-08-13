#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import luigi
import unicodecsv as csv

from pipeline_utils import *


class SelectionInputFile(luigi.Task):
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

    priority = PRIORITY_LOW
    resources = {'cpu-cores': CPU_CORES_ONE}

    def output(self):
        return luigi.LocalTarget("selection/{}-{}-{}.input".format(self.species, self.population, self.modsnp_id))

    def run(self):

        dbc = db_conn(self.species)

        gen_time = GENERATION_TIME[self.species]
        pop_size = POPULATION_SIZE[self.species][self.population]

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
               AND s.status = '{population}'
          GROUP BY bin_high

             UNION

            # add the modern frequency 
            SELECT derived_count, ancestral_count + derived_count, 0, 0
              FROM modern_snps ms
             WHERE ms.id = {modsnp_id}

          ORDER BY bin_high
               """.format(modsnp_id=self.modsnp_id, population=self.population, gen_time=gen_time,
                          pop_size=pop_size), key=None)

        if len(bins) < MCMC_MIN_BINS:
            raise RuntimeError('ERROR: Insufficient time bins to run `selection` (n={})'.format(len(bins)))

        # write the sample input file
        with open(self.output().path, "wb") as tsv_file:

            fields = ['derived_count', 'sample_size', 'bin_high', 'bin_low']
            writer = csv.DictWriter(tsv_file, fieldnames=fields, delimiter='\t')

            # write the data to disk
            for bin in bins:
                writer.writerow(bin)


class SelectionRunMCMC(luigi.Task):
    """
    Run `selection` for the given SNP.
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    modsnp_id = luigi.IntParameter()

    priority = PRIORITY_HIGH
    resources = {'cpu-cores': CPU_CORES_ONE}

    def requires(self):
        return SelectionInputFile(self.species, self.population, self.modsnp_id)

    def output(self):
        return [luigi.LocalTarget("selection/{}-{}-{}.{}".format(self.species, self.population, self.modsnp_id, ext))
                    for ext in ['log', 'param', 'time', 'traj']]

    def run(self):

        # compose the input and output file paths
        input_file = self.input().path
        log_file = self.output()[0].path
        output_prefix = trim_ext(log_file)
        
        # get path of the population history file
        pop_hist = POPULATION_HISTORY[self.species][self.population]

        try:
            # run `selection`
            run_cmd(['sr',
                     '-D', input_file,        # path to data input file
                     '-P', pop_hist,          # path to population size history file
                     '-o', output_prefix,     # output file prefix
                     '-a',                    # flag to infer allele age
                     '-h', MCMC_DOMINANCE,    # derived allele is: 1=dominant, 0=recessive, 0.5=additive
                     '-n', MCMC_CYCLES,       # number of MCMC cycles to run
                     '-f', MCMC_PRINT,        # frequency of printing output to the screen
                     '-s', MCMC_SAMPLE_FREQ,  # frequency of sampling from the posterior
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


class SelectionPlot(luigi.Task):
    """
    Plot the allele trajectory.
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    modsnp_id = luigi.IntParameter()

    priority = PRIORITY_MAX
    resources = {'cpu-cores': CPU_CORES_ONE, 'ram-gb': 64}

    def requires(self):
        return SelectionRunMCMC(self.species, self.population, self.modsnp_id)

    def output(self):
        return luigi.LocalTarget("pdf/{}-{}-{}.pdf".format(self.species, self.population, self.modsnp_id))

    def run(self):

        # compose the input and output file paths
        input_file = "selection/{}-{}-{}.input".format(self.species, self.population, self.modsnp_id)
        output_prefix = trim_ext(input_file)

        gen_time = GENERATION_TIME[self.species]
        pop_size = POPULATION_SIZE[self.species][self.population]

        try:
            # plot the allele trajectory
            run_cmd(['Rscript', 'rscript/plot-selection.R', input_file, output_prefix, gen_time, pop_size, MCMC_BURN_IN])

        except RuntimeError as e:
            # delete the broken PDF
            self.output().remove()

            raise RuntimeError(e)

        # TODO plot the strength of selection (a1 and a2)


class SelectionHorseGWAS(luigi.WrapperTask):
    """
    Run `selection` on all the direct GWAS hits for horses.
    """

    def requires(self):

        dbc = db_conn('horse')

        # get the modsnp_id for every GWAS hit
        modsnps = dbc.get_records_sql("""
            SELECT ms.id
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

        for modsnp_id in modsnps:
            yield SelectionPlot('horse', 'DOM2', modsnp_id)


if __name__ == '__main__':
    luigi.run()
