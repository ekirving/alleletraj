#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
import os

# third party modules
import luigi
import unicodecsv as csv

# local modules
from alleletraj import utils
from alleletraj.ancient.snps import AncientSNPsPipeline
from alleletraj.const import GENERATION_TIME
from alleletraj.modern.demog import DadiDemography
from alleletraj.qtl.analyse import AnalyseQTLsPipeline

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

# derived allele is
MCMC_RECESSIVE = 0
MCMC_ADDITIVE = 0.5
MCMC_DOMINANT = 1

# random number seed
MCMC_RANDOM_SEED = 234395

# minimum number of time bins
MCMC_MIN_BINS = 3


class SelectionInputFile(utils.PipelineTask):
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

    def requires(self):
        yield DadiDemography(self.species, self.population)
        yield AncientSNPsPipeline(self.species)

    def output(self):
        return luigi.LocalTarget('data/selection/{}.input'.format(self.basename))

    def run(self):
        # unpack the inputs
        (_, nref_file), _ = self.input()

        dbc = self.db_conn()

        gen_time = GENERATION_TIME[self.species]

        # get the Nref population size
        with nref_file.open() as fin:
            pop_size = int(fin.read())

        # resolve pseudo-population DOM2WLD
        pop_sql = self.population if self.population != 'DOM2WLD' else "DOM2', 'WILD"

        # noinspection SqlResolve
        bins = dbc.get_records_sql("""
            # get the ancient frequencies in each bin
            SELECT SUM(sr.base = ms.derived) AS derived_count,
                   COUNT(sr.id) AS sample_size,
                   -(age - (age % 500) + 500) / (2 * {pop_size} * {gen_time}) AS bin_high,
                   -(age - (age % 500) + 1) / (2 * {pop_size} * {gen_time}) AS bin_low
              FROM modern_snps ms
              JOIN ancient_sample_reads sr
                ON sr.chrom = ms.chrom
               AND sr.site = ms.site
              JOIN ancient_samples s
                ON s.id = sr.sample_id
             WHERE ms.id = {modsnp_id}
               AND s.age IS NOT NULL
               AND s.population IN ('{pop_sql}')
          GROUP BY bin_high

             UNION

            # add the modern frequency 
            SELECT derived_count, ancestral_count + derived_count, 0, 0
              FROM modern_snps ms
             WHERE ms.id = {modsnp_id}

          ORDER BY bin_high
               """.format(modsnp_id=self.modsnp_id, pop_sql=pop_sql, gen_time=gen_time,
                          pop_size=pop_size), key=None)

        if len(bins) < MCMC_MIN_BINS:
            raise RuntimeError('ERROR: Insufficient time bins to run `selection` (n={})'.format(len(bins)))

        # write the sample input file
        with self.output().open('wb') as tsv_file:

            fields = ['derived_count', 'sample_size', 'bin_high', 'bin_low']
            writer = csv.DictWriter(tsv_file, fieldnames=fields, delimiter='\t')

            # write the data to disk
            for b in bins:
                writer.writerow(b)


class SelectionRunMCMC(utils.PipelineTask):
    """
    Run `selection` for the given SNP.
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    modsnp_id = luigi.IntParameter()
    pop_hist = luigi.Parameter()
    mcmc_cycles = luigi.IntParameter()
    mcmc_freq = luigi.IntParameter()

    def requires(self):
        return SelectionInputFile(self.species, self.population, self.modsnp_id)

    def output(self):
        return [luigi.LocalTarget('data/selection/{}.{}'.format(self.basename, ext))
                for ext in ['log', 'param', 'time', 'traj']]

    def run(self):

        # compose the input and output file paths
        input_file = self.input().path
        log_file, param_file, time_file, traj_file = self.output()
        output_prefix = utils.trim_ext(log_file.path)

        # get path of the population history file
        pop_hist = 'data/selection/{}-{}-{}.pop'.format(self.species, self.population, self.pop_hist)

        try:
            with log_file.open('w') as fout:
                # run `selection`
                cmd = ['sr',
                       '-D', input_file,        # path to data input file
                       '-P', pop_hist,          # path to population size history file
                       '-o', output_prefix,     # output file prefix
                       '-a',                    # flag to infer allele age
                       '-A',                    # ascertainment flag
                       '-h', MCMC_ADDITIVE,     # assume derived allele is additive
                       '-n', self.mcmc_cycles,  # number of MCMC cycles to run
                       '-s', self.mcmc_freq,    # frequency of sampling from the posterior
                       '-f', MCMC_PRINT,        # frequency of printing output to the screen
                       '-F', MCMC_FRACTION,     # fraction of the allele frequency to update during a trajectory move
                       '-e', MCMC_RANDOM_SEED]  # random number seed

                utils.run_cmd(cmd, stdout=fout)

        except RuntimeError as e:
            # delete the unfinished *.time and *.traj files
            time_file.remove()
            traj_file.remove()

            raise RuntimeError(e)

        # TODO gzip the output files
        # TODO measure ESS and enforce threshold
        # https://www.rdocumentation.org/packages/LaplacesDemon/versions/16.1.0/topics/ESS
        # https://cran.r-project.org/web/packages/coda/index.html


class SelectionPlot(utils.PipelineTask):
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

    resources = {'cpu-cores': 1, 'ram-gb': 64}

    def requires(self):
        yield DadiDemography(self.species, self.population)
        yield SelectionRunMCMC(self.species, self.population, self.modsnp_id, self.pop_hist, self.mcmc_cycles,
                               self.mcmc_freq)

    def output(self):
        return luigi.LocalTarget('data/pdf/{}.pdf'.format(self.basename))

    def run(self):
        # unpack the inputs
        (_, nref_file), _ = self.input()

        # compose the input and output file paths
        input_file = 'data/selection/{}-{}-{}.input'.format(self.species, self.population, self.modsnp_id)
        output_prefix = utils.trim_ext(self.input()[0].path)

        gen_time = GENERATION_TIME[self.species]

        # get the Nref population size
        with nref_file.open() as fin:
            pop_size = int(fin.read())

        burn_in = (self.mcmc_cycles / self.mcmc_freq) * MCMC_BURN_IN

        try:
            # plot the allele trajectory
            utils.run_cmd(['Rscript',
                           'rscript/plot-selection.R',
                           input_file,
                           output_prefix,
                           gen_time,
                           pop_size,
                           burn_in])

        except RuntimeError as e:
            # delete the broken PDF
            if os.path.isfile(self.output().path):
                self.output().remove()

            raise RuntimeError(e)


class SelectionGWASSNPs(utils.PipelineWrapperTask):
    """
    Run `selection` on all the direct GWAS hits.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):

        dbc = self.db_conn()

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

        for pop in self.modern_pops:
            for modsnp_id in modsnps:
                yield SelectionPlot(self.species, pop, modsnp_id, MCMC_POP_CONST)


class SelectionBestQTLSNPs(utils.PipelineWrapperTask):
    """
    Run `selection` on all the 'best' QTL SNPs

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        # mark the best SNPs
        yield AnalyseQTLsPipeline(self.species)

        dbc = self.db_conn()

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
            yield SelectionPlot(self.species, 'DOM2WLD', modsnp_id, MCMC_POP_CONST)


class SelectionExportSLURM(utils.PipelineWrapperTask):
    """
    Export a script for batch running `selection` on the HPC, using SLURM.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        # TODO implement this
        pass


if __name__ == '__main__':
    luigi.run()
