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

# number of replicate chains
from alleletraj.qtl.load import MIN_DAF

MCMC_REPLICATES = 4

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

# genetic model
MODEL_RECESSIVE = 0
MODEL_ADDITIVE = 0.5
MODEL_DOMINANT = 1

# minimum number of time bins
MCMC_MIN_BINS = 3


class SelectionInputFile(utils.DatabaseTask):
    """
    Generate the 4-column sample input file for `selection`.

    See https://github.com/Schraiber/selection

    The sample input is a 4-column, white-space-separated file. Each line corresponds to a single sample, which can be
    just one individual or many individuals from approximately the same time period pooled together.

    For each sample, each column is
    1. the number of derived alleles
    2. the sample size (in haploid genomes)
    3. the most ancient end of the possible age of the sample (i.e. the oldest it could be)
    4. the most recent end of the possible age of the sample (i.e. the youngest it could be)

    :type species: str
    :type population: str
    :type modsnp: int
    :type mispolar: bool
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    modsnp = luigi.IntParameter()
    mispolar = luigi.BoolParameter()

    def requires(self):
        yield DadiDemography(self.species, self.population)
        yield AncientSNPsPipeline(self.species)

    def output(self):
        return luigi.LocalTarget('data/selection/{}.input'.format(self.basename))

    def run(self):
        # unpack the inputs
        (_, nref_file), _ = self.input()

        # get the Nref population size
        with nref_file.open() as fin:
            pop_size = int(fin.read())

        # NOTE some SNPs may be mispolarised, so we switch the derived/ancestral alleles in those cases
        derived = 'derived' if not self.mispolar else 'ancestral'

        gen_time = GENERATION_TIME[self.species]

        # noinspection SqlResolve
        bins = self.dbc.get_records_sql("""
            # get the ancient frequencies in each bin
            SELECT SUM(sr.base = ms.{derived}) AS derived_count,
                   COUNT(sr.id) AS sample_size,
                   -(sb.max / (2 * {pop_size} * {gen_time})) AS max,
                   -(sb.min / (2 * {pop_size} * {gen_time})) AS min                   
              FROM modern_snps ms
              JOIN sample_reads sr
                ON sr.chrom = ms.chrom
               AND sr.site = ms.site
              JOIN samples s
                ON s.id = sr.sample_id
              JOIN sample_bins sb
                ON sb.id = s.bin_id  
             WHERE ms.id = {modsnp}
               AND s.population = '{population}'
          GROUP BY sb.id

             UNION

            # add the modern frequency 
            SELECT {derived}_count, ancestral_count + derived_count, 0, 0
              FROM modern_snps ms
              JOIN modern_snp_daf msd
                ON msd.modsnp_id = ms.id
             WHERE ms.id = {modsnp}
               AND msd.population = '{population}'

          ORDER BY max
               """.format(derived=derived, modsnp=self.modsnp, population=self.population, pop_size=pop_size,
                          gen_time=gen_time), key=None)

        if len(bins) < MCMC_MIN_BINS:
            raise RuntimeError('ERROR: Insufficient time bins to run `selection` (n={})'.format(len(bins)))

        # write the sample input file
        with self.output().open('w') as tsv_file:
            fields = ['derived_count', 'sample_size', 'max', 'min']
            writer = csv.DictWriter(tsv_file, fieldnames=fields, delimiter='\t')

            # write the data to disk
            for b in bins:
                writer.writerow(b)


class SelectionRunMCMC(utils.PipelineTask):
    """
    Run `selection` for the given SNP.

    :type species: str
    :type population: str
    :type modsnp: int
    :type chain: int
    :type n: int
    :type s: int
    :type h: float
    :type mispolar: bool
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    modsnp = luigi.IntParameter()
    chain = luigi.IntParameter()
    n = luigi.IntParameter()
    s = luigi.IntParameter()
    h = luigi.FloatParameter()
    mispolar = luigi.BoolParameter()

    def requires(self):
        yield DadiDemography(self.species, self.population)
        yield SelectionInputFile(self.species, self.population, self.modsnp, self.mispolar)

    def output(self):
        return [luigi.LocalTarget('data/selection/{}.{}'.format(self.basename, ext))
                for ext in ['log', 'param', 'time.gz', 'traj.gz']]

    def run(self):
        # compose the input and output file paths
        (pop_file, _), input_file = self.input()
        log_file, param_file, time_file, traj_file = self.output()

        output_prefix = utils.trim_ext(log_file.path)

        # make a deterministic random seed
        seed = int('{}{}'.format(self.modsnp, self.chain))

        try:
            with log_file.open('w') as fout:
                # run `selection`
                cmd = ['sr',
                       '-D', input_file.path,   # path to data input file
                       '-P', pop_file.path,     # path to population size history file
                       '-o', output_prefix,     # output file prefix
                       '-a',                    # flag to infer allele age
                       '-A', MIN_DAF,           # ascertainment in modern individuals
                       '-n', self.n,            # number of MCMC cycles to run
                       '-s', self.s,            # frequency of sampling from the posterior
                       '-h', self.h,            # genetic model (additive, recessive, dominant)
                       '-f', MCMC_PRINT,        # frequency of printing output to the screen
                       '-F', MCMC_FRACTION,     # fraction of the allele frequency to update during a trajectory move
                       '-e', seed]              # random number seed

                utils.run_cmd(cmd, stdout=fout)

        except RuntimeError as e:
            # delete the unfinished *.time and *.traj files
            os.remove(utils.trim_ext(time_file.path))
            os.remove(utils.trim_ext(traj_file.path))

            raise RuntimeError(e)

        else:
            # gzip the output files
            utils.run_cmd(['gzip {}'.format(time_file.path)], shell=True)
            utils.run_cmd(['gzip {}'.format(traj_file.path)], shell=True)

            # TODO measure ESS and enforce threshold
            # https://www.rdocumentation.org/packages/LaplacesDemon/versions/16.1.0/topics/ESS
            # https://cran.r-project.org/web/packages/coda/index.html


class SelectionPlot(utils.PipelineTask):
    """
    Plot the allele trajectory.

    :type species: str
    :type population: str
    :type modsnp: int
    :type chain: int
    :type n: int
    :type s: int
    :type h: float
    :type mispolar: bool
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    modsnp = luigi.IntParameter()
    chain = luigi.IntParameter()
    n = luigi.IntParameter(default=MCMC_CYCLES)
    s = luigi.IntParameter(default=MCMC_SAMPLE_FREQ)
    h = luigi.FloatParameter(default=MODEL_ADDITIVE)
    mispolar = luigi.BoolParameter(default=False)

    resources = {'cpu-cores': 1, 'ram-gb': 64}

    def requires(self):
        yield DadiDemography(self.species, self.population)
        yield SelectionRunMCMC(self.species, self.population, self.modsnp, self.chain, self.n, self.s, self.h,
                               self.mispolar)

    def output(self):
        return luigi.LocalTarget('data/pdf/{}.pdf'.format(self.basename))

    def run(self):
        # unpack the inputs
        (_, nref_file), _ = self.input()

        # compose the input and output file paths
        input_file = 'data/selection/{}-{}-{}.input'.format(self.species, self.population, self.modsnp)
        output_prefix = utils.trim_ext(self.input()[0].path)

        gen_time = GENERATION_TIME[self.species]

        # get the Nref population size
        with nref_file.open() as fin:
            pop_size = int(fin.read())

        burn_in = (self.n / self.s) * MCMC_BURN_IN

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

        # get the modsnp id for every GWAS hit
        modsnps = self.dbc.get_records_sql("""
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

        for pop in self.list_populations(modern=True):
            for modsnp in modsnps:
                for chain in range(MCMC_REPLICATES):
                    yield SelectionPlot(self.species, pop, modsnp, chain)


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
