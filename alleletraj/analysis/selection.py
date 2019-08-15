#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
import glob
import json
import os
import random
from time import time

# third party modules
import luigi
import unicodecsv as csv

# local modules
from alleletraj import utils
from alleletraj.ancient.snps import AncientSNPsPipeline
from alleletraj.const import GENERATION_TIME
from alleletraj.db.conn import Database
from alleletraj.modern.demog import DadiBestModel, DADI_FOLDED
from alleletraj.qtl.load import MIN_DAF

# number of independent MCMC replicates to run
MCMC_NUM_CHAINS = 2  # TODO put back to 4 for production run

# number of MCMC cycles to run
MCMC_CYCLES = int(1e8)

# number of MCMC cycles to benchmark
MCMC_BENCHMARK = int(1e4)

# fraction of MCMC cycles to discard as burn in
MCMC_BURN_PCT = 0.5

# thinning is analytically unnecessary, but makes the MCMC run much faster (as Josh's method writes directly to disk)
# see https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.2041-210X.2011.00131.x
MCMC_THIN = 100

# frequency of printing output to the log
MCMC_PRINT = 1000

# genetic model
MODEL_RECESSIVE = 0
MODEL_ADDITIVE = 0.5
MODEL_DOMINANT = 1

# minimum number of time bins
MCMC_MIN_BINS = 3

# number of DAF paired neutral SNPs to run for every non-neutral SNP
NEUTRAL_REPLICATES = 5

# max number of models to schedule at one time
BENCHMARK_LIMIT = 10000


def selection_fetch_neutral_snps(species, population, modsnp_id, mispolar=False):
    """
    Fetch the 'neutral' SNPs paired to the modsnp in this population.
    """
    dbc = Database(species)

    params = {
        'population': population,
        'modsnp_id':  modsnp_id,
        'mispolar':   int(mispolar)
    }

    modsnps = dbc.get_records('selection_neutrals', params)

    if not modsnps:
        modsnps = selection_pair_neutral_snps(species, population, modsnp_id, mispolar)

    if len(modsnps) != NEUTRAL_REPLICATES:
        print('WARNING: Insufficient neutral SNPs for `selection` {} (n={})'.format([population, modsnp_id, mispolar],
                                                                                    len(modsnps)))

    return modsnps


def selection_pair_neutral_snps(species, population, modsnp_id, mispolar):
    """
    Find 'neutral' SNPs, by pairing the non-neutral SNP based on chromosome, mutation and DAF.

    If the SNP is flagged as mispolar then the allele polarization and DAF are inverted.
    """
    dbc = Database(species)

    bins = dbc.get_records('sample_bins')

    bid = []
    lsq = []

    # we need to ensure that there is a comparable number of calls in each bin, so we sort the neutrals by the least
    # squared error of the differences in bin counts, then randomise
    for bin_id in bins:
        bid.append('SUM(s.bin_id = {id}) AS bin{id}'.format(id=bin_id))
        lsq.append('POW(ABS(SUM(s.bin_id = {id}) - nn.bin{id}), 2)'.format(id=bin_id))

    bin_sql = ','.join(bid)
    sqr_sql = '+'.join(lsq)

    modsnps = dbc.get_records_sql("""
        SELECT nn.population,
               nn.id AS modsnp_id,
               {mispolar} AS mispolar,
               ms.id AS neutral_id
          FROM (

        SELECT ms.id,
               msd.population, 
               ms.chrom, 
               IF({mispolar}, ms.ancestral, ms.derived) AS derived,
               IF({mispolar}, ms.derived, ms.ancestral) AS ancestral,
               IF({mispolar}, 1-msd.daf, msd.daf) AS daf,
               {bin_sql}
          FROM modern_snps ms
          JOIN modern_snp_daf msd
            ON msd.modsnp_id = ms.id 
          JOIN sample_reads sr
            ON sr.chrom = ms.chrom
           AND sr.site = ms.site
          JOIN samples s
            ON s.id = sr.sample_id  
         WHERE msd.population = '{population}'
           AND ms.id = {modsnp}
      GROUP BY ms.id

          ) AS nn
          JOIN modern_snps ms
            ON ms.chrom = nn.chrom
           AND ms.derived = nn.derived
           AND ms.ancestral = nn.ancestral
           AND ms.neutral = 1
           AND ms.mispolar IS NULL
           AND ms.variant_id IS NOT NULL
           AND ms.id != nn.id
          JOIN modern_snp_daf msd
            ON msd.modsnp_id = ms.id
           AND msd.population = nn.population
           AND round(msd.daf, 2) = round(nn.daf, 2)
          JOIN sample_reads sr
            ON sr.chrom = ms.chrom
           AND sr.site = ms.site
          JOIN samples s
            ON s.id = sr.sample_id  
      GROUP BY ms.id
      ORDER BY round(({sqr_sql})/{bins})*{bins}, RAND({modsnp})
         LIMIT {num}
           """.format(sqr_sql=sqr_sql, bins=len(bins), mispolar=int(mispolar), bin_sql=bin_sql, population=population,
                      modsnp=modsnp_id, num=NEUTRAL_REPLICATES), key=None)

    # add them to the table so we don't have to query a second time
    for modsnp in modsnps:
        dbc.save_record('selection_neutrals', modsnp)

    return modsnps


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
    :type no_modern: bool
    :type mispolar: bool
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    modsnp = luigi.IntParameter()
    no_modern = luigi.BoolParameter()
    mispolar = luigi.BoolParameter()

    # do not retry after failure, as this just chews CPU cycles
    retry_count = 0

    def requires(self):
        yield DadiBestModel(self.species, self.population, DADI_FOLDED)
        yield AncientSNPsPipeline(self.species)

    def output(self):
        return luigi.LocalTarget('data/selection/{}.input'.format(self.basename))

    def run(self):
        # unpack the inputs
        (_, nref_file, _), _ = self.input()

        # get the Nref population size
        with nref_file.open() as fin:
            pop_size = int(fin.read())

        # time is measured in diffusion units
        diff_units = 2 * pop_size * GENERATION_TIME[self.species]

        # NOTE some SNPs may be mispolarised, so we switch the derived/ancestral alleles in those cases
        derived = 'derived' if not self.mispolar else 'ancestral'

        params = {
            'derived': derived,
            'modsnp': self.modsnp,
            'population': self.population,
            'units': diff_units,
        }

        sql = """
            # get the ancient frequencies in each bin
            SELECT SUM(sr.base = ms.{derived}) AS derived_count,
                   COUNT(sr.id) AS sample_size,
                   -(sb.max / {units}) AS max,
                   -(sb.min / {units}) AS min                   
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
               """.format(**params)

        # noinspection SqlResolve
        modern_sql = """
            # get the modern frequencies
            SELECT {derived}_count AS derived_count,
                   ancestral_count + derived_count AS sample_size,
                   0 AS max,
                   0 AS min
              FROM modern_snps ms
              JOIN modern_snp_daf msd
                ON msd.modsnp_id = ms.id
             WHERE ms.id = {modsnp}
               AND msd.population = '{population}'
               """.format(**params)

        if not self.no_modern:
            sql += " UNION " + modern_sql

        bins = self.dbc.get_records_sql(sql + " ORDER BY max", key=None)

        if len(bins) < MCMC_MIN_BINS:
            raise RuntimeError('ERROR: Insufficient time bins to run `selection` (n={})'.format(len(bins)))

        # write the sample input file
        with self.output().open('w') as tsv_file:
            fields = ['derived_count', 'sample_size', 'max', 'min']
            writer = csv.DictWriter(tsv_file, fieldnames=fields, delimiter='\t')

            # write the data to disk
            for b in bins:
                writer.writerow(b)


class SelectionRunMCMC(utils.PipelineExternalTask):  # TODO put back to `utils.PipelineTask` when done with SLURM
    """
    Run `selection` for the given SNP.

    :type species: str
    :type population: str
    :type modsnp: int
    :type no_modern: bool
    :type mispolar: bool
    :type n: int
    :type s: int
    :type h: float
    :type chain: int
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    modsnp = luigi.IntParameter()
    no_modern = luigi.BoolParameter()
    mispolar = luigi.BoolParameter()
    n = luigi.IntParameter()
    s = luigi.IntParameter()
    h = luigi.FloatParameter()
    chain = luigi.IntParameter()

    # do not retry after failure, as this just chews CPU cycles
    retry_count = 0

    def requires(self):
        yield DadiBestModel(self.species, self.population, DADI_FOLDED)
        yield SelectionInputFile(self.species, self.population, self.modsnp, self.no_modern, self.mispolar)

    def output(self):
        return [luigi.LocalTarget('data/selection/{}.{}'.format(self.basename, ext))
                for ext in ['param.gz', 'time.gz', 'traj.gz', 'log']]

    # TODO uncomment when done with SLURM
    # def run(self):
    #     # compose the input and output file paths
    #     (pop_file, _, _), input_file = self.input()
    #     param_file, time_file, traj_file, log_file = self.output()
    #
    #     output_prefix = utils.trim_ext(log_file.path)
    #
    #     # make a deterministic random seed (helps keep everything easily reproducible)
    #     seed = int('{}{}'.format(self.modsnp, self.chain))
    #
    #     with log_file.open('w') as fout:
    #         # run `selection`
    #         cmd = ['sr',
    #                '-D', input_file.path,   # path to data input file
    #                '-P', pop_file.path,     # path to population size history file
    #                '-o', output_prefix,     # output file prefix
    #                '-a',                    # flag to infer allele age
    #                '-A', MIN_DAF,           # ascertainment in modern individuals
    #                '-n', self.n,            # number of MCMC cycles to run
    #                '-s', self.s,            # frequency of sampling from the posterior
    #                '-h', self.h,            # genetic model (additive, recessive, dominant)
    #                '-f', MCMC_PRINT,        # frequency of printing output to the screen
    #                '-e', seed]              # random number seed
    #
    #         utils.run_cmd(cmd, stdout=fout)


class SelectionBenchmark(utils.MySQLTask):
    """
    Benchmark the speed of the MCMC.

    Models with long trajectories take a lot longer to run.

    :type species: str
    :type population: str
    :type modsnp: int
    :type no_modern: bool
    :type mispolar: bool
    :type n: int
    :type s: int
    :type h: float
    :type chain: int
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    modsnp = luigi.IntParameter()
    no_modern = luigi.BoolParameter()
    mispolar = luigi.BoolParameter()
    n = luigi.IntParameter()
    s = luigi.IntParameter()
    h = luigi.FloatParameter()
    chain = luigi.IntParameter()

    # do not retry after failure, as this just chews CPU cycles
    retry_count = 0

    def requires(self):
        yield DadiBestModel(self.species, self.population, DADI_FOLDED)
        yield SelectionInputFile(self.species, self.population, self.modsnp, self.no_modern, self.mispolar)

    def queries(self):
        # compose the input and output file paths
        (pop_file, _, _), input_file = self.input()

        # generate a unique prefix for temporary files
        output_prefix = 'data/selection/luigi-tmp-{:010}'.format(random.randrange(0, 1e10))

        # make a deterministic random seed (helps keep everything easily reproducible)
        seed = int('{}{}'.format(self.modsnp, self.chain))

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
               '-e', seed]              # random number seed

        start = time()

        # run selection
        utils.run_cmd(cmd)

        # duration is simply the elapsed time
        duration = time() - start

        # now project this forward to get the estimated run time in hours
        est_hours = duration * MCMC_CYCLES / self.n / 3600

        # log the duration
        benchmark = {
            'population': self.population,
            'modsnp_id':  self.modsnp,
            'length':     self.n,
            'thin':       self.s,
            'model':      self.h,
            'no_modern':  self.no_modern,
            'mispolar':   self.mispolar,
            'chain':      self.chain,
            'duration':   duration,
            'est_hours':  est_hours
        }

        self.dbc.save_record('selection_benchmarks', benchmark)

        # tidy up all the temporary files
        for tmp in glob.glob('{}*'.format(output_prefix)):
            os.remove(tmp)


class SelectionPlot(utils.PipelineTask):
    """
    Plot the allele trajectory.

    :type species: str
    :type population: str
    :type modsnp: int
    :type no_modern: bool
    :type mispolar: bool
    :type n: int
    :type s: int
    :type h: float
    :type chain: int
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    modsnp = luigi.IntParameter()
    no_modern = luigi.BoolParameter()
    mispolar = luigi.BoolParameter()
    n = luigi.IntParameter()
    s = luigi.IntParameter()
    h = luigi.FloatParameter()
    chain = luigi.IntParameter()

    resources = {'cpu-cores': 1, 'ram-gb': 24}

    def requires(self):
        yield DadiBestModel(self.species, self.population, DADI_FOLDED)
        yield SelectionInputFile(self.species, self.population, self.modsnp, self.no_modern, self.mispolar)
        yield SelectionRunMCMC(**self.all_params())

    def output(self):
        return luigi.LocalTarget('data/pdf/selection/{}-traj.pdf'.format(self.basename))

    def run(self):
        # unpack the inputs
        (_, nref_file, _), input_file, (param_file, _, _, _) = self.input()
        pdf_file = self.output()

        # compose the input and output file paths
        output_prefix = utils.trim_ext(param_file.path, 2)

        # get the Nref population size
        with nref_file.open() as fin:
            pop_size = int(fin.read())

        # time is measured in diffusion units
        diff_units = 2 * pop_size * GENERATION_TIME[self.species]

        try:
            # plot the allele trajectory
            utils.run_cmd(['Rscript',
                           'rscript/mcmc_plot_selection.R',
                           input_file.path,
                           output_prefix,
                           diff_units,
                           MCMC_BURN_PCT,
                           pdf_file.path])

        except RuntimeError as e:
            # delete the broken PDF
            if os.path.isfile(pdf_file.path):
                pdf_file.remove()

            raise RuntimeError(e)


class SelectionDiagnostics(utils.PipelineTask):
    """
    Diagnostics for an MCMC run.

    Calculate effective sample size (ESS), and plot ESS vs. burn in, autocorrelation, and traces.

    :type species: str
    :type population: str
    :type modsnp: int
    :type no_modern: bool
    :type mispolar: bool
    :type n: int
    :type s: int
    :type h: float
    :type chain: int
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    modsnp = luigi.IntParameter()
    no_modern = luigi.BoolParameter()
    mispolar = luigi.BoolParameter()
    n = luigi.IntParameter()
    s = luigi.IntParameter()
    h = luigi.FloatParameter()
    chain = luigi.IntParameter()

    resources = {'cpu-cores': 1, 'ram-gb': 4}

    # do not retry after failure, as this just chews CPU cycles
    retry_count = 0

    def requires(self):
        yield DadiBestModel(self.species, self.population, DADI_FOLDED)
        yield SelectionRunMCMC(**self.all_params())

    def output(self):
        yield luigi.LocalTarget('data/selection/{}.ess'.format(self.basename))
        yield luigi.LocalTarget('data/selection/{}.map'.format(self.basename))
        yield luigi.LocalTarget('data/selection/{}.diag'.format(self.basename))
        yield luigi.LocalTarget('data/pdf/selection/{}-ess-burn.pdf'.format(self.basename))
        yield luigi.LocalTarget('data/pdf/selection/{}-autocorr.pdf'.format(self.basename))
        yield luigi.LocalTarget('data/pdf/selection/{}-trace-pt1.png'.format(self.basename))

    def run(self):
        # unpack the params
        (_, nref_file, _), (param_file, _, _, _) = self.input()
        ess_file, map_file, diag_file, ess_pdf, autocorr_pdf, trace_png = self.output()

        # get the Nref population size
        with nref_file.open() as fin:
            pop_size = int(fin.read())

        # time is measured in diffusion units
        diff_units = 2 * pop_size * GENERATION_TIME[self.species]

        with diag_file.temporary_path() as diag_path:
            utils.run_cmd(['Rscript',
                           'rscript/mcmc_diagnostics.R',
                           param_file.path,
                           MCMC_BURN_PCT,
                           self.s,
                           diff_units,
                           pop_size,
                           ess_file.path,
                           map_file.path,
                           ess_pdf.path,
                           autocorr_pdf.path,
                           trace_png.path], stdout=open(diag_path, 'w'))


class LoadSelectionDiagnostics(utils.MySQLTask):
    """
    Load the diagnostics for an MCMC run into the db.

    :type species: str
    :type population: str
    :type modsnp: int
    :type no_modern: bool
    :type mispolar: bool
    :type n: int
    :type s: int
    :type h: float
    :type chain: int
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    modsnp = luigi.IntParameter()
    no_modern = luigi.BoolParameter()
    mispolar = luigi.BoolParameter()
    n = luigi.IntParameter()
    s = luigi.IntParameter()
    h = luigi.FloatParameter()
    chain = luigi.IntParameter()

    def requires(self):
        return SelectionDiagnostics(**self.all_params())

    def queries(self):
        # unpack the params
        ess_file, map_file, _, _, _, _ = self.input()

        selection = {
            'population': self.population,
            'modsnp_id':  self.modsnp,
            'length':     self.n,
            'thin':       self.s,
            'model':      self.h,
            'no_modern':  self.no_modern,
            'mispolar':   self.mispolar,
        }

        # get the parent record
        record = self.dbc.get_record('selection', selection)

        if record:
            selection_id = record['id']
        else:
            selection_id = self.dbc.save_record('selection', selection)

        # load the ESS
        with ess_file.open('r') as fin:
            ess = json.load(fin)

        ess['selection_id'] = selection_id
        ess['chain'] = self.chain

        self.dbc.save_record('selection_ess', ess)

        # load the MAP
        with map_file.open('r') as fin:
            posteriori = json.load(fin)

        posteriori['selection_id'] = selection_id
        posteriori['chain'] = self.chain

        self.dbc.save_record('selection_map', posteriori)


class SelectionPSRF(utils.PipelineTask):
    """
    Calculate the Potential Scale Reduction Factor (PSRF) for all the replicate chains.

    AKA. the Gelman and Rubin's convergence diagnostic.

    https://www.rdocumentation.org/packages/coda/versions/0.19-2/topics/gelman.diag

    :type species: str
    :type population: str
    :type modsnp: int
    :type no_modern: bool
    :type mispolar: bool
    :type n: int
    :type s: int
    :type h: float
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    modsnp = luigi.IntParameter()
    no_modern = luigi.BoolParameter()
    mispolar = luigi.BoolParameter()
    n = luigi.IntParameter()
    s = luigi.IntParameter()
    h = luigi.FloatParameter()

    resources = {'cpu-cores': 1, 'ram-gb': 4}

    # do not retry after failure, as this just chews CPU cycles
    retry_count = 0

    def requires(self):
        params = self.all_params()
        for chain in range(1, MCMC_NUM_CHAINS + 1):
            params['chain'] = chain

            yield SelectionRunMCMC(**params)
            yield SelectionPlot(**params)
            yield LoadSelectionDiagnostics(**params)

    def output(self):
        yield luigi.LocalTarget('data/selection/{}-chainAll.ess'.format(self.basename))
        yield luigi.LocalTarget('data/selection/{}-chainAll.psrf'.format(self.basename))
        yield luigi.LocalTarget('data/selection/{}-chainAll.diag'.format(self.basename))
        yield luigi.LocalTarget('data/pdf/selection/{}-chainAll-trace-pt1.png'.format(self.basename))
        yield luigi.LocalTarget('data/pdf/selection/{}-chainAll-gelman-pt1.png'.format(self.basename))

    def run(self):
        param_paths = [param_file.path for param_file in self.input_targets(ext='param.gz')]
        ess_file, psrf_file, diag_file, trace_png, gelman_png = self.output()

        # TODO calculate the MAP for the combined chains

        with diag_file.temporary_path() as diag_path:
            utils.run_cmd(['Rscript',
                           'rscript/mcmc_gelman.R',
                           MCMC_BURN_PCT,
                           self.s,
                           ess_file.path,
                           psrf_file.path,
                           trace_png.path,
                           gelman_png.path] + param_paths, stdout=open(diag_path, 'w'))


class LoadSelectionPSRF(utils.MySQLTask):
    """
    Load the PSRF diagnostics for for all the replicate chains.

    :type species: str
    :type population: str
    :type modsnp: int
    :type no_modern: bool
    :type mispolar: bool
    :type n: int
    :type s: int
    :type h: float
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    modsnp = luigi.IntParameter()
    no_modern = luigi.BoolParameter(default=False)
    mispolar = luigi.BoolParameter(default=False)
    n = luigi.IntParameter(default=MCMC_CYCLES)
    s = luigi.IntParameter(default=MCMC_THIN)
    h = luigi.FloatParameter(default=MODEL_ADDITIVE)

    def requires(self):
        return SelectionPSRF(**self.all_params())

    def queries(self):
        # unpack the params
        ess_file, psrf_file, _, _, _ = self.input()

        selection = {
            'population': self.population,
            'modsnp_id':  self.modsnp,
            'length':     self.n,
            'thin':       self.s,
            'model':      self.h,
            'no_modern':  self.no_modern,
            'mispolar':   self.mispolar,
        }

        selection_id = self.dbc.get_record('selection', selection).pop('id')

        # load the ESS
        with ess_file.open('r') as fin:
            ess = json.load(fin)

        ess['selection_id'] = selection_id
        self.dbc.save_record('selection_ess', ess)

        # load the PSRF
        with psrf_file.open('r') as fin:
            psrf = json.load(fin)

        psrf['selection_id'] = selection_id
        self.dbc.save_record('selection_psrf', psrf)


class SelectionPairNeutrals(utils.PipelineWrapperTask):
    """
    Pair the given modSNP with some neutral replicates, and run `selection` on them all.

    :type species: str
    :type population: str
    :type modsnp: int
    :type no_modern: bool
    :type mispolar: bool
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    modsnp = luigi.IntParameter()
    no_modern = luigi.BoolParameter(default=False)
    mispolar = luigi.BoolParameter(default=False)

    def requires(self):

        yield LoadSelectionPSRF(self.species, self.population, self.modsnp, self.no_modern)

        # get the neutral controls for this modsnp
        neutral = selection_fetch_neutral_snps(self.species, self.population, self.modsnp)

        for modsnp_id in neutral:
            yield LoadSelectionPSRF(self.species, self.population, modsnp_id, self.no_modern)

        if self.mispolar:
            # also run the modsnp with the ancestral/derived alleles reversed
            yield LoadSelectionPSRF(self.species, self.population, self.modsnp, self.no_modern, mispolar=True)

            # and paired to neutral controls using the onverted DAF
            neutral = selection_fetch_neutral_snps(self.species, self.population, self.modsnp, mispolar=True)

            for modsnp_id in neutral:
                yield LoadSelectionPSRF(self.species, self.population, modsnp_id, self.no_modern)


class SelectionGWASPeakSNPs(utils.PipelineWrapperTask):
    """
    Run `selection` on all the direct GWAS hits.

    :type species: str
    :type population: str
    :type no_modern: bool
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    no_modern = luigi.BoolParameter(default=False)

    def requires(self):

        # get the modsnp id for every GWAS hit
        modsnps = self.dbc.get_records_sql("""
            SELECT DISTINCT ms.id, ms.mispolar
              FROM qtls q
              JOIN modern_snps ms
                ON ms.chrom = q.chrom
               AND ms.site = q.site
             WHERE q.associationType = 'Association'
               AND q.valid = 1
               """, key=None)

        for modsnp in modsnps:
            yield SelectionPairNeutrals(self.species, self.population, modsnp['id'], self.no_modern, modsnp['mispolar'])


class SelectionPipeline(utils.PipelineWrapperTask):
    """
    Run all the `selection` jobs.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        for pop in self.list_populations(ancient=True):
            for no_modern in [True, False]:
                yield SelectionGWASPeakSNPs(self.species, pop, no_modern)


class SelectionBenchmarkGWASNeutrals(utils.PipelineWrapperTask):
    """
    Benchmark all the GWAS hits and their neutrals.

    :type species: str
    :type population: str
    :type no_modern: bool
    :type step: int
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    no_modern = luigi.BoolParameter(default=False)
    step = luigi.IntParameter()

    def requires(self):

        offset = BENCHMARK_LIMIT * self.step

        modsnps = self.dbc.get_records_sql("""
            SELECT DISTINCT modsnp_id AS id, IFNULL(mispolar, 0) mispolar
              FROM selection_neutrals
             WHERE population = '{population}'
               AND mispolar = 0
             
             UNION 
            
            SELECT DISTINCT neutral_id  AS id, IFNULL(mispolar, 0) mispolar
              FROM selection_neutrals
             WHERE population = '{population}'
               AND mispolar = 0
               
             LIMIT {limit} OFFSET {offset} 
               """.format(population=self.population, limit=BENCHMARK_LIMIT, offset=offset), key=None)

        params = self.all_params()
        params.pop('step')

        for modsnp in modsnps:
            params['modsnp'] = modsnp['id']
            params['mispolar'] = bool(modsnp['mispolar'])
            params['n'] = MCMC_BENCHMARK
            params['s'] = MCMC_THIN
            params['h'] = MODEL_ADDITIVE

            for chain in range(1, MCMC_NUM_CHAINS + 1):
                params['chain'] = chain

                yield SelectionBenchmark(**params)


class SelectionTidyPipeline(utils.PipelineWrapperTask):
    """
    Tidy up all the `selection` jobs that were completed by SLURM.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):

        # get all the completed models
        completed = glob.glob('data/selection/{}-*.param.gz'.format(self.species))

        for filename in completed:
            # e.g. data/selection/horse-DOM2-modsnp9876899-n100000000-s100-h0.5-chain1.param.gz
            _, population, modsnp, n, s, h, _ = os.path.basename(filename).split('-')

            yield LoadSelectionPSRF(self.species, population, int(modsnp[6:]), n=int(n[1:]), s=int(s[1:]), h=float(h[1:]))


if __name__ == '__main__':
    luigi.run()
