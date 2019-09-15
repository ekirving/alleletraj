#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
import glob
import itertools
import json
import os
import random

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
MCMC_NUM_CHAINS = 2

# maximum number of MCMC replicates to run in search of converged runs
MCMC_MAX_CHAINS = 6

# number of MCMC cycles to run
MCMC_CYCLES = int(1e7)

# fraction of MCMC cycles to discard as burn in
MCMC_BURN_PCT = 0.5

# thinning is analytically unnecessary, but makes the MCMC run much faster (as Josh's method writes directly to disk)
# see https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.2041-210X.2011.00131.x
MCMC_THIN = 100

# fraction of the path to update (i.e. length/F)
MCMC_FRACTION = 20

# frequency of printing output to the log
MCMC_PRINT = 10000

# the minimum ESS threshold for an MCMC run
MCMC_MIN_ESS = 100

# the maximum MPSRF threshold for a set of MCMC runs
MCMC_MAX_MPSRF = 1.2

# genetic model
MODEL_RECESSIVE = 0
MODEL_ADDITIVE = 0.5
MODEL_DOMINANT = 1

# minimum number of time bins
MCMC_MIN_BINS = 3

# number of DAF paired neutral SNPs to run for every non-neutral SNP
NEUTRAL_REPLICATES = 5


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

    modsnps = dbc.get_records('selection_neutrals', params, key='neutral_id')

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
    :type const_pop: bool
    :type no_age: bool
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    modsnp = luigi.IntParameter()
    no_modern = luigi.BoolParameter()
    mispolar = luigi.BoolParameter()
    const_pop = luigi.BoolParameter()
    no_age = luigi.BoolParameter()

    # do not retry after failure, as this just chews CPU cycles
    retry_count = 0

    def requires(self):
        yield DadiBestModel(self.species, self.population, DADI_FOLDED, self.const_pop)
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
                   -round((sb.max / {units}), 4) AS max,
                   -round(((sb.max + sb.min) / 2 / {units}), 4) AS median,
                   -round((sb.min / {units}), 4) AS min                   
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
                   0.0 AS max,
                   0.0 AS median,
                   0.0 AS min
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
            if self.no_age:
                fields = ['derived_count', 'sample_size', 'median', 'median']
            else:
                fields = ['derived_count', 'sample_size', 'max', 'min']

            writer = csv.DictWriter(tsv_file, fieldnames=fields, delimiter='\t', extrasaction='ignore')

            # write the data to disk
            for b in bins:
                writer.writerow(b)


class SelectionRunMCMC(utils.PipelineTask):
    """
    Run `selection` for the given SNP.

    :type species: str
    :type population: str
    :type modsnp: int
    :type no_modern: bool
    :type mispolar: bool
    :type const_pop: bool
    :type no_age: bool
    :type n: int
    :type s: int
    :type h: float
    :type F: int
    :type chain: int
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    modsnp = luigi.IntParameter()
    no_modern = luigi.BoolParameter()
    mispolar = luigi.BoolParameter()
    const_pop = luigi.BoolParameter()
    no_age = luigi.BoolParameter()
    n = luigi.IntParameter()
    s = luigi.IntParameter()
    h = luigi.FloatParameter()
    F = luigi.IntParameter()
    chain = luigi.IntParameter()

    # TODO remove when done with SLURM jobs on the cluster
    # resources = {'SLURM': 2}

    # do not retry after failure, as this just chews CPU cycles
    retry_count = 0

    def requires(self):
        yield DadiBestModel(self.species, self.population, DADI_FOLDED, self.const_pop)
        yield SelectionInputFile(self.species, self.population, self.modsnp, self.no_modern, self.mispolar,
                                 self.const_pop, self.no_age)

    def output(self):
        return [luigi.LocalTarget('data/selection/{}.{}'.format(self.basename, ext))
                for ext in ['param.gz', 'time.gz', 'traj.gz', 'log']]

    def run(self):
        # compose the input and output file paths
        (pop_file, _, _), input_file = self.input()
        param_file, time_file, traj_file, log_file = self.output()

        output_prefix = utils.trim_ext(log_file.path)

        # make a deterministic random seed (helps keep everything easily reproducible)
        seed = int('{}{}'.format(self.modsnp, self.chain))

        with log_file.open('w') as fout:
            # run `selection`
            cmd = ['sr',
                   '-D', input_file.path,   # path to data input file
                   '-P', pop_file.path,     # path to population size history file
                   '-o', output_prefix,     # output file prefix
                   '-A', MIN_DAF,           # ascertainment in modern individuals
                   '-n', self.n,            # number of MCMC cycles to run
                   '-s', self.s,            # frequency of sampling from the posterior
                   '-h', self.h,            # genetic model (additive, recessive, dominant)
                   '-F', self.F,            # fraction of the path to update (i.e. length/F)
                   '-f', MCMC_PRINT,        # frequency of printing output to the screen
                   '-e', seed]              # random number seed

            if not self.no_age:
                cmd += ['-a']               # flag to infer allele age

            utils.run_cmd(cmd, stdout=fout)


class SelectionPlot(utils.PipelineTask):
    """
    Plot the allele trajectory.

    :type species: str
    :type population: str
    :type modsnp: int
    :type no_modern: bool
    :type mispolar: bool
    :type const_pop: bool
    :type no_age: bool
    :type n: int
    :type s: int
    :type h: float
    :type F: int
    :type chain: int
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    modsnp = luigi.IntParameter()
    no_modern = luigi.BoolParameter()
    mispolar = luigi.BoolParameter()
    const_pop = luigi.BoolParameter()
    no_age = luigi.BoolParameter()
    n = luigi.IntParameter()
    s = luigi.IntParameter()
    h = luigi.FloatParameter()
    F = luigi.IntParameter()
    chain = luigi.IntParameter()

    resources = {'cpu-cores': 1, 'ram-gb': 4}

    # do not retry after failure, as this just chews CPU cycles
    retry_count = 0

    def requires(self):
        yield DadiBestModel(self.species, self.population, DADI_FOLDED, self.const_pop)
        yield SelectionInputFile(self.species, self.population, self.modsnp, self.no_modern, self.mispolar,
                                 self.const_pop, self.no_age)
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
    :type const_pop: bool
    :type no_age: bool
    :type n: int
    :type s: int
    :type h: float
    :type F: int
    :type chain: int
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    modsnp = luigi.IntParameter()
    no_modern = luigi.BoolParameter()
    mispolar = luigi.BoolParameter()
    const_pop = luigi.BoolParameter()
    no_age = luigi.BoolParameter()
    n = luigi.IntParameter()
    s = luigi.IntParameter()
    h = luigi.FloatParameter()
    F = luigi.IntParameter()
    chain = luigi.IntParameter()

    resources = {'cpu-cores': 1, 'ram-gb': 4}

    # do not retry after failure, as this just chews CPU cycles
    retry_count = 0

    def requires(self):
        yield DadiBestModel(self.species, self.population, DADI_FOLDED, self.const_pop)
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
    :type const_pop: bool
    :type no_age: bool
    :type n: int
    :type s: int
    :type h: float
    :type F: int
    :type chain: int
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    modsnp = luigi.IntParameter()
    no_modern = luigi.BoolParameter()
    mispolar = luigi.BoolParameter()
    const_pop = luigi.BoolParameter()
    no_age = luigi.BoolParameter()
    n = luigi.IntParameter()
    s = luigi.IntParameter()
    h = luigi.FloatParameter()
    F = luigi.IntParameter()
    chain = luigi.IntParameter()

    def requires(self):
        return SelectionDiagnostics(**self.all_params())

    def queries(self):
        # unpack the params
        ess_file, map_file, _, _, _, _ = self.input()

        selection = {
            'population': self.population,
            'modsnp_id':  self.modsnp,
            'no_modern':  self.no_modern,
            'mispolar':   self.mispolar,
            'const_pop':  self.const_pop,
            'no_age':     self.no_age,
            'length':     self.n,
            'thin':       self.s,
            'model':      self.h,
            'frac':       self.F,
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


class SelectionCalculateMPSRF(utils.PipelineTask):
    """
    Calculate the MPSRF for the given set of chains.

    If MPSRF > 1.2 then SelectionPSRF will call this task again with a different set of chains to compare.

    :type species: str
    :type population: str
    :type modsnp: int
    :type no_modern: bool
    :type mispolar: bool
    :type const_pop: bool
    :type no_age: bool
    :type n: int
    :type s: int
    :type h: float
    :type F: int
    :type chains: list
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    modsnp = luigi.IntParameter()
    no_modern = luigi.BoolParameter()
    mispolar = luigi.BoolParameter()
    const_pop = luigi.BoolParameter()
    no_age = luigi.BoolParameter()
    n = luigi.IntParameter()
    s = luigi.IntParameter()
    h = luigi.FloatParameter()
    F = luigi.IntParameter()
    chains = luigi.ListParameter()

    def output(self):
        return luigi.LocalTarget('data/selection/{}.mpsrf'.format(self.basename))

    def run(self):
        param_paths = []

        params = self.all_params()
        params.pop('chains')

        for chain in self.chains:
            params['chain'] = chain

            # get the path to the MCMC parameter chain
            mcmc = yield SelectionRunMCMC(**params)
            param_paths.append(mcmc[0].path)

        with self.output().temporary_path() as mpsrf_path:
            utils.run_cmd(['Rscript',
                           'rscript/mcmc_mpsrf.R',
                           MCMC_BURN_PCT,
                           self.s] + param_paths, stdout=open(mpsrf_path, 'w'))


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
    :type const_pop: bool
    :type no_age: bool
    :type n: int
    :type s: int
    :type h: float
    :type F: int
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    modsnp = luigi.IntParameter()
    no_modern = luigi.BoolParameter()
    mispolar = luigi.BoolParameter()
    const_pop = luigi.BoolParameter()
    no_age = luigi.BoolParameter()
    n = luigi.IntParameter()
    s = luigi.IntParameter()
    h = luigi.FloatParameter()
    F = luigi.IntParameter()

    resources = {'cpu-cores': 1, 'ram-gb': 4}

    # do not retry after failure, as this just chews CPU cycles
    retry_count = 0

    def requires(self):
        # setup the basic replicate chain requirements, these are extended at runtime depending on quality metrics
        params = self.all_params()
        for chain in range(1, MCMC_NUM_CHAINS + 1):
            params['chain'] = chain
            yield SelectionDiagnostics(**params)

    def output(self):
        yield luigi.LocalTarget('data/selection/{}-chainAll.ess'.format(self.basename))
        yield luigi.LocalTarget('data/selection/{}-chainAll.psrf'.format(self.basename))
        yield luigi.LocalTarget('data/selection/{}-chainAll.diag'.format(self.basename))
        yield luigi.LocalTarget('data/pdf/selection/{}-chainAll-trace-pt1.png'.format(self.basename))
        yield luigi.LocalTarget('data/pdf/selection/{}-chainAll-gelman-pt1.png'.format(self.basename))

    def run(self):
        ess_file, psrf_file, diag_file, trace_png, gelman_png = self.output()

        mpsrf_good = []
        ess_good = []
        chain = 1

        while not mpsrf_good:
            if chain > MCMC_MAX_CHAINS:
                raise RuntimeError('Failed to converge MCMC chains after {} attempts'.format(MCMC_MAX_CHAINS))

            # get the next MCMC chain
            params = self.all_params()
            params['chain'] = chain
            diag_task = yield SelectionDiagnostics(**params)

            # add the metrics to the db
            yield LoadSelectionDiagnostics(**params)

            # and plot the trajectory
            yield SelectionPlot(**params)

            # get the ESS of the chain
            with next(diag_task).open('r') as fin:
                ess = json.load(fin)

            # get the minimum ESS
            min_ess = min(ess.values())  # TODO consider replacing with a multivariate ESS

            # enforce the threshold
            if min_ess > MCMC_MIN_ESS:
                ess_good.append(chain)

                if len(ess_good) >= MCMC_NUM_CHAINS:
                    # get all possible combinations of the MCMC chains (in sets of size MCMC_NUM_CHAINS)
                    chain_sets = list(itertools.combinations(ess_good, MCMC_NUM_CHAINS))

                    params = self.all_params()
                    for chain_set in chain_sets:
                        params['chains'] = chain_set

                        # calculation the MPSRF for this set of chains
                        mpsrf_file = yield SelectionCalculateMPSRF(**params)

                        with mpsrf_file.open('r') as fin:
                            mpsrf = float(fin.read())

                            if mpsrf <= MCMC_MAX_MPSRF:
                                mpsrf_good = chain_set
                                break
            chain += 1

        param_paths = []

        # now we have a suitable set, let's finish things off
        params = self.all_params()
        for chain in mpsrf_good:
            params['chain'] = chain

            # get the path to the MCMC parameter chain
            mcmc = yield SelectionRunMCMC(**params)
            param_paths.append(mcmc[0].path)

        # TODO this should calculate the maximum a posteriori (MAP) based on the combined chains
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
    :type const_pop: bool
    :type no_age: bool
    :type n: int
    :type s: int
    :type h: float
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    modsnp = luigi.IntParameter()
    no_modern = luigi.BoolParameter(default=False)
    mispolar = luigi.BoolParameter(default=False)
    const_pop = luigi.BoolParameter(default=False)
    no_age = luigi.BoolParameter(default=False)
    n = luigi.IntParameter(default=MCMC_CYCLES)
    s = luigi.IntParameter(default=MCMC_THIN)
    h = luigi.FloatParameter(default=MODEL_ADDITIVE)
    F = luigi.IntParameter(default=MCMC_FRACTION)

    def requires(self):
        return SelectionPSRF(**self.all_params())

    def queries(self):
        # unpack the params
        ess_file, psrf_file, _, _, _ = self.input()

        selection = {
            'population': self.population,
            'modsnp_id':  self.modsnp,
            'no_modern':  self.no_modern,
            'mispolar':   self.mispolar,
            'const_pop':  self.const_pop,
            'no_age':     self.no_age,
            'length':     self.n,
            'thin':       self.s,
            'model':      self.h,
            'frac':       self.F,
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

        for param in psrf:
            if psrf[param] in ('Inf', 'NA', 'NaN'):
                psrf[param] = None

        psrf['selection_id'] = selection_id
        self.dbc.save_record('selection_psrf', psrf)


class SelectionPairNeutrals(utils.MySQLTask):  # TODO this task type is misleading
    """
    Pair the given modSNP with some neutral replicates, and run `selection` on them all.

    But, wait until the target SNP has been successfully modelled before launching the neutral replicates.

    :type species: str
    :type population: str
    :type modsnp: int
    :type no_modern: bool
    :type mispolar: bool
    :type const_pop: bool
    :type no_age: bool
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    modsnp = luigi.IntParameter()
    no_modern = luigi.BoolParameter(default=False)
    mispolar = luigi.BoolParameter(default=False)
    const_pop = luigi.BoolParameter(default=False)
    no_age = luigi.BoolParameter(default=False)

    def requires(self):
        return LoadSelectionPSRF(**self.all_params())

    def run(self):
        # get the neutral controls for this modsnp
        neutrals = selection_fetch_neutral_snps(self.species, self.population, self.modsnp, self.mispolar)

        # TODO only run neutrals if the target SNP looks interesting
        params = self.all_params()
        for modsnp_id in neutrals:
            params['modsnp'] = modsnp_id
            yield LoadSelectionPSRF(**self.all_params())


class SelectionGWASPeakSNPs(utils.PipelineWrapperTask):
    """
    Run `selection` on all the direct GWAS hits (except any flagged as mispolar)

    :type species: str
    :type population: str
    :type no_modern: bool
    :type mispolar: bool
    :type const_pop: bool
    :type no_age: bool
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    no_modern = luigi.BoolParameter(default=False)
    mispolar = luigi.BoolParameter(default=False)
    const_pop = luigi.BoolParameter(default=False)
    no_age = luigi.BoolParameter(default=False)

    def requires(self):

        # get the modsnp id for every GWAS hit
        modsnps = self.dbc.get_records_sql("""
            SELECT DISTINCT ms.id
              FROM qtls q
              JOIN modern_snps ms
                ON ms.chrom = q.chrom
               AND ms.site = q.site
             WHERE q.associationType = 'Association'
               AND q.valid = 1
               AND IFNULL(ms.mispolar, 0) = {mispolar}
               """.format(mispolar=int(self.mispolar)))

        params = self.all_params()
        for modsnp in modsnps:
            params['modsnp'] = modsnp
            yield SelectionPairNeutrals(**params)


class SelectionPipeline(utils.PipelineWrapperTask):
    """
    Run all the `selection` jobs.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        for pop in self.list_populations(ancient=True):
            for no_modern in [True, False]:
                for const_pop in [True, False]:
                    yield SelectionGWASPeakSNPs(self.species, pop, no_modern, const_pop)


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
            # TODO this will break with const_pop
            # e.g. data/selection/horse-DOM2-modsnp9876899-n100000000-s100-h0.5-chain1.param.gz
            population, modsnp, n, s, h = os.path.basename(filename).split('-')[1:6]

            yield LoadSelectionPSRF(self.species, population, int(modsnp[6:]), n=int(n[1:]), s=int(s[1:]),
                                    h=float(h[1:]))


if __name__ == '__main__':
    luigi.run()
