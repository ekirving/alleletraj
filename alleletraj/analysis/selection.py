#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
import json
import os

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
MCMC_REPLICATES = 4

# number of MCMC cycles to run
MCMC_CYCLES = int(5e7)

# fraction of MCMC cycles to discard as burn in
MCMC_BURN_PCT = 0.2

# thinning in analytically unnecessary, but makes the MCMC run much faster
# see https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.2041-210X.2011.00131.x
MCMC_THIN = 1000

# frequency of printing output to the screen
MCMC_PRINT = 1000

# fraction of the allele frequency to update during a trajectory update move
MCMC_FRACTION = 20  # TODO what units is this in?

# genetic model
MODEL_RECESSIVE = 0
MODEL_ADDITIVE = 0.5
MODEL_DOMINANT = 1

# minimum number of time bins
MCMC_MIN_BINS = 3

# number of DAF paired neutral SNPs to run for every non-neutral SNP
NEUTRAL_REPLICATES = 5


def selection_neutral_snps(species, population, modsnp_id):
    """
    Find 'neutral' SNPs, by pairing the non-neutral SNP based on chromosome, mutation type and DAF.

    If the SNP is flagged as mispolar then the DAF is inverted.
    """
    dbc = Database(species)

    # TODO should we also consider the SNP coverage?
    modsnps = dbc.get_records_sql("""
        SELECT ms.id
          FROM (

        SELECT msd.population, ms.chrom, type, IF(ms.mispolar, 1-msd.daf, msd.daf) AS daf
          FROM modern_snps ms
          JOIN modern_snp_daf msd
            ON msd.modsnp_id = ms.id 
         WHERE msd.population = '{population}'
           AND ms.id = {modsnp}

          ) AS nn
          JOIN modern_snps ms
            ON ms.chrom = nn.chrom
           AND ms.type = nn.type
           AND ms.neutral = 1
           AND ms.mispolar IS NULL
           AND ms.variant_id IS NOT NULL
          JOIN modern_snp_daf msd
            ON msd.modsnp_id = ms.id
           AND msd.population = nn.population
           AND round(msd.daf, 2) = round(nn.daf, 2)
      ORDER BY RAND({modsnp})
         LIMIT {num}
           """.format(population=population, modsnp=modsnp_id, num=NEUTRAL_REPLICATES)).keys()

    if len(modsnps) != NEUTRAL_REPLICATES:
        raise RuntimeError('ERROR: Insufficient neutral SNPs to run `selection` (n={})'.format(len(modsnps)))

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
            SELECT IFNULL(SUM(base = derived), 0) AS derived_count,
                   COUNT(read_id) AS sample_size,
                   -(sb.max / {units}) AS max,
                   -(sb.min / {units}) AS min                    
              FROM sample_bins sb
         LEFT JOIN (
                    SELECT s.bin_id,
                           sr.id AS read_id,
                           sr.base,
                           ms.{derived} AS derived
                      FROM modern_snps ms
                      JOIN sample_reads sr
                        ON sr.chrom = ms.chrom
                       AND sr.site = ms.site
                      JOIN samples s
                        ON s.id = sr.sample_id
                     WHERE ms.id = {modsnp}
                       AND s.population = '{population}'
                   
                   ) AS rec ON rec.bin_id = sb.id
          GROUP BY sb.id
               """.format(**params)

        # noinspection SqlResolve
        modern_sql = """
            # get the modern frequencies
            SELECT {derived}_count, ancestral_count + derived_count, 0, 0
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


class SelectionRunMCMC(utils.PipelineTask):
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

    def requires(self):
        yield DadiBestModel(self.species, self.population, DADI_FOLDED)
        yield SelectionInputFile(self.species, self.population, self.modsnp, self.no_modern, self.mispolar)

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

        param_path = utils.trim_ext(param_file.path)
        time_path = utils.trim_ext(time_file.path)
        traj_path = utils.trim_ext(traj_file.path)

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
            # delete the unfinished .time and .traj files (as these can be very large)
            os.remove(time_path)
            os.remove(traj_path)

            # but keep the param file, as this may be useful for diagnosing the error
            utils.run_cmd(['gzip', param_path])

            raise RuntimeError(e)

        else:
            # gzip the output files
            utils.run_cmd(['gzip', param_path])
            utils.run_cmd(['gzip', time_path])
            utils.run_cmd(['gzip', traj_path])


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

    resources = {'cpu-cores': 1, 'ram-gb': 64}  # TODO need to refactor Josh's code to fix massive memory requirement

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

    def requires(self):
        yield DadiBestModel(self.species, self.population, DADI_FOLDED)
        yield SelectionRunMCMC(**self.all_params())

    def output(self):
        yield luigi.LocalTarget('data/selection/{}.ess'.format(self.basename))
        yield luigi.LocalTarget('data/selection/{}.map'.format(self.basename))
        yield luigi.LocalTarget('data/selection/{}.diag'.format(self.basename))
        yield luigi.LocalTarget('data/pdf/selection/{}-ess-burn.pdf'.format(self.basename))
        yield luigi.LocalTarget('data/pdf/selection/{}-autocorr.pdf'.format(self.basename))
        yield luigi.LocalTarget('data/pdf/selection/{}-trace.pdf'.format(self.basename))

    def run(self):
        # unpack the params
        (_, nref_file, _), (param_file, _, _, _) = self.input()
        ess_file, map_file, diag_file, ess_pdf, autocorr_pdf, trace_pdf = self.output()

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
                           trace_pdf.path], stdout=open(diag_path, 'w'))


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
            map = json.load(fin)

        map['selection_id'] = selection_id
        map['chain'] = self.chain

        self.dbc.save_record('selection_map', map)


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

    def requires(self):
        params = self.all_params()
        for chain in range(1, MCMC_REPLICATES + 1):
            params['chain'] = chain

            yield SelectionRunMCMC(**params)
            yield SelectionPlot(**params)
            yield LoadSelectionDiagnostics(**params)

    def output(self):
        yield luigi.LocalTarget('data/selection/{}-chainAll.ess'.format(self.basename))
        yield luigi.LocalTarget('data/selection/{}-chainAll.psrf'.format(self.basename))
        yield luigi.LocalTarget('data/selection/{}-chainAll.diag'.format(self.basename))
        yield luigi.LocalTarget('data/pdf/selection/{}-chainAll-trace.pdf'.format(self.basename))
        yield luigi.LocalTarget('data/pdf/selection/{}-chainAll-gelman.pdf'.format(self.basename))

    def run(self):
        param_paths = [param_file.path for param_file in self.input_targets(ext='param.gz')]
        ess_file, psrf_file, diag_file, trace_pdf, gelman_pdf = self.output()

        with diag_file.temporary_path() as diag_path:
            utils.run_cmd(['Rscript',
                           'rscript/mcmc_gelman.R',
                           MCMC_BURN_PCT,
                           self.s,
                           ess_file.path,
                           psrf_file.path,
                           trace_pdf.path,
                           gelman_pdf.path] + param_paths, stdout=open(diag_path, 'w'))


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

        selection_id = self.dbc.get_record('selection', selection)['id']

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


class SelectionGWASSNPs(utils.PipelineWrapperTask):
    """
    Run `selection` on all the direct GWAS hits.

    :type species: str
    :type population: str

    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    no_modern = luigi.BoolParameter(default=False)

    def requires(self):

        # get the modsnp id for every GWAS hit
        modsnps = self.dbc.get_records_sql("""
            SELECT DISTINCT ms.id, ms.mispolar
              FROM qtls q
              JOIN qtl_snps qs
                ON qs.qtl_id = q.id
              JOIN modern_snps ms
                ON qs.modsnp_id = ms.id 
              JOIN ensembl_variants ev              
                ON ms.variant_id = ev.id
               AND ev.rsnumber = q.peak
             WHERE q.associationType = 'Association'
               AND q.valid = 1""", key=None)

        for modsnp in modsnps:

            # get the neutral controls for this modsnp
            neutral = selection_neutral_snps(self.species, self.pop, modsnp['id'])

            yield LoadSelectionPSRF(self.species, self.population, modsnp['id'], self.no_modern)

            if modsnp['mispolar']:
                # run the modsnp with the ancestral/derived alleles reversed
                yield LoadSelectionPSRF(self.species, self.population, modsnp['id'], self.no_modern, mispolar=True)

            for modsnp_id in neutral:
                yield LoadSelectionPSRF(self.species, self.population, modsnp_id, self.no_modern)


class SelectionPipeline(utils.PipelineWrapperTask):
    """
    Run all the `selection` jobs.


    :type species: str
    :type population: str
    """
    species = luigi.Parameter()
    population = luigi.OptionalParameter(default=None)

    def requires(self):

        # use all ancient populations, if population is not set
        populations = [self.population] if self.population else self.list_populations(ancient=True)

        # TODO run all SNPs without the modern frequency to see how different the inference is

        for pop in populations:
            yield SelectionGWASSNPs(self.species, pop)


if __name__ == '__main__':
    luigi.run()
