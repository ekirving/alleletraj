#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.plink import PlinkTask, PlinkHighGeno

# the maximum number of migration events in Treemix
TREEMIX_MAX_M = 5

# number of SNPs to group for LD
TREEMIX_K = 1000

# what level should Treemix group by, pops OR samples
GROUP_BY_POPS = 'grp-pops'
GROUP_BY_SMPL = 'grp-smpl'


class PlinkBedToFreq(PlinkTask):
    """
    Convert a BED file into a minor allele frequency report, needed for input into Treemix.

    :type species: str
    :type groupby: str
    """
    species = luigi.Parameter()
    groupby = luigi.Parameter()

    def requires(self):
        return PlinkHighGeno(self.species)

    def output(self):
        return [luigi.LocalTarget('data/plink/{}.{}'.format(self.basename, ext)) for ext in ['frq.strat.gz', 'log']]

    def run(self):
        # unpack the inputs/outputs
        bed_path, bim_path, fam_path, _ = [in_file.path for in_file in self.input()]
        _, log_file = self.output()

        # to group by samples, we need to reassign them to their own families
        if self.groupby == GROUP_BY_SMPL:
            # replace family with sample code
            fam = utils.run_cmd(["awk '{$1=$2}$0' " + fam_path.path], shell=True)

            # make a new fam file
            fam_path = utils.insert_suffix(fam_path, GROUP_BY_SMPL)
            with open(fam_path, 'w') as fout:
                fout.write(fam)

        utils.run_cmd(['plink',
                       '--chr-set', self.chrset,
                       '--freq', 'gz',  # make a gzipped MAF report
                       '--family',      # group by population
                       '--bed', bed_path,
                       '--bim', bim_path,
                       '--fam', fam_path,
                       '--out', utils.trim_ext(log_file.path)])


class TreemixPlinkFreq(utils.PipelineTask):
    """
    Convert a Plink MAF frequency file into Treemix format

    :type species: str
    :type groupby: str
    """
    species = luigi.Parameter()
    groupby = luigi.Parameter()

    def requires(self):
        return PlinkBedToFreq(self.species, self.groupby)

    def output(self):
        return luigi.LocalTarget('treemix/{}.frq.gz'.format(self.basename))

    def run(self):
        # unpack the inputs/outputs
        frq_input, _ = self.input()
        frq_output = self.output()

        with frq_output.temporary_path() as frq_path:
            # convert the file
            utils.run_cmd(['python',
                           'plink2treemix.py',
                           frq_input.path,
                           frq_path])


class TreemixM(utils.PipelineTask):
    """
    Run Treemix with the given number of migration events.

    :type species: str
    :type groupby: str
    :type m: int
    """
    species = luigi.Parameter()
    groupby = luigi.Parameter()
    m = luigi.IntParameter(default=0)

    def requires(self):
        return TreemixPlinkFreq(self.species, self.groupby)

    def output(self):
        return [luigi.LocalTarget('treemix/{}.{}'.format(self.basename, ext)) for ext in ['cov.gz', 'log']]

    def run(self):
        # unpack the inputs/outputs
        frq_file = self.input()
        covgz_file, log_file, = self.output()

        # get the stub name of the output files
        out_stub = utils.trim_ext(covgz_file.path, 2)

        # run treemix
        cmd = ['treemix',
               '-i',    frq_file.path,
               '-root', self.outgroup,
               '-k',    TREEMIX_K,   # group together "k" SNPs to account for linkage disequilibrium
               '-m',    self.m,      # build the ML graph with "m" migration events
               '-o',    out_stub]

        log = utils.run_cmd(cmd)

        # save the log file
        with log_file.open('w') as fout:
            fout.write(log)


class TreemixPlotM(utils.PipelineTask):
    """
    Plot the Treemix output for the given number of migration events.

    :type species: str
    :type groupby: str
    :type m: int
    """
    species = luigi.Parameter()
    groupby = luigi.Parameter()
    m = luigi.IntParameter(default=0)

    def requires(self):
        return TreemixM(self.species, self.groupby, self.m)

    def output(self):
        yield luigi.LocalTarget('treemix/{}.tips'.format(self.basename))
        yield luigi.LocalTarget('pdf/treemix/{}.treemix.pdf'.format(self.basename))

    def run(self):
        # unpack the inputs/outputs
        covgz_file, log_file, = self.input()
        tips_file, pdf_file = self.output()

        # extract the tip names
        tips = utils.run_cmd(['gunzip -c {} | head -n1'.format(covgz_file.path)], shell=True).split()

        # save them to a file, for use by the r script
        with tips_file.open('w') as fout:
            fout.write('\n'.join(tips))

        # get the stub name of the output files
        out_stub = utils.trim_ext(covgz_file.path, 2)

        # TODO fix errors here
        with pdf_file.temporary_path() as pdf_path:
            # plot the treemix tree
            utils.run_cmd(['Rscript',
                           'rscript/treemix-plot.R',
                           out_stub,
                           self.groupby,
                           tips_file.path,
                           pdf_path])


class PipelineTreemix(utils.PipelineWrapperTask):
    """
    Run all the Treemix plots

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        for m in range(TREEMIX_MAX_M + 1):
            yield TreemixPlotM(self.species, GROUP_BY_SMPL, m)


if __name__ == '__main__':
    luigi.run()
