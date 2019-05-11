#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
import os
import re
from collections import defaultdict
from shutil import copyfile

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.const import CPU_CORES_MED
from alleletraj.plink import PlinkPruneBed

# the number of bootstrap replicates to run
ADMIXTURE_BOOTSTRAP = 0  # TODO put this back to 100

# the sort order options for plotting admixture column charts
ADMIXTURE_SORT_K = 1
ADMIXTURE_SORT_ORDER = 0

# the largest value of K to model
ADMIXTURE_MAX_K = 20


class AdmixtureK(utils.PipelineTask):
    """
    Run admixture, with K ancestral populations, on the pruned reference data BED file.

    :type species: str
    :type k: int
    """
    species = luigi.Parameter()
    k = luigi.IntParameter()

    resources = {'cpu-cores': CPU_CORES_MED}

    def requires(self):
        return PlinkPruneBed(self.species)

    def output(self):
        # admixture does not let us control the output name
        basename = utils.trim_path_ext(self.input()[0].path)
        return [luigi.LocalTarget('data/admix/{}.{}.{}'.format(basename, self.k, ext)) for ext in ['P', 'Q', 'log']]

    def run(self):
        # unpack the inputs/outputs
        bed_input, _, _, _ = self.input()
        _, _, log_file = self.output()

        # CPU threads to use
        threads = self.resources['cpu-cores']

        # admixture only outputs to the current directory, so we need to change the working directory
        os.chdir('data/admixture')

        bed_file = '../{}'.format(bed_input.path)

        cmd = ['admixture',
               '-j{}'.format(threads),              # use multi-threading
               '-B{}'.format(ADMIXTURE_BOOTSTRAP),  # the number of bootstrap replicates to run
               '--cv=10',                           # generate cross-validation estimates
               bed_file,                            # using the pruned data file
               self.k]                              # for K ancestral populations

        log = utils.run_cmd(cmd)

        # restore previous working directory
        os.chdir('../..')

        # save the log file
        with log_file.open('w') as fout:
            fout.write(log)


class AdmixtureSortK(utils.PipelineTask):
    """
    Sort the admixture Q matrix for the given value of K, such that the column orders maintain optimal relative
    relationships.

    :type species: str
    :type k: int
    """
    species = luigi.Parameter()
    k = luigi.IntParameter()

    def requires(self):
        yield AdmixtureK(self.species, self.k)

        # sorting of k requires that k-1 also be sorted (because we need something to compare against)
        yield AdmixtureSortK(self.species, self.k - 1) if self.k > 1 else []

    def output(self):
        return luigi.LocalTarget('data/admix/{}.sorted.Q'.format(self.basename))

    def run(self):
        # unpack the inputs/outputs
        (_, unsorted_q, _), sorted_in = self.input()
        sorted_out = self.output()

        # no need to sort a file with only one k
        if self.k == 1:
            copyfile(unsorted_q.path, sorted_out.path)
            return

        with sorted_out.temporary_path() as out_path:
            # sort all the matrices
            utils.run_cmd(['Rscript',
                           'rscript/admix-sort-k.R',
                           unsorted_q.path,  # input Q file (unsorted)
                           sorted_in.path,   # sorted k-1 Q file
                           out_path])        # output Q file (sorted)


class AdmixturePlotK(utils.PipelineTask):
    """
    Use ggplot to plot the admixture Q stats

    :type species: str
    :type maxk: int
    :type k: int
    """
    species = luigi.Parameter()
    maxk = luigi.IntParameter(significant=False)
    k = luigi.IntParameter()

    def requires(self):
        yield AdmixtureSortK(self.species, self.k)
        yield PlinkPruneBed(self.species)

    def output(self):
        yield luigi.LocalTarget('data/admix/{}.data'.format(self.basename))
        yield luigi.LocalTarget('data/pdf/admix/{}-admix.pdf'.format(self.basename))

    def run(self):
        # unpack the inputs/outputs
        sorted_q, (_, _, fam_file, _) = self.input()
        data_file, pdf_file = self.output()

        # use awk and paste to add population and sample names, needed for the plot
        awk = "awk '{ print $1 \" \" $2 }' " + fam_file.path + " | paste - " + sorted_q.path

        # get the admix data
        data = utils.run_cmd([awk], shell=True)

        # parse the data into a dict, so it can be output in a specific order
        datadict = defaultdict(list)
        for line in data.strip().split('\n'):
            datadict[line.split()[0]].append(line)

        # compose the header row
        header = ['Pop{}'.format(i) for i in range(1, self.k + 1)]
        header.insert(0, 'Sample')
        header.insert(0, 'Population')

        # save the labeled file
        with data_file.open('w') as fout:
            # output the header row
            fout.write('\t'.join(header) + '\n')

            # TODO replace with sorting in R
            # output the populations, in the chosen order
            for pop in DATASETS[self.species]:
                for line in datadict[pop]:
                    fout.write(line + '\n')

        with pdf_file.temporary_path() as pdf_path:
            # generate a PDF of the admixture stacked column chart
            utils.run_cmd(['Rscript',
                           'rscript/admix-plot-k.R',
                           data_file.path,
                           pdf_path,
                           self.maxk,
                           ADMIXTURE_SORT_ORDER])


class AdmixtureCV(utils.PipelineTask):
    """
    Run admixture for the given population, determine the optimal K value, and plot the graphs

    :type species: str
    :type maxk: int
    """
    species = luigi.Parameter()
    maxk = luigi.IntParameter(default=None, significant=False)

    def requires(self):
        if not self.maxk:
            # set maximum K to be the number of actual populations
            self.maxk = min(len(DATASETS[self.species]), ADMIXTURE_MAX_K)

        # run admixture for each value of K
        for k in range(1, self.maxk + 1):
            yield AdmixturePlotK(self.species, self.maxk, k)

    def output(self):
        yield luigi.LocalTarget('data/admix/{}.cv.data'.format(self.basename))
        yield luigi.LocalTarget('data/pdf/admix/{}-admix.cv.pdf'.format(self.basename))

    def run(self):
        # unpack the outputs
        data_file, pdf_file = self.output()

        # use grep to extract the cross-validation scores from all the log files
        cvs = utils.run_cmd(['grep -h CV admix/{}.*.log'.format(self.basename)], shell=True)

        # extract the K value and CV score
        # e.g. "CV error (K=1): 1.20340"
        data = [tuple([re.sub(r'[^\d.]+', '', val) for val in cv.split(':')]) for cv in cvs.splitlines()]

        # write the scores to a data file
        with data_file.open('w') as fout:
            fout.write('\t'.join(['K', 'CV']) + '\n')
            for row in data:
                fout.write('\t'.join(str(datum) for datum in row) + '\n')

        with pdf_file.temporary_path() as pdf_path:
            # plot the CV values as a line graph
            utils.run_cmd(['Rscript',
                           'rscript/plot-line-graph.R',
                           data_file.path,
                           pdf_path,
                           'Ancestral populations (K)',
                           'Cross-validation Error'])


if __name__ == '__main__':
    luigi.run()
