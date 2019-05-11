#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.const import CPU_CORES_MED
from alleletraj.eigenstrat import ConvertfBedToEigenstrat

# which pairs of PCA components should we print
PCA_COMPONENTS = [(1, 2), (3, 4), (5, 6)]


class SmartPCA(utils.PipelineTask):
    """
    Calculate the eigenvectors for the given dataset

    :type species: str
    :type project: list
    """
    species = luigi.Parameter()
    project = luigi.ListParameter(default=None)

    resources = {'cpu-cores': CPU_CORES_MED}

    def requires(self):
        return ConvertfBedToEigenstrat(self.species)

    def output(self):
        return [luigi.LocalTarget('data/smartpca/{}.{}'.format(self.basename, ext)) for ext in
                ['pca.evec', 'eval', 'poplist', 'par']]

    def run(self):
        # unpack the inputs/outputs
        geno_file, snp_file, ind_file, _, _ = self.input()
        evec_file, eval_file, pop_list, par_file = self.output()

        # TODO fix this
        # default to all ancient populations
        project = ANCIENT_POPS if self.project is None else self.project

        # tell smartpca which pops to use for calculating the eigenvectors, and by inference, which to project
        with pop_list.open('w') as fout:
            fout.write('\n'.join([pop for pop in self.populations if pop not in project]))

        # compose the config settings for smartpca
        config = [
            'genotypename:   {}'.format(geno_file.path),
            'snpname:        {}'.format(snp_file.path),
            'indivname:      {}'.format(ind_file.path),
            'evecoutname:    {}'.format(evec_file.path),
            'evaloutname:    {}'.format(eval_file.path),
            'numthreads:     {}'.format(CPU_CORES_MED),  # number of threads to use
            'poplistname:    {}'.format(pop_list.path),  # the list of pops to calculate the eigenvectors from
            'lsqproject:     YES',                       # use least squares projection, best for missing data
            'numoutlieriter: 0',                         # don't exclude outliers
            'inbreed:        YES',                       # compute inbreeding stats
        ]

        with par_file.open('w') as fout:
            fout.write('\n'.join(config))

        # calculate the PCA and project the ancient samples on top
        utils.run_cmd(['smartpca', '-p', par_file.path])


class SmartPCAPlot(utils.PipelineTask):
    """
    Use ggplot to plot the PCA

    :type species: str
    :type project: list
    :type components: list
    """
    species = luigi.Parameter()
    project = luigi.ListParameter(default=None)
    components = luigi.ListParameter(default=PCA_COMPONENTS, significant=False)

    def requires(self):
        return SmartPCA(self.species, self.project)

    def output(self):
        for ext in ['pve', 'calc', 'proj']:
            yield luigi.LocalTarget('smartpca/{}.{}'.format(self.basename, ext))

        for pc1, pc2 in self.components:
            yield luigi.LocalTarget('pdf/smartpca/{}-pca({},{}).pdf'.format(self.basename, pc1, pc2))

    def run(self):
        # unpack the inputs/outputs
        evec_file, eval_file, pop_list, par_file = self.input()
        (pve_file, calc_file, proj_file), pdf_files = list(self.output())[:3], list(self.output())[3:]

        # calculate the percentage of variance explained by each PC, by dividing each eigenvalue by the sum total
        var_sum = utils.run_cmd(["awk '{ sum+=$1 } END {print sum}' " + eval_file.path], shell=True)
        var_pve = utils.run_cmd(["awk '{print $1/" + var_sum.strip() + "}' " + eval_file.path], shell=True)

        # save the percentage variance
        with pve_file.open('w') as fout:
            fout.write(var_pve)

        # copy the population name into the first column, skip the header row (to make things easier for the Rscript)
        awk = "awk 'NR>1 {print $NF $0}' " + evec_file.path

        # default to all ancient populations
        project = ANCIENT_POPS if self.project is None else self.project

        # make a regex to match only the ancient pops
        pop_regex = '|'.join(['^{} '.format(pop) for pop in project])

        # split the data into calculated and projected
        calc = utils.run_cmd([awk + ' | grep -Ev "{}"'.format(pop_regex)], shell=True)
        proj = utils.run_cmd([awk + ' | grep -E  "{}"'.format(pop_regex)], shell=True)

        # save the pca data
        with calc_file.open('w') as fout:
            fout.write(calc)

        with proj_file.open('w') as fout:
            fout.write(proj)

        # plot the first 6 components
        for (pc1, pc2), pdf_target in zip(self.components, pdf_files):

            # generate both labeled and unlabeled PDFs
            pdfs = [(0, pdf_target), (1, luigi.LocalTarget(utils.insert_suffix(pdf_target.path, 'labeled')))]

            for labeled, pdf_file in pdfs:

                with pdf_file.temporary_path() as pdf_path:
                    # generate a PDF of the PCA plot
                    utils.run_cmd(['Rscript',
                                   'rscript/smartpca-plot.R',
                                   calc_file.path,  # pca data, used for calculating the eigenvectors
                                   proj_file.path,  # projected pca data
                                   pve_file.path,   # pve data (% variance)
                                   pdf_path,        # location to save the pdf file
                                   pc1,             # component num for x-axis
                                   pc2,             # component num for y-axis
                                   labeled])        # show point labels (0/1)


if __name__ == '__main__':
    luigi.run()
