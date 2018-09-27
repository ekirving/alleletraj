#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

from populate_coverage import *
from discover_modern_snps import discover_modern_snps, links_modern_snps
from discover_snps import discover_snps
from analyse_qtls import analyse_qtls
from populate_samples import populate_pig_samples, populate_horse_samples
from ascertainment import perform_ascertainment
from graph_derived import graph_derived
from pipeline_utils import *

from luigi.contrib.mysqldb import MySqlTarget

class LoadEnsemblGenes(PipelineTask):
    """
    Load the GTF (General Transfer Format) data from Ensembl.

    GTF provides access to all annotated transcripts which make up an Ensembl gene set.

    See https://www.ensembl.org/info/website/upload/gff.html#fields
    """
    def output(self):
        return luigi.LocalTarget("selection/{}.input".format(self.basename))

    def run(self):
        load_ensembl_genes()


class LoadEnsemblVariants(PipelineTask):
    """
    Load the GVF (Genome Variation Format) data from Ensembl.

    Contains all germline variants from the current Ensembl release for this species.

    See http://www.sequenceontology.org/gvf.html
    """
    def output(self):
        return luigi.LocalTarget("selection/{}.input".format(self.basename))

    def run(self):
        load_ensembl_variants()

        self.output().touch()


class LoadSNPChipVariants(PipelineTask):
    """
    Load the SNP chip variants from SNPchimp

    See http://bioinformatics.tecnoparco.org/SNPchimp
    """
    def output(self):
        return luigi.LocalTarget("selection/{}.input".format(self.basename))

    def run(self):
        load_snpchip_variants()


class SomeTask(PipelineTask):
    """

    """
    species = luigi.Parameter()

    def requires(self):
        return SomeTask(self.species)

    def output(self):
        return luigi.LocalTarget("selection/{}.input".format(self.basename))

    def run(self):

        pass

class SomeTask(PipelineTask):
    """

    """
    species = luigi.Parameter()

    def requires(self):
        return SomeTask(self.species)

    def output(self):
        return luigi.LocalTarget("selection/{}.input".format(self.basename))

    def run(self):

        pass

class SomeTask(PipelineTask):
    """

    """
    species = luigi.Parameter()

    def requires(self):
        return SomeTask(self.species)

    def output(self):
        return luigi.LocalTarget("selection/{}.input".format(self.basename))

    def run(self):

        pass

class SomeTask(PipelineTask):
    """

    """
    species = luigi.Parameter()

    def requires(self):
        return SomeTask(self.species)

    def output(self):
        return luigi.LocalTarget("selection/{}.input".format(self.basename))

    def run(self):

        pass


class SelectionHorseTest(luigi.WrapperTask):
    """
    Run `selection` on all a sub-set of SNPs to test MCMC params.
    """

    def requires(self):

        pass


if __name__ == '__main__':
    luigi.run()
