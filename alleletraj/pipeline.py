#!/usr/bin/env python
# -*- coding: utf-8 -*-

import luigi

from alleletraj import utils
from alleletraj.database.create import CreateDatabase
from alleletraj.modern.load_snps import ModernSNPsPipeline
from alleletraj.modern.demography import DadiPipeline
from alleletraj.ensembl import EnsemblPipeline
from alleletraj.snpchip import SNPChipPipeline
from alleletraj.qtl.qtls import QTLPipeline
from alleletraj.ancient.samples import SamplesPipeline
from alleletraj.ancient.load_snps import AncientSNPsPipeline
from alleletraj.qtl.analyse import AnalyseQTLsPipeline


class RunAll(utils.PipelineWrapperTask):
    """
    Build the species database

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):

        # create a new database and add all the empty tables
        yield CreateDatabase(self.species)

        # ascertain SNPs in modern whole genome data
        yield ModernSNPsPipeline(self.species)

        # find the best fitting ∂a∂i model
        yield DadiPipeline(self.species)

        # link all the Ensembl data to the modern SNPs
        yield EnsemblPipeline(self.species)

        # link the SNPChip data to the modern SNPs
        yield SNPChipPipeline(self.species)

        # load the QTLs from the AnimalQTL database, and other regions of interest
        yield QTLPipeline(self.species)

        # load the sample metadata
        yield SamplesPipeline(self.species)

        # load the sample reads for each ascertained SNP
        yield AncientSNPsPipeline(self.species)

        # analyse the coverage and quality for SNPs in each QTLs
        yield AnalyseQTLsPipeline(self.species)

        # if self.species == 'pig':
        #     # pick the best SNPs to target for a capture array
        #     yield AscertainmentPipeline(self.species)
        #
        # # run `selection` on all the 'best' QTL SNPs
        # yield SelectionBestQTLSNPs(self.species)

        # TODO make a plotting pipeline
        # yield GraphsPipeline(self.species)


if __name__ == '__main__':
    luigi.run()
