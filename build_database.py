#!/usr/bin/env python
# -*- coding: utf-8 -*-

import luigi

from pipeline_utils import PipelineWrapperTask

from pipeline_database import CreateDatabase
from pipeline_modern_snps import ModernSNPsPipeline
from pipeline_ensembl import EnsemblPipeline
from pipeline_snpchip import SNPChipPipeline
from pipeline_qtls import QTLPipeline
from pipeline_samples import SamplesPipeline
from pipeline_sample_reads import SampleReadsPipeline
from pipeline_discover_snps import DiscoverSNPsPipeline
from pipeline_analyse_qtls import AnalyseQTLsPipeline
from pipeline_ascertainment import AscertainmentPipeline
from pipeline_selection import SelectionBestQTLSNPs


class BuildDatabase(PipelineWrapperTask):
    """
    Build the species database

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):

        # create a new database and add all the empty tables
        yield CreateDatabase(self.species)

        # TODO rethink the use of 'population' in modern SNPs (perhaps replace with pop column in table?)
        # ascertain SNPs in modern whole genome data
        yield ModernSNPsPipeline(self.species)

        # link all the Ensembl data to the modern SNPs
        yield EnsemblPipeline(self.species)

        # link the SNPChip data to the modern SNPs
        yield SNPChipPipeline(self.species)

        # load the QTLs from the AnimalQTL database, and other regions of interest
        yield QTLPipeline(self.species)

        # load the sample metadata
        yield SamplesPipeline(self.species)

        # load the sample reads for each ascertained SNP
        yield SampleReadsPipeline(self.species)

        # apply quality filters to the sample reads
        yield DiscoverSNPsPipeline(self.species)

        # analyse the coverage and quality for SNPs in each QTLs
        yield AnalyseQTLsPipeline(self.species)

        # if self.species == 'pig':
        #     # pick the best SNPs to target for a capture array
        #     yield AscertainmentPipeline(self.species)
        #
        # # run `selection` on all the 'best' QTL SNPs
        # yield SelectionBestQTLSNPs(self.species)
        #
        # # TODO make a plotting pipeline


if __name__ == '__main__':
    luigi.run()

