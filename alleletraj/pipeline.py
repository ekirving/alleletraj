#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.ancient.graphs import GraphsPipeline
from alleletraj.ancient.samples import SamplesPipeline
from alleletraj.ancient.selection import SelectionBestQTLSNPs
from alleletraj.ancient.snps import AncientSNPsPipeline
from alleletraj.ascertain import AscertainmentPipeline
from alleletraj.database.load import CreateDatabase
from alleletraj.ensembl.link import EnsemblLinkPipeline
from alleletraj.modern.demog import DadiPipeline
from alleletraj.modern.snps import ModernSNPsPipeline
from alleletraj.qtl.analyse import AnalyseQTLsPipeline
from alleletraj.qtl.qtls import QTLPipeline
from alleletraj.snpchip import SNPChipPipeline


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

        # link all the Ensembl data to the modern SNPs
        yield EnsemblLinkPipeline(self.species)

        # link the SNPChip data to the modern SNPs
        yield SNPChipPipeline(self.species)

        # find the best fitting ∂a∂i model
        yield DadiPipeline(self.species)

        # load the QTLs from the AnimalQTL database, and other regions of interest
        yield QTLPipeline(self.species)

        # load the sample metadata
        yield SamplesPipeline(self.species)

        # load the sample reads for each ascertained SNP
        yield AncientSNPsPipeline(self.species)

        # analyse the coverage and quality for SNPs in each QTLs
        yield AnalyseQTLsPipeline(self.species)

        if self.species == 'pig':
            # pick the best SNPs to target for a capture array
            yield AscertainmentPipeline(self.species)

        # run `selection` on all the 'best' QTL SNPs
        yield SelectionBestQTLSNPs(self.species)

        # make all the plots
        yield GraphsPipeline(self.species)


if __name__ == '__main__':
    luigi.run()
