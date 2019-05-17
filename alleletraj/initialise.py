#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.db.load import CreateDatabase
from alleletraj.ensembl.load import EnsemblLoadPipeline
from alleletraj.qtl.load import PopulateAllLoci
from alleletraj.samples import LoadAllSamples
from alleletraj.snpchip.load import SNPChipLoadPipeline


class Setup(utils.PipelineWrapperTask):
    """
    Setup the DB and load the samples.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        # create a new db and add all the empty tables
        yield CreateDatabase(self.species)

        # populate the ensembl genes and variants tables
        yield EnsemblLoadPipeline(self.species)

        # populate the snpchip tables
        yield SNPChipLoadPipeline(self.species)

        # populate all the QTLs from the AnimalQTL db, and other regions of interest
        yield PopulateAllLoci(self.species)

        # load all samples
        yield LoadAllSamples(self.species)


if __name__ == '__main__':
    luigi.run()
