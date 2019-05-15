#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.samples import LoadSamplesPipeline
from alleletraj.db.load import CreateDatabase


class Setup(utils.PipelineWrapperTask):
    """
    Setup the DB and load the samples.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        # create a new db and add all the empty tables
        yield CreateDatabase(self.species)

        # TODO load the modern samples

        # load all the ancient samples
        yield LoadSamplesPipeline(self.species)


if __name__ == '__main__':
    luigi.run()
