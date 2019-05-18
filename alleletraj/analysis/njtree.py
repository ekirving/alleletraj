#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.plink import PlinkDistMatrix

# types of NJ tree
NJTREE_PHYLO = 'phylogram'
NJTREE_FAN = 'fan'


class NeighborJoiningTree(utils.DatabaseTask):
    """
    Create a neighbor joining phylogenetic tree from a distance matrix.

    :type species: str
    :type tree: str
    """
    species = luigi.Parameter()
    tree = luigi.Parameter(default=NJTREE_PHYLO)

    def requires(self):
        return PlinkDistMatrix(self.species)

    def output(self):
        return [luigi.LocalTarget('data/njtree/{}.{}'.format(self.basename, ext)) for ext in ['tree', 'pdf']]

    def run(self):
        # unpack the inputs/outputs
        data_file, _ = self.input()
        tree_file, pdf_file = self.output()

        with pdf_file.temporary_path() as pdf_path:
            # generate a tree from the labeled data
            utils.run_cmd(['Rscript',
                           'rscript/plot_njtree.R',
                           data_file.path,
                           self.tree,
                           self.outgroup,
                           tree_file.path,
                           pdf_path])


if __name__ == '__main__':
    luigi.run()
