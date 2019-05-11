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


class NeighborJoiningTree(utils.PipelineTask):
    """
    Create a neighbor joining phylogenetic tree from a distance matrix

    :type dataset: str
    :type treetype: str
    """
    dataset = luigi.Parameter()
    treetype = luigi.Parameter()

    def requires(self):
        return PlinkDistMatrix(self.dataset)

    def output(self):
        yield luigi.LocalTarget('data/njtree/{}.tree'.format(self.basename))
        yield luigi.LocalTarget('data/pdf/njtree/{}.njtree.pdf'.format(self.basename))

    def run(self):
        # unpack the inputs/outputs
        data_file, _ = self.input()
        tree_file, pdf_file = self.output()

        with pdf_file.temporary_path() as pdf_path:
            # generate a tree from the labeled data
            utils.run_cmd(['Rscript',
                           'rscript/plot-phylo-tree.R',
                           data_file.path,
                           self.treetype,
                           self.outgroup,
                           tree_file.path,
                           pdf_path])


if __name__ == '__main__':
    luigi.run()
