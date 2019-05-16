#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.ensembl.load import LoadEnsemblGenes, LoadEnsemblVariants
from alleletraj.modern.snps import ModernSNPsPipeline


class LinkEnsemblGenes(utils.DatabaseTask):
    """
    Link modern SNPs to their Ensembl genes

    :type species: str
    :type chrom: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    db_lock_tables = ['modern_snps_{chrom}']

    def requires(self):
        yield LoadEnsemblGenes(self.species)
        yield ModernSNPsPipeline(self.species)

    def output(self):
        return luigi.LocalTarget('data/db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        exec_time = self.dbc.execute_sql("""
            UPDATE modern_snps ms
              JOIN ensembl_genes eg
                ON eg.chrom = ms.chrom
               AND ms.site BETWEEN eg.start AND eg.end
               SET ms.gene_id = eg.id
             WHERE ms.chrom = '{chrom}'""".format(chrom=self.chrom))

        with self.output().open('w') as fout:
            fout.write('Execution took {}'.format(exec_time))


class LinkEnsemblVariants(utils.DatabaseTask):
    """
    Link modern SNPs to their Ensembl dbsnp variants

    :type species: str
    :type chrom: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    db_lock_tables = ['modern_snps_{chrom}']

    def requires(self):
        yield LoadEnsemblVariants(self.species)
        yield ModernSNPsPipeline(self.species)

    def output(self):
        return luigi.LocalTarget('data/db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        exec_time = self.dbc.execute_sql("""
            UPDATE modern_snps ms
              JOIN ensembl_variants v
                ON ms.chrom = v.chrom
               AND ms.site = v.start
               SET ms.variant_id = v.id
             WHERE ms.chrom = '{chrom}'
               AND v.type = 'SNV'
               AND CHAR_LENGTH(alt) = 1
               AND v.ref IN (ms.derived, ms.ancestral)
               AND v.alt IN (ms.derived, ms.ancestral)""".format(chrom=self.chrom))

        with self.output().open('w') as fout:
            fout.write('Execution took {}'.format(exec_time))


class EnsemblLinkPipeline(utils.PipelineWrapperTask):
    """
    Link modern_snps records to their ensembl genes and dbsnp variants.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        for chrom in self.chromosomes:
            yield LinkEnsemblGenes(self.species, chrom)
            yield LinkEnsemblVariants(self.species, chrom)


if __name__ == '__main__':
    luigi.run()
