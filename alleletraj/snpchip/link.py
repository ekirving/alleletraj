#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.modern.snps import LoadModernSNPs
from alleletraj.snpchip.load import SNPChipLoadPipeline


class LinkSNPChipVariants(utils.MySQLTask):
    """
    Link modern SNPs to their SNPchip variants

    :type species: str
    :type chrom: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    db_lock_tables = ['modern_snps_{chrom}']

    def requires(self):
        yield SNPChipLoadPipeline(self.species)
        yield LoadModernSNPs(self.species, self.chrom)

    def queries(self):
        self.dbc.execute_sql("""
            UPDATE modern_snps ms
              JOIN ensembl_variants ev
                ON ev.id = ms.variant_id
              JOIN snpchip sc
                ON sc.rsnumber = ev.rsnumber 
               SET ms.snpchip_id = sc.id
             WHERE ms.chrom = '{chrom}'
               """.format(chrom=self.chrom))


class SNPChipLinkPipeline(utils.PipelineWrapperTask):
    """
    Populate the snpchip_* tables and link modern_snps records to their snpchip variants.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        for chrom in self.chromosomes:
            yield LinkSNPChipVariants(self.species, chrom)


if __name__ == '__main__':
    luigi.run()
