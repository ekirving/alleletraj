#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi
import unicodecsv as csv

# local modules
from alleletraj import utils
from alleletraj.ancient.snps import AncientSNPsPipeline


class GraphDerivedVersusAge(utils.PipelineTask):
    """
    Produce a scatter plot of the oldest observed derived allele vs. the oldest sample at that site.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        yield AncientSNPsPipeline(self.species)

    def output(self):
        return [luigi.LocalTarget('data/pdf/{}-snps-ages.{}'.format(self.species, ext)) for ext in ['pdf', 'png']]

    # noinspection SqlResolve, SqlAmbiguousColumn
    def run(self):
        pdf_file, png_file = self.output()

        dbc = self.db_conn()

        # sql fragment to calculate the median age of each sample
        median = "COALESCE(sd.median, (COALESCE(c14.lower, sd.lower)+COALESCE(c14.upper, sd.upper))/2)"

        # get the age of every covered snp
        reads = dbc.get_records_sql("""
            SELECT sr.id AS read_id,
                   {median} AS median
               FROM (
                        # get a unique list of GWAS peaks
                        SELECT min(id) AS id
                          FROM qtls
                         WHERE associationType = 'Association'
                           AND valid = 1
                      GROUP BY peak
                          
                    ) AS q
               JOIN qtl_snps qs
                 ON q.id = qs.qtl_id
               JOIN modern_snps ms
                 ON ms.id = qs.modsnp_id
               JOIN sample_reads sr
                 ON sr.chrom = ms.chrom
                AND sr.site = ms.site
               JOIN samples s
                 ON s.id = sr.sample_id
          LEFT JOIN sample_dates sd
                 ON s.age = sd.age
          LEFT JOIN sample_dates_c14 c14
                 ON c14.accession = s.accession
            HAVING median IS NOT NULL""".format(median=median), key=None)

        with open("tsv/{}-snps-counts.tsv".format(self.species), "wb") as tsv_file:

            fields = ['read_id', 'median']
            writer = csv.DictWriter(tsv_file, fieldnames=fields, delimiter='\t')
            writer.writeheader()

            # write the data to disk
            for read in reads:
                writer.writerow(read)

        # get the dates of the oldest sample and oldest observation of the derived allele for every QTL SNP
        snps = dbc.get_records_sql("""
             SELECT qs.modsnp_id,
                    t.class AS trait,
                    ms.ancestral,
                    ms.derived,
                    MAX({median}) AS oldest_sample,
                    MAX(IF(sr.base = ms.derived, {median}, NULL)) AS oldest_derived
               FROM (
                        # get a unique list of GWAS peaks
                        SELECT min(id) AS id, trait_id
                          FROM qtls
                         WHERE associationType = 'Association'
                           AND valid = 1
                      GROUP BY peak
                          
                    ) AS q
               JOIN qtl_snps qs
                 ON q.id = qs.qtl_id
               JOIN traits t
                 ON t.id = q.trait_id
               JOIN modern_snps ms
                 ON ms.id = qs.modsnp_id
               JOIN sample_reads sr
                 ON sr.chrom = ms.chrom
                AND sr.site = ms.site
               JOIN samples s
                 ON s.id = sr.sample_id
          LEFT JOIN sample_dates sd
                 ON s.age = sd.age
          LEFT JOIN sample_dates_c14 c14
                 ON c14.accession = s.accession
              WHERE qs.num_reads IS NOT NULL
           GROUP BY qs.modsnp_id""".format(median=median), key=None)

        with open("tsv/{}-snps-ages.tsv".format(self.species), "wb") as tsv_file:

            fields = ['modsnp_id', 'trait', 'ancestral', 'derived', 'oldest_sample', 'oldest_derived']
            writer = csv.DictWriter(tsv_file, fieldnames=fields, delimiter='\t')
            writer.writeheader()

            # write the data to disk
            for snp in snps:
                writer.writerow(snp)

        # now generate the plot
        utils.run_cmd(['Rscript', 'rscript/plot-age-derived.R', utils.trim_ext(pdf_file.path),
                       'Derived Allele vs. Sample Age'])

        # convert to PNG
        utils.run_cmd(["convert -density 300 {pdf} -quality 100 {png}".format(pdf=pdf_file.path, png=png_file.path)],
                      shell=True)


class GraphsPipeline(utils.PipelineWrapperTask):
    """
    Print all the graphs.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        yield GraphDerivedVersusAge(self.species)


if __name__ == '__main__':
    luigi.run()
