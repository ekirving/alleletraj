#!/usr/bin/env python
# -*- coding: utf-8 -*-

import luigi
import os
import gzip

# import my custom modules
# TODO make these into PipelineTask properties
from pipeline_consts import SAMPLES, CHROM_SIZE, MAX_INSERT_SIZE
from pipeline_modern_snps import ProcessSNPs
from pipeline_utils import PipelineTask, db_conn, curl_download

# NOTE these must be paired with the reference genome used for each species
ENSEMBL_URLS = {

    'pig': {
        'gtf': 'ftp://ftp.ensembl.org/pub/release-89/gtf/sus_scrofa/Sus_scrofa.Sscrofa10.2.89.gtf.gz',
        'gvf': 'ftp://ftp.ensembl.org/pub/release-89/variation/gvf/sus_scrofa/Sus_scrofa.gvf.gz'
    },

    'horse': {
        'gtf': 'ftp://ftp.ensembl.org/pub/release-92/gtf/equus_caballus//Equus_caballus.EquCab2.92.gtf.gz',
        'gvf': 'ftp://ftp.ensembl.org/pub/release-92/variation/gvf/equus_caballus/equus_caballus.gvf.gz'
    },

    'goat':  {
        'gtf': 'ftp://ftp.ensembl.org/pub/release-92/gtf/capra_hircus/Capra_hircus.ARS1.92.gtf.gz',
        'gvf': 'ftp://ftp.ensembl.org/pub/release-92/variation/gvf/capra_hircus/capra_hircus.gvf.gz'
    },
}


class DownloadEnsemblData(PipelineTask):
    """
    Fetches the data from the Ensembl FTP site.

    :type species: str
    :type type: str
    """
    species = luigi.Parameter()
    type = luigi.Parameter()

    @property
    def url(self):
        return ENSEMBL_URLS[self.species][self.type]

    def output(self):
        return luigi.LocalTarget('ensembl/{}'.format(os.path.basename(self.url)))

    def run(self):
        with self.output().temporary_path() as tmp_path:
            curl_download(self.url, tmp_path)


class LoadEnsemblGenes(PipelineTask):
    """
    Load the GTF (General Transfer Format) data from Ensembl.

    GTF provides access to all annotated transcripts which make up an Ensembl gene set.

    See https://www.ensembl.org/info/website/upload/gff.html#fields

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return DownloadEnsemblData(self.species, 'gtf')

    def output(self):
        return luigi.LocalTarget('ensembl/{}-genes.log'.format(self.species))

    def run(self):
        # open a db connection
        dbc = db_conn(self.species)

        # the column headers for batch inserting into the db
        fields = ('source', 'gene_id', 'version', 'biotype', 'chrom', 'start', 'end')

        # open the GTF file
        with gzip.open(self.input().path, 'r') as fin:

            # buffer the records for bulk insert
            records = []

            for line in fin:
                # strip flanking whitespace
                line = line.strip()

                # skip comments
                if not line.startswith('#'):

                    # grab the GFF columns
                    chrom, source, feature, start, end, score, strand, frame, attributes = line.split('\t')

                    # we only want genes
                    if feature == 'gene':
                        # unpack the 'key1 "value1";key2 "value2";' attribute pairs
                        atts = dict([pair.strip().split(' ') for pair in attributes.split(';') if pair != ''])

                        # remove quoting
                        atts.update({k: v.strip('"') for k, v in atts.items()})

                        # setup the record to insert, in this order
                        gene = (atts['gene_source'], atts['gene_id'], atts['gene_version'], atts['gene_biotype'], chrom,
                                start, end)

                        # add to the buffer
                        records.append(gene)

            # bulk insert all the records
            dbc.save_records('ensembl_genes', fields, records)

        with self.output().open('w') as fout:
            fout.write('Inserted {} Ensembl gene records'.format(len(records)))


class LoadEnsemblVariants(PipelineTask):
    """
    Load the GVF (Genome Variation Format) data from Ensembl.

    Contains all germline variants from the current Ensembl release for this species.

    See http://www.sequenceontology.org/gvf.html

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return DownloadEnsemblData(self.species, 'gvf')

    def output(self):
        return luigi.LocalTarget('ensembl/{}-variants.log'.format(self.species))

    def run(self):
        dbc = db_conn(self.species)

        # the column headers for batch inserting into the db
        fields = ('dbxref', 'rsnumber', 'type', 'chrom', 'start', 'end', 'ref', 'alt')

        # open the GVF file
        with gzip.open(self.input().path, 'r') as fin:

            num_recs = 0

            # buffer the records for bulk insert
            records = []

            for line in fin:
                # strip flanking whitespace
                line = line.strip()

                # skip comments
                if not line.startswith('#'):
                    num_recs += 1

                    # grab the GFF columns
                    chrom, source, snptype, start, end, score, strand, phase, attributes = line.split('\t')

                    # unpack the 'key1=value1;key2=value2;' attribute pairs
                    atts = dict([pair.strip().split('=') for pair in attributes.split(';') if pair != ''])

                    # unpack dbxref
                    dbxref, rsnumber = atts['Dbxref'].split(':')

                    # setup the record to insert, in this order
                    variant = (dbxref, rsnumber, snptype, chrom, start, end, atts['Reference_seq'], atts['Variant_seq'])

                    # add to the buffer
                    records.append(variant)

                    # bulk insert in chunks
                    if num_recs % MAX_INSERT_SIZE == 0:
                        dbc.save_records('ensembl_variants', fields, records)
                        records = []

            # insert any remaining records
            if records:
                dbc.save_records('ensembl_variants', fields, records)

        with self.output().open('w') as fout:
            fout.write('Inserted {} Ensembl variant records'.format(num_recs))


class LinkEnsemblGenes(PipelineTask):
    """
    Link modern SNPs to their Ensembl genes

    :type species: str
    :type population: str
    :type chrom: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        yield LoadEnsemblGenes(self.species)
        yield ProcessSNPs(self.species, self.population, self.chrom)

    def output(self):
        return luigi.LocalTarget('db/{}-ensembl_genes.log'.format(self.basename))

    def run(self):
        dbc = db_conn(self.species)

        exec_time = dbc.execute_sql("""
            UPDATE modern_snps ms
              JOIN ensembl_genes eg
                ON eg.chrom = ms.chrom
               AND ms.site BETWEEN eg.start AND eg.end
               SET ms.gene_id = eg.id
             WHERE ms.population = '{pop}'
               AND ms.chrom = '{chrom}'""".format(pop=self.population, chrom=self.chrom))

        with self.output().open('w') as fout:
            fout.write('Execution took {}'.format(exec_time))


class LinkEnsemblVariants(PipelineTask):
    """
    Link modern SNPs to their Ensembl dbsnp variants

    :type species: str
    :type population: str
    :type chrom: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        yield LoadEnsemblVariants(self.species)
        yield ProcessSNPs(self.species, self.population, self.chrom)

    def output(self):
        return luigi.LocalTarget('db/{}-ensembl_vars.log'.format(self.basename))

    def run(self):
        dbc = db_conn(self.species)

        exec_time = dbc.execute_sql("""
            UPDATE modern_snps ms
              JOIN ensembl_variants v
                ON ms.chrom = v.chrom
               AND ms.site = v.start
               SET ms.variant_id = v.id
             WHERE ms.population = '{pop}'
               AND ms.chrom = '{chrom}'
               AND v.type = 'SNV'        
               AND CHAR_LENGTH(alt) = 1   
               AND v.ref IN (ms.derived, ms.ancestral)
               AND v.alt IN (ms.derived, ms.ancestral)""".format(pop=self.population, chrom=self.chrom))

        with self.output().open('w') as fout:
            fout.write('Execution took {}'.format(exec_time))


class EnsemblPipeline(luigi.WrapperTask):
    """
    Populate the ensembl_* tables and link modern_snps records to genes, dbsnp variants.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):

        # process all the populations in chromosome chunks
        for pop in SAMPLES[self.species]:
            for chrom in CHROM_SIZE[self.species]:
                yield LinkEnsemblGenes(self.species, pop, chrom)
                yield LinkEnsemblVariants(self.species, pop, chrom)


if __name__ == '__main__':
    luigi.run()
