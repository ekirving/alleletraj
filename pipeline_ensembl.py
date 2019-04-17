#!/usr/bin/env python
# -*- coding: utf-8 -*-

import luigi
import gzip

from collections import Counter

# import my custom modules
from pipeline_consts import ENSEMBL_DATA, MAX_INSERT_SIZE  # TODO make these into PipelineTask properties
from pipeline_utils import PipelineTask, db_conn


class LoadEnsemblGenes(PipelineTask):
    """
    Load the GTF (General Transfer Format) data from Ensembl.

    GTF provides access to all annotated transcripts which make up an Ensembl gene set.

    See https://www.ensembl.org/info/website/upload/gff.html#fields

    :type species: str
    """
    species = luigi.Parameter()

    # TODO
    # def requires(self):
    #     return ProcessSNPs(self.species, self.population, self.chrom)

    # TODO
    # def output(self):
    #     return luigi.LocalTarget('db/{}-snpchip.log'.format(self.basename))

    def run(self):
        # open a db connection
        dbc = db_conn(self.species)

        print("INFO: Loading Ensembl genes")

        # the column headers for batch inserting into the db
        fields = ('source', 'gene_id', 'version', 'biotype', 'chrom', 'start', 'end')

        # open the GTF file
        with gzip.open(ENSEMBL_DATA[self.species]['gtf'], 'r') as fin:

            # buffer the records for bulk insert
            records = []

            for line in fin:
                # strip flanking whitespace
                line = line.strip()

                # skip comments
                if not line.startswith("#"):

                    # grab the GFF columns
                    chrom, source, feature, start, end, score, strand, frame, attributes = line.split("\t")

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

            # bulk insert all the reads for this sample
            if records:
                dbc.save_records('ensembl_genes', fields, records)


class LoadEnsemblVariants(PipelineTask):
    """
    Load the GVF (Genome Variation Format) data from Ensembl.

    Contains all germline variants from the current Ensembl release for this species.

    See http://www.sequenceontology.org/gvf.html

    :type species: str
    """
    species = luigi.Parameter()

    # TODO
    # def requires(self):
    #     return ProcessSNPs(self.species, self.population, self.chrom)

    # TODO
    # def output(self):
    #     return luigi.LocalTarget('db/{}-snpchip.log'.format(self.basename))

    def run(self):
        dbc = db_conn(self.species)

        print("INFO: Loading Ensembl variants")

        # the column headers for batch inserting into the db
        fields = ('dbxref', 'rsnumber', 'type', 'chrom', 'start', 'end', 'ref', 'alt')

        # open the GVF file
        with gzip.open(ENSEMBL_DATA[self.species]['gvf'], 'r') as fin:

            # buffer the records for bulk insert
            records = []

            for line in fin:
                # strip flanking whitespace
                line = line.strip()

                # skip comments
                if not line.startswith("#"):

                    # grab the GFF columns
                    chrom, source, type, start, end, score, strand, phase, attributes = line.split("\t")

                    # unpack the 'key1=value1;key2=value2;' attribute pairs
                    atts = dict([pair.strip().split('=') for pair in attributes.split(';') if pair != ''])

                    # unpack dbxref
                    dbxref, rsnumber = atts['Dbxref'].split(":")

                    # setup the record to insert, in this order
                    variant = (dbxref, rsnumber, type, chrom, start, end, atts['Reference_seq'], atts['Variant_seq'])

                    # add to the buffer
                    records.append(variant)

                    # bulk insert in chunks
                    if len(records) == MAX_INSERT_SIZE:
                        dbc.save_records('ensembl_variants', fields, records)
                        records = []

            # insert any the remaing variants
            if records:
                dbc.save_records('ensembl_variants', fields, records)


if __name__ == '__main__':
    luigi.run()
