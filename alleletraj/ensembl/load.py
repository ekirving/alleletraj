#!/usr/bin/env python
# -*- coding: utf-8 -*-

import luigi
import os
import gzip

from datetime import timedelta
from time import time

# import my custom modules
from alleletraj.database.create import CreateDatabase
from alleletraj import utils

# the most recent Ensembl releases for a given genome assembly
ENSEMBL_RELEASES = {
    # cattle -  https://www.ensembl.org/Bos_taurus/Info/Index
    'UMD3.1': 94,       # Ensembl release 94 - October 2018
    'ARS-UCD1.2': 96,   # Ensembl release 96 - April 2019

    # goat -  https://www.ensembl.org/Capra_hircus/Info/Index
    'ARS1': 96,         # Ensembl release 96 - April 2019

    # horse - https://www.ensembl.org/Equus_caballus/Info/Index
    'EquCab2': 94,      # Ensembl release 94 - October 2018
    'EquCab3.0': 96,    # Ensembl release 96 - April 2019

    # pig - https://www.ensembl.org/Sus_scrofa/Info/Index
    'Sscrofa10.2': 89,  # Ensembl release 89 - May 2017
    'Sscrofa11.1': 96,  # Ensembl release 96 - April 2019
}

# minimum distance +/- from an INDEL
INDEL_BUFFER = 10


class DownloadEnsemblData(utils.PipelineTask):
    """
    Fetches the data from the Ensembl FTP site.

    :type species: str
    :type type: str
    :type bgzip: bool
    """
    species = luigi.Parameter()
    type = luigi.Parameter()
    bgzip = luigi.BoolParameter(default=False)

    @property
    def params(self):
        return {
            'rel': ENSEMBL_RELEASES[self.assembly],
            'Bin': self.binomial,
            'bin': self.binomial.lower(),
            'ref': self.assembly,
        }

    @property
    def url(self):
        if self.type == 'fasta':
            return 'ftp://ftp.ensembl.org/pub/release-{rel}/fasta/{bin}/dna/{Bin}.{ref}.dna.toplevel.fa.gz'.format(
                **self.params)
        elif self.type == 'gtf':
            return 'ftp://ftp.ensembl.org/pub/release-{rel}/gtf/{bin}/{Bin}.{ref}.{rel}.gtf.gz'.format(**self.params)
        elif self.type == 'gvf':
            return 'ftp://ftp.ensembl.org/pub/release-{rel}/variation/gvf/{bin}/{bin}.gvf.gz'.format(**self.params)

    def output(self):
        # handle non-unique gvf filenames
        filename = '{Bin}.{ref}.{rel}.gvf.gz'.format(**self.params) \
            if self.type == 'gvf' else os.path.basename(self.url)

        return luigi.LocalTarget('data/ensembl/{}'.format(filename))

    def run(self):
        with self.output().temporary_path() as tmp_path:
            utils.curl_download(self.url, tmp_path)

            if self.bgzip:
                # convert from gzip to bgzip (so we can index the file)
                utils.run_cmd(["mv {gz} {gz}-tmp; gunzip -c {gz}-tmp | bgzip > {gz}; rm {gz}-tmp".format(gz=tmp_path)],
                              shell=True)


class LoadEnsemblGenes(utils.PipelineTask):
    """
    Load the GTF (General Transfer Format) data from Ensembl.

    GTF provides access to all annotated transcripts which make up an Ensembl gene set.

    See https://www.ensembl.org/info/website/upload/gff.html#fields

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        yield CreateDatabase(self.species)
        yield DownloadEnsemblData(self.species, 'gtf')

    def output(self):
        return luigi.LocalTarget('data/db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        # unpack the inputs
        _, gtf_file = self.input()

        # open a db connection
        dbc = self.db_conn()

        # the column headers for batch inserting into the db
        fields = ('source', 'gene_id', 'gene_name', 'version', 'biotype', 'chrom', 'start', 'end')

        # open the GTF file
        with gzip.open(gtf_file.path, 'r') as fin:

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
                        gene = (atts['gene_source'], atts['gene_id'], atts.get('gene_name', ''), atts['gene_version'],
                                atts['gene_biotype'], chrom, start, end)

                        # add to the buffer
                        records.append(gene)

            # bulk insert all the records
            dbc.save_records('ensembl_genes', fields, records)

        with self.output().open('w') as fout:
            fout.write('Inserted {:,} Ensembl gene records'.format(len(records)))


class LoadEnsemblVariants(utils.PipelineTask):
    """
    Load the GVF (Genome Variation Format) data from Ensembl.

    Contains all germline variants from the current Ensembl release for this species.

    See http://www.sequenceontology.org/gvf.html

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        yield CreateDatabase(self.species)
        yield DownloadEnsemblData(self.species, 'gvf')

    def output(self):
        return luigi.LocalTarget('data/db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        # unpack the inputs
        _, gvf_file = self.input()

        # open a db connection
        dbc = self.db_conn()

        # the column headers for batch inserting into the db
        fields = ('dbxref', 'rsnumber', 'type', 'chrom', 'start', 'end', 'ref', 'alt')

        # open the GVF file
        with gzip.open(gvf_file.path, 'r') as fin:

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
                    if num_recs % dbc.max_insert_size == 0:
                        dbc.save_records('ensembl_variants', fields, records)
                        records = []

            # insert any remaining records
            if records:
                dbc.save_records('ensembl_variants', fields, records)

        with self.output().open('w') as fout:
            fout.write('Inserted {:,} Ensembl variant records'.format(num_recs))


class FlagSNPsNearIndels(utils.PipelineTask):
    """
    Flag any dbsnp variants that are within a given range of an INDEL

    :type species: str
    :type chrom: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    db_lock_tables = ['ensembl_variants_{chrom}']

    def requires(self):
        return LoadEnsemblVariants(self.species)

    def output(self):
        return luigi.LocalTarget('data/db/{}-{}.log'.format(self.basename, self.classname))

    # noinspection SqlResolve
    def run(self):
        dbc = self.db_conn()

        start = time()

        indels = dbc.get_records_sql("""
            SELECT ev.start, ev.end
              FROM ensembl_variants ev
             WHERE ev.chrom = '{chrom}' 
               AND ev.type IN ('insertion', 'deletion')
               """.format(chrom=self.chrom), key=None)

        loci = []

        for indel in indels:
            loci.append((int(indel['start']) - INDEL_BUFFER, int(indel['end']) + INDEL_BUFFER))

        # merge overlapping loci
        loci = list(utils.merge_intervals(loci))

        # process the INDELs in chunks
        for i in range(0, len(loci), dbc.max_query_size):

            # convert each locus into sql conditions
            conds = ["start BETWEEN {} AND {}".format(start, end) for start, end in loci[i:i + dbc.max_query_size]]

            dbc.execute_sql("""
                UPDATE ensembl_variants
                   SET indel = 1
                 WHERE chrom = '{chrom}'
                   AND ({conds})
                   AND type = 'SNV'
                   """.format(chrom=self.chrom, conds=" OR ".join(conds)))

        with self.output().open('w') as fout:
            fout.write('Execution took {}'.format(timedelta(seconds=time() - start)))


class EnsemblLoadPipeline(utils.PipelineWrapperTask):
    """
    Populate the ensembl genes and variants tables.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        yield LoadEnsemblGenes(self.species)

        for chrom in self.chromosomes:
            yield FlagSNPsNearIndels(self.species, chrom)


if __name__ == '__main__':
    luigi.run()
