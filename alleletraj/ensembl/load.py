#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
import gzip
import os

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.db.load import CreateDatabase

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
            if ENSEMBL_RELEASES[self.assembly] >= 91:
                return 'ftp://ftp.ensembl.org/pub/release-{rel}/variation/gvf/{bin}/{bin}.gvf.gz'.format(**self.params)
            else:
                return 'ftp://ftp.ensembl.org/pub/release-{rel}/variation/gvf/{bin}/{Bin}.gvf.gz'.format(**self.params)

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


class LoadEnsemblGenes(utils.MySQLTask):
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

    def queries(self):
        # unpack the params
        _, gtf_file = self.input()

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
            self.dbc.save_records('ensembl_genes', fields, records)


class LoadEnsemblVariants(utils.MySQLTask):
    """
    Load the GVF (Genome Variation Format) data from Ensembl.

    Contains all germline variants from the current Ensembl release for this species.

    See http://www.sequenceontology.org/gvf.html

    :type species: str
    """
    species = luigi.Parameter()

    db_lock_tables = ['ensembl_variants']

    def requires(self):
        yield DownloadEnsemblData(self.species, 'gvf')
        yield CreateDatabase(self.species)

    def queries(self):
        # unpack the params
        gvf_file, _ = self.input()

        # make a temp file to buffer the input data
        tmp_file = luigi.LocalTarget(is_tmp=True)

        with gzip.open(gvf_file.path, 'r') as fin, tmp_file.open('w') as fout:

            for line in fin:
                if line.startswith('#'):
                    # skip comments
                    continue

                # unpack the GFF columns
                chrom, source, snptype, start, end, _, _, _, attributes = line.strip().split('\t')

                # split the 'key1=value1;key2=value2;' attribute pairs
                atts = dict([pair.strip().split('=') for pair in attributes.split(';') if pair != ''])

                # split the 'key=value1:value2' pair
                dbxref, rsnumber = atts['Dbxref'].split(':')

                record = [dbxref, rsnumber, snptype, chrom, start, end, atts['Reference_seq'], atts['Variant_seq']]

                # save to temp file
                fout.write('\t'.join(record) + '\n')

        # load the temp file
        self.dbc.execute_sql("""
            LOAD DATA
         LOCAL INFILE '{tmp}'
           INTO TABLE ensembl_variants (dbxref, rsnumber, type, chrom, start, end, ref, alt)
                  """.format(tmp=tmp_file.path))


class FlagSNPsNearIndels(utils.MySQLTask):
    """
    Flag any dbsnp variants that are within a given range of an INDEL.

    :type species: str
    :type chrom: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    db_lock_tables = ['ensembl_variants_{chrom}']

    def requires(self):
        return LoadEnsemblVariants(self.species)

    # noinspection SqlResolve
    def queries(self):
        indels = self.dbc.get_records_sql("""
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
        for i in range(0, len(loci), self.dbc.max_query_size):
            # convert each locus into sql conditions
            conds = ["start BETWEEN {} AND {}".format(start, end) for start, end in loci[i:i + self.dbc.max_query_size]]

            self.dbc.execute_sql("""
                UPDATE ensembl_variants
                   SET indel = 1
                 WHERE chrom = '{chrom}'
                   AND ({conds})
                   AND type = 'SNV'
                   """.format(chrom=self.chrom, conds=" OR ".join(conds)))


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
