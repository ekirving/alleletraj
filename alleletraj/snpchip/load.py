#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
import os
import random

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.db.load import CreateDatabase
from alleletraj.ensembl.load import LoadEnsemblVariants

AXIOM_URL = 'http://media.affymetrix.com/analysis/downloads/lf/genotyping/Axiom_MNEc670/r2/' \
            'Axiom_MNEc670_Annotation.r2.csv.zip'


class ExternalSNPchimp(utils.PipelineExternalTask):
    """
    External task dependency for SNPchimp data.

    NOTE: We use the rsnumber to anchor the snpchip to the assembly, not the chrom/pos details in the gzip file.

    See http://bioinformatics.tecnoparco.org/SNPchimp

    :type species: str
    """
    species = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget('data/snpchip/SNPchimp_{}.tsv.gz'.format(self.species))


class DownloadAxiomEquineHD(utils.PipelineTask):
    """
    Fetches the data for the Axiom EquineHD array.
    """

    @property
    def url(self):
        return AXIOM_URL

    def output(self):
        return luigi.LocalTarget('data/snpchip/{}'.format(os.path.basename(self.url)))

    def run(self):
        with self.output().temporary_path() as tmp_path:
            utils.curl_download(self.url, tmp_path)


class LoadSNPChipVariants(utils.DatabaseTask):
    """
    Load the SNP chip variants from SNPchimp

    See http://bioinformatics.tecnoparco.org/SNPchimp

    :type species: str
    """
    species = luigi.Parameter()

    db_lock_tables = ['snpchip']

    def requires(self):
        yield ExternalSNPchimp(self.species)
        yield CreateDatabase(self.species)

    def output(self):
        return luigi.LocalTarget('data/db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        # get the input file
        gzip_file, _ = self.input()

        # unzip the archive into a named pipe
        tmp_file = luigi.LocalTarget(is_tmp=True)
        utils.run_cmd(['gzip --stdout -d  {gz} > {tmp}'.format(gz=gzip_file.path, tmp=tmp_file.path)], shell=True,
                      background=True)

        # load the data into the db
        self.dbc.execute_sql("""
            LOAD DATA 
         LOCAL INFILE '{tmp}'
           INTO TABLE snpchip 
               IGNORE 1 LINES (chip_name, rsnumber, @dummy, @dummy, snp_name)
                  """.format(tmp=tmp_file.path))

        # tidy up NULL values which get imported as the string 'NULL'
        self.dbc.execute_sql("""
            UPDATE snpchip
               SET rsnumber = NULL
             WHERE rsnumber = 'NULL'""")

        with self.output().open('w') as fout:
            fout.write('Loaded SNPchimp records')


class LoadAxiomEquineHD(utils.DatabaseTask):
    """
    SNPchimp doesn't have the details for the Affymetrix Axiom EquineHD array, so we have to load the data separately.

    :type species: str
    """
    species = luigi.Parameter()

    db_lock_tables = ['snpchip']

    def requires(self):
        yield DownloadAxiomEquineHD()
        yield LoadSNPChipVariants(self.species)
        yield LoadEnsemblVariants(self.species)

    def output(self):
        return luigi.LocalTarget('data/db/{}-{}.log'.format(self.basename, self.classname))

    # noinspection SqlWithoutWhere
    def run(self):
        # get the input file
        axiom_file, _, _ = self.input()

        # just get the relevant columns
        awk = "awk -F ',' 'NR>1 {print $1 \"\\t\" $4 \"\\t\" $5}'"

        # unzip the dump into a temp file
        tmp_file = luigi.LocalTarget(is_tmp=True)
        utils.run_cmd(["unzip -p {axiom} | grep -vP '^#' | {awk} > {tmp}"
                      .format(axiom=axiom_file.path, awk=awk, tmp=tmp_file.path)], shell=True, background=True)

        # load the data into the db
        self.dbc.execute_sql("""
                LOAD DATA 
             LOCAL INFILE '{tmp}'
               INTO TABLE snpchip_axiom
                   FIELDS 
              ENCLOSED BY '"' (snp_name, chrom, site)""".format(tmp=tmp_file.path))

        # fix the missing rsnumbers
        self.dbc.execute_sql("""
            UPDATE snpchip sc
              JOIN snpchip_axiom sa
                ON sa.snp_name = sc.snp_name
              JOIN ensembl_variants ev
                ON sa.chrom = ev.chrom 
               AND sa.site = ev.start
               SET sc.rsnumber = ev.rsnumber
             WHERE ev.type = 'SNV'""")

        with self.output().open('w') as fout:
            fout.write('Loaded Axiom EquineHD records')


class SNPChipLoadPipeline(utils.PipelineWrapperTask):
    """
    Populate the snpchip tables.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        yield LoadSNPChipVariants(self.species)

        if self.species == 'horse':
            # handle special case for horse data
            yield LoadAxiomEquineHD(self.species)


if __name__ == '__main__':
    luigi.run()
