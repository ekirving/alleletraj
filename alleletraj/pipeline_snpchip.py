#!/usr/bin/env python
# -*- coding: utf-8 -*-

import luigi
import os
import random

# import my custom modules
from alleletraj.database.setup import CreateDatabase
from alleletraj.pipeline_ensembl import LoadEnsemblVariants
from alleletraj.modern.load_snps import LoadModernSNPs
from alleletraj.utils import PipelineTask, PipelineExternalTask, PipelineWrapperTask, run_cmd, curl_download, trim_ext

AXIOM_URL = 'http://media.affymetrix.com/analysis/downloads/lf/genotyping/Axiom_MNEc670/r2/' \
            'Axiom_MNEc670_Annotation.r2.csv.zip'


class ExternalSNPchimp(PipelineExternalTask):
    """
    External task dependency for SNPchimp data.

    See http://bioinformatics.tecnoparco.org/SNPchimp

    :type species: str
    """
    species = luigi.Parameter()

    def output(self):
        # TODO goat is using the wrong assembly, but that might not matter because it has no chrom-pos entries
        return luigi.LocalTarget('snpchip/SNPchimp_{}.tsv.gz'.format(self.assembly))


class DownloadAxiomEquineHD(PipelineTask):
    """
    Fetches the data for the Axiom EquineHD array.
    """
    @property
    def url(self):
        return AXIOM_URL

    def output(self):
        return luigi.LocalTarget('snpchip/{}'.format(os.path.basename(self.url)))

    def run(self):
        with self.output().temporary_path() as tmp_path:
            curl_download(self.url, tmp_path)


class LoadSNPChipVariants(PipelineTask):
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
        return luigi.LocalTarget('db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        # get the input file
        gzip_file, _ = self.input()

        # open a db connection
        dbc = self.db_conn()

        # unzip the archive into a named pipe
        pipe = '{}-luigi-tmp-{:010}'.format(trim_ext(gzip_file.path), random.randrange(0, 1e10))
        run_cmd(['mkfifo -m0666 {pipe}'.format(pipe=pipe)], shell=True)
        run_cmd(['gzip --stdout -d  {gz} > {pipe}'.format(gz=gzip_file.path, pipe=pipe)], shell=True, background=True)

        # load the data into the db
        dbc.execute_sql("""
            LOAD DATA 
         LOCAL INFILE '{pipe}'
           INTO TABLE snpchip 
               IGNORE 1 LINES (chip_name, rsnumber, chrom, site, snp_name)
                  """.format(pipe=pipe))

        # tidy up NULL values which get imported as the string 'NULL'
        dbc.execute_sql("""
            UPDATE snpchip
               SET rsnumber = NULL
             WHERE rsnumber = 'NULL'""")

        # remove the named pipe
        run_cmd(['rm -f {pipe}'.format(pipe=pipe)], shell=True)

        with self.output().open('w') as fout:
            fout.write('Loaded SNPchimp records')


class LoadAxiomEquineHD(PipelineTask):
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
        return luigi.LocalTarget('db/{}-{}.log'.format(self.basename, self.classname))

    # noinspection SqlWithoutWhere
    def run(self):
        # get the input file
        axiom_file, _, _ = self.input()

        # open a db connection
        dbc = self.db_conn()

        # just get the relevant columns
        awk = "awk -F ',' 'NR>1 {print $1 \"\\t\" $4 \"\\t\" $5}'"

        # unzip the dump into a named pipe
        pipe = '{}-luigi-tmp-{:010}'.format(trim_ext(axiom_file.path), random.randrange(0, 1e10))
        run_cmd(['mkfifo -m0666 {pipe}'.format(pipe=pipe)], shell=True)
        run_cmd(["unzip -p {axiom} | grep -vP '^#' | {awk} > {pipe}".format(axiom=axiom_file.path, awk=awk, pipe=pipe)],
                shell=True, background=True)

        # load the data into the db
        dbc.execute_sql("""
                LOAD DATA 
             LOCAL INFILE '{pipe}'
               INTO TABLE snpchip_axiom
                   FIELDS 
              ENCLOSED BY '"' (snp_name, chrom, site)""".format(pipe=pipe))

        # remove the named pipe
        run_cmd(['rm -f {pipe}'.format(pipe=pipe)], shell=True)

        # fix the missing chrom/site data
        dbc.execute_sql("""
            UPDATE snpchip sc
              JOIN snpchip_axiom sa
                ON sa.snp_name = sc.snp_name
               SET sc.chrom = sa.chrom,
                   sc.site = sa.site""")

        # fix the missing rsnumbers
        dbc.execute_sql("""
            UPDATE snpchip sc
              JOIN ensembl_variants ev
                ON ev.chrom = sc.chrom
               AND ev.start = sc.site
               SET sc.rsnumber = ev.rsnumber
             WHERE ev.type = 'SNV'""")

        with self.output().open('w') as fout:
            fout.write('Loaded Axiom EquineHD records')


class LinkSNPChipVariants(PipelineTask):
    """
    Link modern SNPs to their SNPchip variants

    :type species: str
    :type chrom: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    db_lock_tables = ['modern_snps_{chrom}']

    def requires(self):
        yield LoadSNPChipVariants(self.species)

        if self.species == 'horse':
            # handle special case for horse data
            yield LoadAxiomEquineHD(self.species)

        yield LoadModernSNPs(self.species, self.chrom)

    def output(self):
        return luigi.LocalTarget('db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        dbc = self.db_conn()

        exec_time = dbc.execute_sql("""
            UPDATE modern_snps ms
              JOIN snpchip sc
                ON sc.chrom = ms.chrom
               AND sc.site = ms.site
               SET ms.snpchip_id = sc.id
             WHERE ms.chrom = '{chrom}'""".format(chrom=self.chrom))

        with self.output().open('w') as fout:
            fout.write('Execution took {}'.format(exec_time))


class SNPChipPipeline(PipelineWrapperTask):
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
