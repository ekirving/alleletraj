#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
import os
import re
from collections import defaultdict, OrderedDict

# third party modules
import luigi
from natsort import natsorted

# local modules
from alleletraj import utils
from alleletraj.const import QTLDB_RELEASE, SWEEP_DATA
from alleletraj.db.load import CreateDatabase
from alleletraj.ensembl.load import LoadEnsemblVariants, LoadEnsemblGenes
from alleletraj.modern.snps import LoadModernSNPs
from alleletraj.qtl.qtldb_api import QTLdbAPI
from alleletraj.ref import ReferenceFASTA

# offset to use for the QTL window (+/- 50 Kb)
QTL_WINDOW = 50000

# offset all genes by 100 Kb to preclude linkage with our 'neutral' SNPs
GENE_OFFSET = 100000

# the minimum derived allele frequency of modern SNPs to include
MIN_DAF = 0.05

# Melanocortin 1 receptor
MC1R_GENE = 'MC1R'


def extract_qtl_fields(tsv_file, fields):
    """
    Helper function for extracting specified fields from the AnimalQTLdb data file.
    """
    data = defaultdict(list)

    # get all the QTL IDs for this scope
    with open(tsv_file, 'rU') as fin:

        # get the column headers
        header = fin.readline().split('\t')

        # get the index of the target fields
        columns = [(idx, field) for idx, field in enumerate(header) if field in fields]

        for line in fin:
            try:
                line = line.split('\t')
                for idx, field in columns:
                    datum = line[idx].strip()
                    if datum:
                        data[field].append(datum)
            except IndexError:
                # ignore badly formatted lines
                pass

    return data


class ExternalAnimalQTLdb(utils.PipelineExternalTask):
    """
    External task dependency for AnimalQTLdb data.

    See https://www.animalgenome.org/cgi-bin/QTLdb/index

    :type species: str
    """
    species = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget('data/qtldb/{}_cM_{}.txt'.format(self.species, QTLDB_RELEASE))


class PopulateQTLs(utils.MySQLTask):
    """
    Fetch all the QTLs from the QTLdb API and populate the local database.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        yield ExternalAnimalQTLdb(self.species)
        yield CreateDatabase(self.species)

    def queries(self):
        qtl_file, _ = self.input()

        api = QTLdbAPI()

        # get a list of all the QTLDb IDs from the tsv dump
        data = extract_qtl_fields(qtl_file.path, ['QTL_ID'])

        # convert all the IDs to int
        qtl_ids = [int(qtl_id) for qtl_id in data['QTL_ID'] if qtl_id.isdigit()]

        # get all the QTLs already in the DB
        qtls = self.dbc.get_records_sql("""
            SELECT qtldb_id 
              FROM qtls 
             WHERE qtldb_id IS NOT NULL""", key='qtldb_id')

        # find the new IDs in the list
        new_ids = list(set(qtl_ids) - set(qtls.keys()))

        # fields to rename
        key_map = {
            'id': 'qtldb_id',
            'pubmedID': 'pubmed_id',
            'geneId': 'gene_id',
            'chromosome': 'chrom'
        }

        # get all the new records
        for record in api.get_qtls(self.species, new_ids):

            # extract the nested trait record
            trait = record.pop('trait')
            trait['name'] = record.pop('name')

            # set the tait foreign key on the main record
            record['trait_id'] = trait['traitID']

            # does the trait exist
            if not self.dbc.exists_record('traits', {'id': record['trait_id']}):
                # setup the trait record
                trait = dict((field.replace('trait', '').lower(), trait[field]) for field in trait)
                # TODO this is broken!!
                trait['type'] = ''  # api.get_trait_type(self.species, trait['id'], trait['name'])

                self.dbc.save_record('traits', trait, insert=True)

            # does the publication exist
            if not self.dbc.exists_record('pubmeds', {'id': record['pubmedID']}):
                # setup the pubmed record
                pubmed = {}  # api.get_publication(self.species, record['pubmedID'])  # TODO this is broken!!

                if pubmed:
                    pubmed['id'] = pubmed.pop('pubmed_ID')
                    pubmed['year'] = re.search(r'\(([0-9]{4})\)', pubmed['authors']).group(1)
                    pubmed['journal'] = pubmed['journal']['#text'][:-5]

                    self.dbc.save_record('pubmeds', pubmed, insert=True)
                else:
                    # TODO some records have a bogus pubmed ID, but these appear to work on the website
                    record['pubmedID'] = None

            # flatten the other nested records
            for field in record:

                if type(record[field]) is OrderedDict:
                    nested = record.pop(field)

                    for name in nested:
                        # for doubly nested fields, use the parent name as a prefix
                        if type(nested[name]) is OrderedDict:
                            for key in nested[name]:
                                record[name + '_' + key] = nested[name][key]
                        elif field in ['gene']:
                            record[field + name.title()] = nested[name]
                        else:
                            record[name] = nested[name]

            # drop any lingering malformed fields
            record.pop('source', None)
            record.pop('breeds', None)
            record.pop('effects', None)  # TODO check this out
            record.pop('statTests', None)

            # handle malformed data
            for field in ['linkageLoc_end', 'linkageLoc_peak', 'linkageLoc_start']:
                if field in record and record[field] is not None:
                    record[field] = re.sub('[^0-9.]', '', record[field])

            # rename some fields
            for key in record:
                if key in key_map:
                    record[key_map[key]] = record.pop(key)

            # filter out any empty values
            qtl = OrderedDict((key, record[key]) for key in record if record[key] != '-')

            self.dbc.save_record('qtls', qtl, insert=True)


class SetQTLWindows(utils.MySQLTask):
    """
    Calculate the QTL window sizes.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        yield ReferenceFASTA(self.species)
        yield PopulateQTLs(self.species)
        yield LoadEnsemblVariants(self.species)

    def queries(self):
        # unpack the inputs
        (_, fai_file), _, _ = self.input()

        # get all the valid QTL windows
        results = self.dbc.get_records_sql("""
            SELECT q.id, v.chrom, v.start AS site
              FROM qtls q
              JOIN ensembl_variants v
                ON v.rsnumber = q.peak                # only QTLs with a dbsnp peak
               AND v.type = 'SNV'                     # only single nucleotide variants
               AND CHAR_LENGTH(v.alt) = 1             # only biallelic sites
             WHERE q.associationType = 'Association'  # only GWAS studies
               AND q.significance = 'Significant'     # only significant hits
               """, key=None)

        # get the sizes of the chromosomes, to bound the QTL windows
        sizes = utils.get_chrom_sizes(fai_file)

        for result in results:
            # get the size of the current chrom
            chrom_size = sizes[result['chrom']]

            # calculate the bounded window size
            start = result['site'] - QTL_WINDOW if result['site'] > QTL_WINDOW else 1
            end = result['site'] + QTL_WINDOW if result['site'] + QTL_WINDOW < chrom_size else chrom_size

            if end <= start:
                raise Exception(
                    'ERROR: Window size for QTL #{} is negative ({:,} bp)'.format(result['id'], end - start))

            # update the QTL record
            qtl = {
                'id': result['id'],
                'chrom': result['chrom'],
                'valid': 1,
                'site': result['site'],
                'start': start,
                'end': end
            }

            self.dbc.save_record('qtls', qtl)


class PopulateSweepLoci(utils.MySQLTask):
    """
    Populate the db with any selective sweep regions.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        yield CreateDatabase(self.species)

    def queries(self):
        # get the files containing the sweep data
        loci_file = SWEEP_DATA[self.species]['loci']  # TODO make SWEEP_DATA into external dependencies
        snps_file = SWEEP_DATA[self.species]['snps']

        with open(loci_file, 'r') as loci_fin:

            for locus in loci_fin:
                # extract the locus info from the BED file
                chrom, start, end = locus.split()

                # setup a dummy QTL record
                qtl = {
                    'associationType': 'Sweep',
                    'chrom': chrom,
                    'significance': 'Significant',
                    'valid': 1,
                    'start': start,
                    'end': end,
                }

                qtl_id = self.dbc.save_record('qtls', qtl)

                # get the all the SNPs from this locus
                snps = utils.run_cmd(["printf '{}' | bedtools intersect -a {} -b stdin"
                                     .format(locus.strip(), snps_file)], shell=True, verbose=False)

                if not snps:
                    raise Exception('ERROR: Found no SNPs for sweep region {}:{}-{}'.format(chrom, start, end))

                for snp in snps.splitlines():
                    # extract the SNP info
                    chrom, start, end, cdf, p = snp.split()

                    sweep_snp = {
                        'qtl_id': qtl_id,
                        'chrom': chrom,
                        'site': end,
                        'cdf': cdf,
                        'p': p
                    }

                    self.dbc.save_record('sweep_snps', sweep_snp)


class PopulateMC1RLocus(utils.MySQLTask):
    """
    Populate a dummy QTL for the MC1R gene.

    # TODO make this a generic method for genes
    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return LoadEnsemblGenes(self.species)

    def queries(self):
        # get the MC1R gene record
        mc1r = self.dbc.get_record('ensembl_genes', {'gene_name': MC1R_GENE})

        # setup a dummy QTL record
        qtl = {
            'associationType': 'MC1R',
            'chrom': mc1r['chrom'],
            'valid': 1,
            'start': mc1r['start'],
            'end': mc1r['end'],
        }

        self.dbc.save_record('qtls', qtl)


class PopulatePigMummyLoci(utils.MySQLTask):
    """
    Balancing selection on a recessive lethal deletion with pleiotropic effects on two neighboring genes in the porcine
    genome.

    # TODO why can't we make this a 'sweep' locus?

    See https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007661

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        yield ReferenceFASTA(self.species)
        yield LoadEnsemblVariants(self.species)

    # noinspection SqlResolve
    def queries(self):
        # unpack the inputs
        (_, fai_file), _ = self.input()

        # get the sizes of the chromosomes
        sizes = utils.get_chrom_sizes(fai_file)

        # compose a CASE statement to cap the upper bound of the QTLs by the size of the chromosome
        max_chrom = ' '.join(["WHEN '{}' THEN {}".format(chrom, sizes[chrom]) for chrom in sizes])

        # get all pig mummy SNPs and turn them into pseudo-loci
        results = self.dbc.get_records_sql("""
            SELECT ev.rsnumber, 
                   ev.chrom, 
                   GREATEST(ev.start - {offset}, 1) AS `start`,
                   LEAST(ev.end + {offset}, CASE ev.chrom {max_chrom} END) AS `end`
              FROM ensembl_variants ev
             WHERE rsnumber IN ('rs321688936','rs324480289','rs331589427','rs333966168','rs334806226','rs337374395',
                                'rs340056046','rs340481593','rs341362223','rs342358904','rs81238716','rs81255938',
                                'rs81270496','rs81317835','rs81469247','rs81469256','rs81469273','rs81469281',
                                'rs81469291', 'rs81469298','rs81469311','rs81469316','rs81469337','rs81469339',
                                'rs81469341','rs81469348')
                                """.format(offset=GENE_OFFSET, max_chrom=max_chrom), key=None)

        for result in results:
            # setup a dummy QTL record
            qtl = {
                'associationType': 'Mummy',
                'chrom': result['chrom'],
                'peak': result['rsnumber'],
                'valid': 1,
                'start': result['start'],
                'end': result['end'],
            }

            self.dbc.save_record('qtls', qtl)


class PopulateTraitLoci(utils.PipelineWrapperTask):
    """
    Wrapper task to populate all the QTLs and other trait loci of interest.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):

        # load the QTLs from the AnimalQTL db, and set the window size
        yield SetQTLWindows(self.species)

        # load pseudo-QTLs from other sources
        if self.species == 'pig':
            yield PopulateSweepLoci(self.species)

        yield PopulateMC1RLocus(self.species)

        if self.species == 'pig':
            yield PopulatePigMummyLoci(self.species)


class PopulateNeutralLoci(utils.MySQLTask):
    """
    Populate dummy QTLs for all the 'neutral' loci (i.e. regions outside of all QTLs and gene regions +/- a buffer).

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        yield ReferenceFASTA(self.species)
        yield PopulateTraitLoci(self.species)

    # noinspection SqlResolve
    def queries(self):
        # unpack the inputs
        (_, fai_file), _ = self.input()

        # get the sizes of the chromosomes
        sizes = utils.get_chrom_sizes(fai_file)

        # compose a CASE statement to cap the upper bound of the QTLs by the size of the chromosome
        max_chrom = ' '.join(["WHEN '{}' THEN {}".format(chrom, sizes[chrom]) for chrom in sizes])

        # get all the non-neutral regions (defined here as all valid QTLs and gene regions +/- offset)
        results = self.dbc.get_records_sql("""

            # get all the QTLs with a valid rsnumber
            SELECT ev.chrom,
                   GREATEST(ev.start - {offset}, 1) AS `start`,
                   LEAST(ev.end + {offset}, CASE ev.chrom {max_chrom} END) AS `end`
              FROM qtls q
              JOIN ensembl_variants ev
                ON ev.rsnumber = q.peak
             WHERE q.associationType = 'Association'

            UNION

            # and all the sweep regions
            SELECT q.chrom, 
                   GREATEST(q.start - {offset}, 1) AS `start`, 
                   LEAST(q.end + {offset}, CASE q.chrom {max_chrom} END) AS `end`
              FROM qtls q
             WHERE q.associationType = 'Sweep'

            UNION

            # and all the genes
            SELECT eg.chrom,
                   GREATEST(eg.start - {offset}, 1) AS `start`,
                   LEAST(eg.end + {offset}, CASE eg.chrom {max_chrom} END) AS `end`
              FROM ensembl_genes eg

          ORDER BY chrom, start, end
               """.format(offset=GENE_OFFSET, max_chrom=max_chrom), key=None)

        intervals = defaultdict(list)

        # group the intervals by chrom
        for result in results:
            intervals[result['chrom']].append((result['start'], result['end']))

        # make some temp files
        all_regions = luigi.LocalTarget(is_tmp=True)
        non_neutral = luigi.LocalTarget(is_tmp=True)

        # write a BED file for the whole genome
        with all_regions.open('w') as fout:
            for chrom in self.chromosomes:
                fout.write('{}\t{}\t{}\n'.format(chrom, 1, sizes[chrom]))

        # write all the non-neutral regions to a BED file
        with non_neutral.open('w') as fout:
            for chrom in natsorted(intervals.keys()):
                # merge overlapping intervals
                for start, stop in utils.merge_intervals(intervals[chrom], capped=False):
                    fout.write('{}\t{}\t{}\n'.format(chrom, start, stop))

        # subtract the non-neutral regions from the whole genome
        loci = utils.run_cmd(['bedtools', 'subtract', '-a', all_regions.path, '-b', non_neutral.path])

        for locus in loci.splitlines():
            chrom, start, end = locus.split()

            # setup a dummy QTL record
            qtl = {
                'associationType': 'Neutral',
                'chrom': chrom,
                'valid': 1,
                'start': start,
                'end': end,
            }

            self.dbc.save_record('qtls', qtl)


class PopulateAllLoci(utils.PipelineWrapperTask):
    """
    Wrapper task to populate all the QTL and pseudo-QTL windows.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        # load all the QTLs and neutral loci
        yield PopulateTraitLoci(self.species)
        yield PopulateNeutralLoci(self.species)


class MergeAllLoci(utils.DatabaseTask):
    """
    Merge all the overlapping QTL and pseudo-QTL windows.

    :type species: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        yield PopulateAllLoci(self.species)

    def output(self):
        return luigi.LocalTarget('data/bed/{}-loci.bed'.format(self.basename))

    def run(self):
        # get all the unique QTL windows
        loci = self.dbc.get_records_sql("""
            SELECT DISTINCT q.chrom, q.start, q.end
              FROM qtls q
             WHERE q.chrom = '{chrom}' 
               AND q.valid = 1
          ORDER BY q.start, q.end""".format(chrom=self.chrom), key=None)

        tmp_loci = 'data/bed/{}-tmp-loci.bed'.format(self.basename)

        # write all the QTL regions to a BED file
        with open(tmp_loci, 'w') as fout:
            for locus in loci:
                # NOTE that BED starts are zero-based and BED ends are one-based
                fout.write('{}\t{}\t{}\n'.format(locus['chrom'], locus['start'] - 1, locus['end']))

        # now merge overlapping loci
        bed = utils.run_cmd(['bedtools', 'merge', '-i', tmp_loci])

        # tidy up the tmp file
        os.remove(tmp_loci)

        # save the merged loci
        with self.output().open('w') as fout:
            fout.write(bed)


class PopulateQTLSNPs(utils.MySQLTask):
    """
    Now we have ascertained all the modern SNPs, let's find those that intersect with the QTLs.

    :type species: str
    :type chrom: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        yield PopulateAllLoci(self.species)
        yield LoadModernSNPs(self.species, self.chrom)

    # noinspection SqlResolve
    def queries(self):
        self.dbc.execute_sql("""
            INSERT INTO qtl_snps (qtl_id, modsnp_id)
            SELECT DISTINCT q.id, ms.id
              FROM qtls q
              JOIN modern_snps ms
                ON ms.chrom = q.chrom
               AND ms.site BETWEEN q.start AND q.end
              JOIN modern_snp_daf msd
                ON msd.modsnp_id = ms.id
             WHERE q.chrom = '{chrom}'
               AND q.valid = 1
               AND msd.daf >= {daf}
               """.format(chrom=self.chrom, daf=MIN_DAF))


class MarkNeutralSNPs(utils.MySQLTask):
    """
    Mark neutral SNPs (i.e. SNPs outside of all QTLs and gene regions).

    :type species: str
    :type chrom: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    db_lock_tables = ['modern_snps_{chrom}']

    def requires(self):
        return PopulateQTLSNPs(self.species, self.chrom)

    def queries(self):
        self.dbc.execute_sql("""
            UPDATE modern_snps ms
              JOIN qtl_snps qs
                ON qs.modsnp_id = ms.id
              JOIN qtls q
                ON q.id = qs.qtl_id
               SET ms.neutral = 1
             WHERE q.chrom = '{chrom}'
               AND q.associationType = 'Neutral'
               AND q.valid = 1
               """.format(chrom=self.chrom))


class QTLPipeline(utils.PipelineWrapperTask):
    """
    Load all the QTLs, pseudo-QTLs and Neutral regions.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        # process all the chromosomes in chunks
        for chrom in self.chromosomes:
            # flag the modern SNPs which fall into 'neutral' regions
            yield MarkNeutralSNPs(self.species, chrom)


if __name__ == '__main__':
    luigi.run()
