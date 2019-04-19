#!/usr/bin/env python
# -*- coding: utf-8 -*-

import luigi
import re

from natsort import natsorted
from collections import defaultdict, OrderedDict

# import my custom modules
# TODO make these into PipelineTask properties
from pipeline_consts import CHROM_SIZE, MIN_DAF
from pipeline_ensembl import LoadEnsemblVariants, LoadEnsemblGenes
from pipeline_utils import PipelineTask, PipelineExternalTask, PipelineWrapperTask, run_cmd, merge_intervals
from dbconn import DBConn

from qtldb_api import QTLdbAPI

# QTLdb settings
QTLDB_RELEASE = 'rel37'

# offset to use for the QTL window (+/- 50 Kb)
QTL_WINDOW = 50000

# offset all genes by 100 Kb to preclude linkage with our 'neutral' SNPs
GENE_OFFSET = 100000

# genomic regions of selective sweeps ascertained in other papers
SWEEP_DATA = {
    # 'cattle': {}, # TODO add other species
    # 'goat': {},
    # 'pig': {},

    # see https://www.nature.com/articles/ng.3394
    'pig': {
        'loci': 'data/sweep/EUD_Sweep_p001_FINAL_cutoff_MERGED10kb.bed',
        'snps': 'data/sweep/EUD_Sweep_p001_FINAL_cutoff.bed'
    }
}

# the Ensembl gene ID
MC1R_GENE_ID = {
    'cattle': 'ENSBTAG00000023731',  # see http://www.ensembl.org/Bos_taurus/Gene/Summary?g=ENSBTAG00000023731
    'goat':   'ENSCHIG00000010476',  # see http://www.ensembl.org/Capra_hircus/Gene/Summary?db=core;g=ENSCHIG00000010476
    'pig':    'ENSSSCG00000020924',  # see http://www.ensembl.org/Sus_scrofa/Gene/Summary?g=ENSSSCG00000020924
    'horse':  'ENSECAG00000000900',  # see http://www.ensembl.org/Equus_caballus/Gene/Summary?g=ENSECAG00000000900
}


def extract_qtl_fields(dbfile, fields):
    """
    Helper function for extracting specified fields from the AnimalQTLdb data file.
    """
    data = defaultdict(list)

    # get all the QTL IDs for this scope
    with open(dbfile, 'rU') as fin:

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


class ExternalAnimalQTLdb(PipelineExternalTask):
    """
    External task dependency for AnimalQTLdb data.

    See https://www.animalgenome.org/cgi-bin/QTLdb/index

    :type species: str
    """
    species = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget('qtldb/{}_cM_{}.txt'.format(self.species, QTLDB_RELEASE))


class PopulateQTLs(PipelineTask):
    """
    Fetch all the QTLs from the QTLdb API and populate the local database.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return ExternalAnimalQTLdb(self.species)

    def output(self):
        return luigi.LocalTarget('db/{}-qtls.log'.format(self.species))

    def run(self):
        # get the file containing all the QTL IDs
        qtl_file = self.input()

        dbc = DBConn(self.species)
        api = QTLdbAPI()

        # get a list of all the QTLDb IDs from the tsv dump
        data = extract_qtl_fields(qtl_file.path, ['QTL_ID'])

        # convert all the IDs to int
        qtl_ids = [int(qtl_id) for qtl_id in data['QTL_ID'] if qtl_id.isdigit()]

        with self.output().open('w') as fout:
            fout.write("INFO: Processing {:,} QTLs from '{}'".format(len(qtl_ids), qtl_file.path))

            # get all the QTLs already in the DB
            qtls = dbc.get_records('qtls')

            # find the new IDs in the list
            new_ids = list(set(qtl_ids) - set(qtls.keys()))

            fout.write('INFO: Found {:,} new QTLs to add'.format(len(new_ids)))

            # rename these fields
            key_map = {
                'pubmedID':   'pubmed_id',
                'geneId':     'gene_id',
                'chromosome': 'chrom'
            }

            added = 0

            # get all the new records
            for record in api.get_qtls(self.species, new_ids):

                # TODO when resultset is len() = 1 then this throws an error
                # extract the nested trait record
                trait = record.pop('trait')
                trait['name'] = record.pop('name')

                # set the tait foreign key on the main record
                record['trait_id'] = trait['traitID']

                # does the trait exist
                if not dbc.exists_record('traits', {'id': record['trait_id']}):

                    # setup the trait record
                    trait = dict((field.replace('trait', '').lower(), trait[field]) for field in trait)
                    trait['type'] = ''  # api.get_trait_type(self.species, trait['id'], trait['name'])  # TODO this is broken!!

                    dbc.save_record('traits', trait, insert=True)

                # does the publication exist
                if not dbc.exists_record('pubmeds', {'id': record['pubmedID']}):

                    # setup the pubmed record
                    pubmed = False  # api.get_publication(self.species, record['pubmedID'])  # TODO this is broken!!

                    if pubmed:
                        pubmed['id'] = pubmed.pop('pubmed_ID')
                        pubmed['year'] = re.search('\(([0-9]{4})\)', pubmed['authors']).group(1)
                        pubmed['journal'] = pubmed['journal']['#text'][:-5]

                        dbc.save_record('pubmeds', pubmed, insert=True)
                    else:
                        # TODO some records have a bogus pubmed ID, but these appear to work on the website
                        record['pubmedID'] = None

                # flatten the other nested records
                for field, value in record.iteritems():

                    if type(value) is OrderedDict:
                        nested = record.pop(field)

                        for nested_name, nested_value in nested.iteritems():
                            # for doubly nested fields, use the parent name as a prefix
                            if type(nested_value) is OrderedDict:
                                for key in nested_value:
                                    record[nested_name + '_' + key] = nested_value[key]
                            elif field in ['gene']:
                                record[field + nested_name.title()] = nested_value
                            else:
                                record[nested_name] = nested_value

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
                qtl = OrderedDict((key, value) for key, value in record.iteritems() if value != '-')

                dbc.save_record('qtls', qtl, insert=True)

                added += 1

                fout.write('INFO: Added {:5d} new QTLs'.format(added))

            fout.write('INFO: Finished adding {} new QTLs'.format(len(new_ids)))


class SetQTLWindows(PipelineTask):
    """
    Calculate the QTL window sizes.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        yield PopulateQTLs(self.species)
        yield LoadEnsemblVariants(self.species)

    def output(self):
        return luigi.LocalTarget('db/{}-qtl_windows.log'.format(self.species))

    def run(self):
        # open a db connection
        dbc = DBConn(self.species)

        # get all the valid QTL windows
        results = dbc.get_records_sql("""
            SELECT q.id, v.chrom, v.start AS site
              FROM qtls q
              JOIN ensembl_variants v
                ON v.rsnumber = q.peak                # only QTLs with a dbsnp peak
               AND v.type = 'SNV'                     # only single nucleotide variants
               AND CHAR_LENGTH(v.alt) = 1             # only biallelic sites
             WHERE q.associationType = 'Association'  # only GWAS studies
               AND q.significance = 'Significant'     # only significant hits
               """, key=None)

        for result in results:
            # get the size of the current chrom
            chom_size = CHROM_SIZE[self.assembly][result['chrom']]

            # calculate the bounded window size
            start = result['site'] - QTL_WINDOW if result['site'] > QTL_WINDOW else 1
            end = result['site'] + QTL_WINDOW if result['site'] + QTL_WINDOW < chom_size else chom_size

            if end <= start:
                raise Exception('ERROR: Window size for QTL #{} is negative ({:,} bp)'.format(result['id'], end-start))

            # update the QTL record
            qtl = {
                'id': result['id'],
                'chrom': result['chrom'],
                'valid': 1,
                'site': result['site'],
                'start': start,
                'end': end
            }

            dbc.save_record('qtls', qtl)

        with self.output().open('w') as fout:
            fout.write('INFO: Set window sizes for {:,} QTLs'.format(len(results)))


class PopulateSweepLoci(PipelineTask):
    """
    Populate the db with any selective sweep regions.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        pass  # return BlahBlah(self.species) # TODO fix me

    def output(self):
        return luigi.LocalTarget('db/{}-sweep_loci.log'.format(self.species))

    def run(self):
        dbc = DBConn(self.species)

        # get the files containing the sweep data
        loci_file = SWEEP_DATA[self.species]['loci']  # TODO make SWEEP_DATA into external dependencies
        snps_file = SWEEP_DATA[self.species]['snps']

        num_loci = 0
        num_snps = 0

        with open(loci_file, 'r') as loci_fin:

            for locus in loci_fin:
                # extract the locus info from the BED file
                chrom, start, end = locus.split()

                # setup a dummy QTL record
                qtl = {
                    'associationType': 'Sweep',
                    'chrom':           chrom,
                    'significance':    'Significant',
                    'valid':           1,
                    'start':           start,
                    'end':             end,
                }

                qtl_id = dbc.save_record('qtls', qtl)

                num_loci += 1

                # get the all the SNPs from this locus
                snps = run_cmd(["printf '{}' | bedtools intersect -a {} -b stdin".format(locus.strip(), snps_file)],
                               shell=True)

                if not snps:
                    raise Exception('ERROR: Found no SNPs for sweep region {}:{}-{}'.format(chrom, start, end))

                for snp in snps.splitlines():
                    # extract the SNP info
                    chrom, start, end, cdf, p = snp.split()

                    sweep_snp = {
                        'qtl_id': qtl_id,
                        'chrom':  chrom,
                        'site':   end,
                        'cdf':    cdf,
                        'p':      p
                    }

                    dbc.save_record('sweep_snps', sweep_snp)

                    num_snps += 1

        with self.output().open('w') as fout:
            fout.write('INFO: Loaded {} selective sweep loci (inc. {} SNPs)'.format(num_loci, num_snps))


class PopulateMC1RLocus(PipelineTask):
    """
    Populate a dummy QTL for the MC1R gene.

    # TODO we want to model ALL 18 dbsnp variants, not just the 8 which are variable in EUD

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return LoadEnsemblGenes(self.species)

    def output(self):
        return luigi.LocalTarget('db/{}-mc1r_locus.log'.format(self.species))

    def run(self):
        dbc = DBConn(self.species)

        # get the MC1R gene details
        mc1r = dbc.get_record('ensembl_genes', {'gene_id': MC1R_GENE_ID[self.species]})

        # setup a dummy QTL record
        qtl = {
            'associationType': 'MC1R',
            'chrom': mc1r['chrom'],
            'valid': 1,
            'start': mc1r['start'],
            'end': mc1r['end'],
        }

        dbc.save_record('qtls', qtl)

        with self.output().open('w') as fout:
            fout.write('INFO: Added the MC1R gene locus')


class PopulatePigMummyLoci(PipelineTask):
    """
    Balancing selection on a recessive lethal deletion with pleiotropic effects on two neighboring genes in the porcine
    genome.

    # TODO why can't we make this a 'sweep' locus?

    See https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007661

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return LoadEnsemblVariants(self.species)

    def output(self):
        return luigi.LocalTarget('db/{}-mummy_loci.log'.format(self.species))

    def run(self):
        dbc = DBConn(self.species)

        sizes = CHROM_SIZE[self.assembly]

        # compose a CASE statement to cap the upper bound of the QTLs by the size of the chromosome
        max_chrom = ' '.join(["WHEN '{}' THEN {}".format(chrom, sizes[chrom]) for chrom in sizes])

        # get all pig mummy SNPs and turn them into pseudo-loci
        results = dbc.get_records_sql("""
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
                'peak':  result['rsnumber'],
                'valid': 1,
                'start': result['start'],
                'end':   result['end'],
            }

            dbc.save_record('qtls', qtl)

        with self.output().open('w') as fout:
            fout.write('INFO: Added {:,} pig mummy loci'.format(len(results)))


class PopulateTraitLoci(PipelineWrapperTask):
    """
    Wrapper task to populate all the QTLs and other trait loci of interest

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):

        # load the QTLs from the AnimalQTL database
        yield PopulateQTLs(self.species)
        yield SetQTLWindows(self.species)

        # load pseudo-QTLs from other sources
        yield PopulateSweepLoci(self.species)
        yield PopulateMC1RLocus(self.species)

        if self.species == 'pig':
            yield PopulatePigMummyLoci(self.species)


class PopulateNeutralLoci(PipelineTask):
    """
    Populate dummy QTLs for all the 'neutral' loci (i.e. regions outside of all QTLs and gene regions +/- a buffer)

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return PopulateTraitLoci(self.species)

    def output(self):
        return luigi.LocalTarget('db/{}-neutral_loci.log'.format(self.species))

    def run(self):
        dbc = DBConn(self.species)

        sizes = CHROM_SIZE[self.assembly]

        # compose a CASE statement to cap the upper bound of the QTLs by the size of the chromosome
        max_chrom = ' '.join(["WHEN '{}' THEN {}".format(chrom, sizes[chrom]) for chrom in sizes])

        # get all the non-neutral regions (defined here as all valid QTLs and gene regions +/- offset)
        results = dbc.get_records_sql("""

            # get all the QTLs with a valid rsnumber
            SELECT ev.chrom,
                   GREATEST(ev.start - {offset}, 1) AS `start`,
                   LEAST(ev.end + {offset}, CASE ev.chrom {max_chrom} END) AS `end`
              FROM qtls q
              JOIN ensembl_variants ev
                ON ev.rsnumber = q.peak
             WHERE q.associationType NOT IN ('Neutral', 'Sweep', 'MC1R')  # TODO improve clarity

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

        allregions = 'bed/{}_allregions.bed'.format(self.species)
        nonneutral = 'bed/{}_nonneutral.bed'.format(self.species)

        # write a BED file for the whole genome
        with open(allregions, 'w') as fout:
            for chrom in self.chromosomes:
                fout.write('{}\t{}\t{}\n'.format(chrom, 1, CHROM_SIZE[self.assembly][chrom]))

        # write all the non-neutral regions to a BED file
        with open(nonneutral, 'w') as fout:
            for chrom in natsorted(intervals.keys()):
                # merge overlapping intervals
                for start, stop in merge_intervals(intervals[chrom], capped=False):
                    fout.write('{}\t{}\t{}\n'.format(chrom, start, stop))

        # subtract the non-neutral regions from the whole genome
        loci = run_cmd(['bedtools', 'subtract', '-a', allregions, '-b', nonneutral])

        num_loci = 0

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

            # check if this QTL already exists
            if dbc.get_record('qtls', qtl):
                continue

            dbc.save_record('qtls', qtl)

            num_loci += 1

        with self.output().open('w') as fout:
            fout.write('INFO: Added {:,} neutral loci'.format(num_loci))


class PopulateAllLoci(PipelineWrapperTask):
    """
    Wrapper task to populate all the QTL and pseudo-QTL windows.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        # load all the QTLs and
        yield PopulateTraitLoci(self.species)
        yield PopulateNeutralLoci(self.species)


class PopulateQTLSNPs(PipelineTask):
    """
    Now we have ascertained all the modern SNPs, let's find those that intersect with the QTLs.

    :type species: str
    :type population: str
    :type chrom: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        yield PopulateAllLoci(self.species)

    def output(self):
        return luigi.LocalTarget('db/{}-qtl_snps.log'.format(self.basename))

    def run(self):
        dbc = DBConn(self.species)

        # TODO make sure this work with DOM and DOM2
        # insert linking records to make future queries much quicker
        exec_time = dbc.execute_sql("""
            INSERT INTO qtl_snps (qtl_id, modsnp_id)
                 SELECT q.id, ms.id
                   FROM qtls q
                   JOIN modern_snps ms
                     ON ms.population = '{population}'
                    AND ms.chrom = q.chrom
                    AND ms.site BETWEEN q.start AND q.end
                  WHERE q.chrom = '{chrom}'
                    AND q.valid = 1
                    AND ms.daf >= {daf}
                    """.format(population=self.population, chrom=self.chrom, daf=MIN_DAF))

        with self.output().open('w') as fout:
            fout.write('INFO: Execution took {}'.format(exec_time))


class MarkNeutralSNPs(PipelineTask):
    """
    Mark neutral SNPs (i.e. SNPs outside of all QTLs and gene regions)

    :type species: str
    :type population: str
    :type chrom: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        return PopulateQTLSNPs(self.species, self.population, self.chrom)

    def output(self):
        return luigi.LocalTarget('db/{}-neutral_snps.log'.format(self.basename))

    def run(self):
        dbc = DBConn(self.species)

        exec_time = dbc.execute_sql("""
            UPDATE modern_snps ms
              JOIN qtl_snps qs
                ON qs.modsnp_id = ms.id
              JOIN qtls q
                ON q.id = qs.qtl_id
               SET ms.neutral = 1
             WHERE q.chrom = '{chrom}'
               AND q.associationType = 'Neutral'
               AND q.valid = 1
               """.format(chrom=self.chrom))  # TODO add self.population

        with self.output().open('w') as fout:
            fout.write('INFO: Execution took {}'.format(exec_time))


class QTLPipeline(PipelineWrapperTask):
    """
    Call SNPs using the bcftools `mpileup | call` workflow.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):

        # process all the populations in chromosome chunks
        for pop in self.populations:
            for chrom in self.chromosomes:
                # flag the modern SNPs which fall into 'neutral' regions
                yield MarkNeutralSNPs(self.species, pop, chrom)


if __name__ == '__main__':
    luigi.run()
