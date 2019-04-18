#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

from natsort import natsorted
import gzip

from qtldb_api import *
from pipeline_utils import *



class PopulateQTLs(PipelineTask):
    """
    Fetch all the QTLs from the QTLdb API and populate the local database.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return BlahBlah(self.species)

    def output(self):
        return luigi.LocalTarget('db/{}-qtls.log'.format(self.species))

    def run(self):
        dbc = db_conn()
        api = qtldb_api()

        # get the file containing all the QTL IDs
        qtldb_file = QTL_FILES[SPECIES]

        # get a list of all the QTLDb IDs from the tsv dump
        data = extract_qtl_fields(qtldb_file, ['QTL_ID'])

        # convert all the IDs to int
        qtl_ids = [int(qtl_id) for qtl_id in data['QTL_ID'] if qtl_id.isdigit()]

        print("INFO: Processing {:,} QTLs from '{}'".format(len(qtl_ids), qtldb_file))

        # get all the QTLs already in the DB
        qtls = dbc.get_records('qtls')

        # find the new IDs in the list
        new_ids = list(set(qtl_ids) - set(qtls.keys()))

        print("INFO: Found {:,} new QTLs to add".format(len(new_ids)))

        # rename these fields
        key_map = {
            'pubmedID':   'pubmed_id',
            'geneId':     'gene_id',
            'chromosome': 'chrom'
        }

        added = 0

        # get all the new records
        for record in api.get_qtls(SPECIES, new_ids):

            # TODO when resultset is len() = 1 then this throws an error
            # extract the nested trait record
            trait = record.pop('trait')
            trait['name'] = record.pop('name')

            # set the tait foreign key on the main record
            record['trait_id'] = trait['traitID']

            # does the trait exist
            if not dbc.exists_record('traits', {'id': record['trait_id']}):

                # setup the trait record
                trait = dict((field.replace('trait', '').lower(), value) for field, value in trait.iteritems())
                trait['type'] = ''  # api.get_trait_type(SPECIES, trait['id'], trait['name'])  # TODO this is broken!!

                dbc.save_record('traits', trait, insert=True)

            # does the publication exist
            if not dbc.exists_record('pubmeds', {'id': record['pubmedID']}):

                # setup the pubmed record
                pubmed = False # api.get_publication(SPECIES, record['pubmedID'])  # TODO this is broken!!

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
            record.pop('effects', None)
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

            if VERBOSE and added % 100 == 0:
                print("INFO: Added {:5d} new QTLs".format(added))

        # calculate the QTL windows
        ComputeQTLWindows()

        print("INFO: Finished adding {} new QTLs".format(len(new_ids)))


class PopulateSweeps(PipelineTask):
    """
    Populate the db with any selective sweep regions.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return BlahBlah(self.species)

    def output(self):
        return luigi.LocalTarget('db/{}-qtls.log'.format(self.species))

    def run(self):

        # open a db connection
        dbc = db_conn()

        try:
            # get the files containing the sweep data
            loci_file = SWEEP_DATA[SPECIES]['loci']
            snps_file = SWEEP_DATA[SPECIES]['snps']

        except KeyError:

            print("WARNING: No selective sweep loci for {}".format(SPECIES))
            return

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

                # check if this QTL already exists
                if dbc.get_record('qtls', qtl):
                    continue

                qtl_id = dbc.save_record('qtls', qtl)

                num_loci += 1

                # get the all the SNPs from this locus
                snps = run_cmd(["printf '{locus}' | bedtools intersect -a {snps_file} -b stdin"
                               .format(locus=locus.strip(), snps_file=snps_file)], shell=True)

                if not snps:
                    raise Exception("ERROR: Found no SNPs for sweep region {}:{}-{}".format(chrom, start, end))

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

        print("INFO: Loaded {} selective sweep loci (inc. {} SNPs)".format(num_loci, num_snps))


class PopulateMC1RLocus(PipelineTask):
    """
    Populate a dummy QTLs for the MC1R gene.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return BlahBlah(self.species)

    def output(self):
        return luigi.LocalTarget('db/{}-qtls.log'.format(self.species))

    def run(self):

        # TODO we want to model ALL 18 dbsnp variants, not just the 8 which are variable in EUD

        dbc = db_conn()

        # get the MC1R gene details
        mc1r = dbc.get_record('ensembl_genes', {'gene_id': MC1R_GENE_ID[SPECIES]})

        # setup a dummy QTL record
        qtl = {
            'associationType': 'MC1R',
            'chrom': mc1r['chrom'],
            'valid': 1,
            'start': mc1r['start'],
            'end': mc1r['end'],
        }

        # check if this QTL already exists
        if not dbc.get_record('qtls', qtl):
            dbc.save_record('qtls', qtl)

            print("INFO: Added the MC1R gene locus")


class PopulateNeutralLoci(PipelineTask):
    """
    Populate dummy QTLs for all the "neutral" loci (i.e. regions outside of all QTLs and gene regions +/- a buffer)

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return BlahBlah(self.species)

    def output(self):
        return luigi.LocalTarget('db/{}-qtls.log'.format(self.species))

    def run(self):
        dbc = db_conn()

        # compose a CASE statement to cap the upper bound of the QTLs by the size of the chromosome
        max_chrom = " ".join(["WHEN '{}' THEN {}".format(chrom, size) for chrom, size in CHROM_SIZE[SPECIES].iteritems()])

        # get all the non-neutral regions (defined here as all valid QTLs and gene regions +/- offset)
        results = dbc.get_records_sql("""
            
            # get all the QTLs with a valid rsnumber
            SELECT ev.chrom, 
                   GREATEST(ev.start - {offset}, 1) AS `start`,
                   LEAST(ev.end + {offset}, CASE ev.chrom {max_chrom} END) AS `end`
              FROM qtls q
              JOIN ensembl_variants ev
                ON ev.rsnumber = q.peak
             WHERE q.associationType NOT IN ('Neutral', 'Sweep', 'MC1R')
               
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

        allregions = 'bed/{}_allregions.bed'.format(SPECIES)
        nonneutral = 'bed/{}_nonneutral.bed'.format(SPECIES)

        # write a BED file for the whole genome
        with open(allregions, 'w') as fout:
            for chrom in CHROM_SIZE[SPECIES].keys():
                fout.write("{}\t{}\t{}\n".format(chrom, 1, CHROM_SIZE[SPECIES][chrom]))

        # write all the non-neutral regions to a BED file
        with open(nonneutral, 'w') as fout:
            for chrom in natsorted(intervals.keys()):
                # merge overlapping intervals
                for start, stop in merge_intervals(intervals[chrom], capped=False):
                    fout.write("{}\t{}\t{}\n".format(chrom, start, stop))

        # subtract the non-neutral regions from the whole genome
        loci = run_cmd(["bedtools", "subtract", "-a", allregions, "-b", nonneutral])

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

        print("INFO: Added {:,} neutral loci".format(num_loci))


class PopulatePigMummyLoci(PipelineTask):
    """
    Balancing selection on a recessive lethal deletion with pleiotropic effects on two neighboring genes in the porcine
    genome.

    See https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007661

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return BlahBlah(self.species)

    def output(self):
        return luigi.LocalTarget('db/{}-qtls.log'.format(self.species))

    def run(self):

        dbc = db_conn()

        # compose a CASE statement to cap the upper bound of the QTLs by the size of the chromosome
        max_chrom = " ".join(["WHEN '{}' THEN {}".format(chrom, size) for chrom, size in CHROM_SIZE[SPECIES].iteritems()])

        # get all pig mummy SNPs and turn them into pseudo-loci
        results = dbc.get_records_sql("""
            SELECT ev.rsnumber, ev.chrom, 
                   GREATEST(ev.start - {offset}, 1) AS `start`,
                   LEAST(ev.end + {offset}, CASE ev.chrom {max_chrom} END) AS `end`
              FROM ensembl_variants ev
             WHERE rsnumber IN ('rs321688936','rs324480289','rs331589427','rs333966168','rs334806226','rs337374395',
                                'rs340056046','rs340481593','rs341362223','rs342358904','rs81238716','rs81255938',
                                'rs81270496','rs81317835','rs81469247','rs81469256','rs81469273','rs81469281','rs81469291',
                                'rs81469298','rs81469311','rs81469316','rs81469337','rs81469339','rs81469341','rs81469348')
               """.format(offset=GENE_OFFSET, max_chrom=max_chrom), key=None)

        num_loci = 0

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

            # check if this QTL already exists
            if dbc.get_record('qtls', qtl):
                continue

            dbc.save_record('qtls', qtl)

            num_loci += 1

        print("INFO: Added {:,} pig mummy loci".format(num_loci))


class PopulateQTLSNPs(population):
    """
    Now we have ascertained all the modern SNPs, let's find those that intersect with the QTLs.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return BlahBlah(self.species)

    def output(self):
        return luigi.LocalTarget('db/{}-qtls.log'.format(self.species))

    def run(self):

        dbc = db_conn()

        start = time()

        print("INFO: Populating QTL SNPs... ", end='')

        chroms = CHROM_SIZE[SPECIES].keys()

        for chrom in chroms:
            # insert linking records to make future queries much quicker
            dbc.execute_sql("""
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
                        """.format(population=population, chrom=chrom, daf=MIN_DAF))

        print("({}).".format(timedelta(seconds=time() - start)))


class MarkNeutralSNPs(PipelineTask):
    """
    Mark neutral SNPs (i.e. SNPs outside of all QTLs and gene regions)

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return BlahBlah(self.species)

    def output(self):
        return luigi.LocalTarget('db/{}-qtls.log'.format(self.species))

    def run(self):

        dbc = db_conn()

        start = time()

        print("INFO: Marking neutral SNPs... ", end='')

        chroms = CHROM_SIZE[SPECIES].keys()

        for chrom in chroms:
            dbc.execute_sql("""
                UPDATE modern_snps ms
                  JOIN qtl_snps qs
                    ON qs.modsnp_id = ms.id
                  JOIN qtls q
                    ON q.id = qs.qtl_id
                   SET ms.neutral = 1
                 WHERE q.chrom = '{chrom}'
                   AND q.associationType = 'Neutral'
                   AND q.valid = 1
                   """.format(chrom=chrom))

        print("({}).".format(timedelta(seconds=time() - start)))


class ComputeQTLWindows(PipelineTask):
    """

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return BlahBlah(self.species)

    def output(self):
        return luigi.LocalTarget('db/{}-qtls.log'.format(self.species))

    def run(self):
        # open a db connection
        dbc = db_conn()

        # reset the existing QTL flags
        dbc.execute_sql("""
          UPDATE qtls
             SET valid = NULL,
                  site = NULL,
                 start = NULL,
                   end = NULL""")

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

        print("INFO: Computing window sizes for {:,} QTLs".format(len(results)))

        for result in results:

            # get the size of the current chrom
            chom_size = CHROM_SIZE[SPECIES][result['chrom']]

            # calculate the bounded window size
            start = result['site'] - QTL_WINDOW if result['site'] > QTL_WINDOW else 1
            end = result['site'] + QTL_WINDOW if result['site'] + QTL_WINDOW < chom_size else chom_size

            if end <= start:
                raise Exception("ERROR: Window size for QTL #{} is negative ({:,} bp)".format(result['id'], end-start))

            # update the QTL record
            qtl = {'id': result['id'], 'chrom': result['chrom'], 'valid': 1, 'site': result['site'], 'start': start,
                   'end': end}

            dbc.save_record('qtls', qtl)


class QTLPipeline(luigi.WrapperTask):
    """
    Call SNPs using the bcftools `mpileup | call` workflow.
    """
    def requires(self):

        # load the QTLs from the AnimalQTL database
        PopulateQTLs()

        # load pseudo-QTLs from other sources
        PopulateSweeps()
        PopulateMC1RLocus()

        if SPECIES == 'pig':
            PopulatePigMummyLoci()

        PopulateNeutralLoci()

        # TODO make this work with DOM and DOM2
        # link each QTL to the ascertained modern SNPs
        PopulateQTLSNPs(POPULATION)

        # flag the modern SNPs which fall into "neutral" regions
        MarkNeutralSNPs()


if __name__ == '__main__':
    luigi.run()
