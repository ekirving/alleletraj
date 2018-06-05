#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from db_conn import *
from qtldb_api import *
from qtldb_lib import *

from pipeline_utils import *

import gzip
import shlex


def populate_qtls(species):
    """
    Fetch all the QTLs from the QTLdb API and populate the local database.
    """
    dbc = db_conn()
    api = qtldb_api()

    # get the file containing all the QTL IDs
    qtldb_file = QTL_FILES[species]

    # get a list of all the QTLDb IDs from the tsv dump
    data = extract_qtl_fields(qtldb_file, ['QTL_ID'])

    # convert all the IDs to int
    qtl_ids = [int(qtl_id) for qtl_id in data['QTL_ID'] if qtl_id.isdigit()]

    print("INFO: Processing %s %s QTLs from '%s'" % (len(qtl_ids), species, qtldb_file))

    # get all the QTLs already in the DB
    qtls = dbc.get_records('qtls')

    # find the new IDs in the list
    new_ids = list(set(qtl_ids) - set(qtls.keys()))

    if VERBOSE:
        print("INFO: Found %s new %s QTLs to add" % (len(new_ids), species))

    # rename these fields
    key_map = {
        'pubmedID':   'pubmed_id',
        'geneId':     'gene_id',
        'chromosome': 'chrom'
    }

    added = 0

    # get all the new records
    for record in api.get_qtls(species, new_ids):

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
            trait['type'] = api.get_trait_type(species, trait['id'], trait['name'])  # TODO this is broken!!

            dbc.save_record('traits', trait, insert=True)

        # does the publication exist
        if not dbc.exists_record('pubmeds', {'id': record['pubmedID']}):

            # setup the pubmed record
            pubmed = api.get_publication(species, record['pubmedID'])

            if pubmed:
                pubmed['id'] = pubmed.pop('pubmed_ID')
                pubmed['year'] = re.search('\(([0-9]{4})\)', pubmed['authors']).group(1)
                pubmed['journal'] = pubmed['journal']['#text'][:-5]

                dbc.save_record('pubmeds', pubmed, insert=True)
            else:
                # some records have a bogus pubmed ID
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

        qtl['species'] = species

        dbc.save_record('qtls', qtl, insert=True)

        added += 1

        if VERBOSE and added % 100 == 0:
            print("INFO: Added %5d new QTLs" % added)

    if VERBOSE:
        print("INFO: Finished adding %s new %s QTLs" % (len(new_ids), species))


def populate_sweeps(species):
    """
    Populate the db with selective sweep regions for the given species.
    """

    # open a db connection
    dbc = db_conn()

    try:
        # get the files containing the sweep data
        loci_file = SWEEP_DATA[species]['loci']
        snps_file = SWEEP_DATA[species]['snps']

    except KeyError:

        print("INFO: No selective sweep loci for {}".format(species))
        return

    num_loci = 0
    num_snps = 0

    with open(loci_file, 'r') as loci_fin:

        for locus in loci_fin:
            # extract the locus info from the BED file
            chrom, start, end = locus.split()

            # setup a dummy QTL record
            qtl = {
                'species':         species,
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
                           .format(locus=locus.strip(),snps_file=snps_file)], shell=True)

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

    print("INFO: Loaded {} selective sweep loci (inc. {} SNPs) for {}.".format(num_loci, num_snps, species))


def populate_mc1r_locus(species):
    """
    Populate a dummy QTLs for the MC1R gene.
    """

    dbc = db_conn()

    # get the MC1R gene details
    mc1r = dbc.get_record('ensembl_genes', {'gene_id': MC1R_GENE_ID})

    # setup a dummy QTL record
    qtl = {
        'species': species,
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


def populate_neutral_loci(species):
    """
    Populate dummy QTLs for all the "neutral" loci (i.e. regions outside of all QTLs and gene regions +/- a buffer)
    """

    dbc = db_conn()

    # get all the non-neutral regions (defined here as all valid QTLs and gene regions +/- offset)
    results = dbc.get_records_sql("""
        SELECT chrom, GREATEST(start - {offset}, 0) AS start, end + {offset} AS end 
          FROM qtls
         WHERE associationType != 'Neutral'
           AND valid = 1
           
         UNION
         
        SELECT chrom, GREATEST(start - {offset}, 0) AS start, end + {offset} AS end
          FROM ensembl_genes
          
      ORDER BY chrom, start, end
           """.format(offset=GENE_OFFSET), key=None)

    intervals = defaultdict(list)

    # group the intervals by chrom
    for result in results:
        intervals[result['chrom']].append((result['start'], result['end']))

    allregions = 'bed/{}_allregions.bed'.format(species)
    nonneutral = 'bed/{}_nonneutral.bed'.format(species)

    # write a BED file for the whole genome
    with open(allregions, 'w') as fout:
        for chrom in natsorted(CHROM_SIZE[species].keys()):
            fout.write("{}\t{}\t{}\n".format(chrom, 0, CHROM_SIZE[species][chrom]))

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
            'species': species,
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


def mark_neutral_snps():
    """
    Mark neutral SNPs (i.e. SNPs outside of all QTLs and gene regions)
    """

    dbc = db_conn()

    start = time()

    print("INFO: Marking neutral SNPs... ", end='')

    dbc.execute_sql("""
        UPDATE modern_snps ms
          JOIN qtl_snps qs
            ON qs.modsnp_id = ms.id
          JOIN qtls q
            ON q.id = qs.qtl_id
           SET ms.neutral = 1
         WHERE q.associationType = 'Neutral'
           AND q.valid = 1""")

    print("({}).".format(timedelta(seconds=time() - start)))



def load_ensembl_variants(species):
    """
    Load the GVF (Genome Variation Format) data from Ensembl.

    Contains all germline variations from the current Ensembl release for this species.

    See http://www.sequenceontology.org/gvf.html
    """

    # open a db connection
    dbc = db_conn()

    print("INFO: Loading Ensembl {} variants".format(species))

    # the column headers for batch inserting into the db
    fields = ('species', 'dbxref', 'rsnumber', 'type', 'chrom', 'start', 'end', 'ref', 'alt')

    # open the GVF file
    with gzip.open(ENSEMBL_DATA[species]['gvf'], 'r') as fin:

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
                variant = (species, dbxref, rsnumber, type, chrom, start, end, atts['Reference_seq'], atts['Variant_seq'])

                # add to the buffer
                records.append(variant)

                # bulk insert in chunks
                if len(records) == MAX_INSERT_SIZE:
                    dbc.save_records('ensembl_variants', fields, records)
                    records = []

        # insert any the remaing variants
        if records:
            dbc.save_records('ensembl_variants', fields, records)




def load_ensembl_genes(species):
    """
    Load the GTF (General Transfer Format) data from Ensembl.

    GTF provides access to all annotated transcripts which make up an Ensembl gene set.

    See https://www.ensembl.org/info/website/upload/gff.html#fields
    """

    # open a db connection
    dbc = db_conn()

    print("INFO: Loading Ensembl {} genes".format(species))

    # the column headers for batch inserting into the db
    fields = ('species', 'source', 'gene_id', 'version', 'biotype', 'chrom', 'start', 'end')

    # open the GTF file
    with gzip.open(ENSEMBL_DATA[species]['gtf'], 'r') as fin:

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
                    atts.update({k:v.strip('"') for k, v in atts.items()})

                    # setup the record to insert, in this order
                    gene = (species, atts['gene_source'], atts['gene_id'], atts['gene_version'], atts['gene_biotype'], chrom, start, end)

                    # add to the buffer
                    records.append(gene)

        # bulk insert all the reads for this sample
        if records:
            dbc.save_records('ensembl_genes', fields, records)


def load_snpchip_variants(species):
    """
    Load the SNP chip variants from SNPchimp

    See http://bioinformatics.tecnoparco.org/SNPchimp
    """

    # open a db connection
    dbc = db_conn()

    print("INFO: Loading SNP chip {} variants".format(species))

    pipe = "/tmp/SNPchimp_{}".format(species)

    # unzip the dump into a named pipe
    run_cmd(["mkfifo --mode=0666 {pipe}".format(pipe=pipe)], shell=True)
    run_cmd(["gzip --stdout -d  {gz} > {pipe}".format(gz=SNP_CHIP_DATA[species], pipe=pipe)],
            shell=True, background=True)

    dbc.execute_sql("""
        LOAD DATA 
     LOCAL INFILE '{pipe}'
       INTO TABLE dbsnp_snpchip 
           IGNORE 1 LINES (chip_name, rsnumber, chrom, site, snp_name)
              """.format(pipe=pipe))


def compute_qtl_windows(species):

    # open a db connection
    dbc = db_conn()

    # reset the existing QTL flags
    dbc.execute_sql("""
      UPDATE qtls
         SET valid = NULL,
              site = NULL,
             start = NULL,
               end = NULL
       WHERE q.species = '{species}'
         """.format(species=species))

    # get all the valid QTL windows
    results = dbc.get_records_sql("""
        SELECT q.id, v.chrom, v.start AS site
          FROM qtls q
          JOIN ensembl_variants v
            ON v.rsnumber = q.peak                # only QTLs with a dbsnp peak
           AND v.type = 'SNV'                     # only single nucleotide variants
           AND CHAR_LENGTH(v.alt) = 1             # only biallelic sites
         WHERE q.species = '{species}'             
           AND q.associationType = 'Association'  # only GWAS studies
           AND q.significance = 'Significant'     # only significant hits
           """.format(species=species), key=None)

    print("INFO: Computing window sizes for {:,} {} QTLs".format(len(results), species))

    for result in results:

        # get the size of the current chrom
        chom_size = CHROM_SIZE[species][result['chrom']]

        # calculate the bounded window size
        start = result['site'] - QTL_WINDOW if result['site'] > QTL_WINDOW else 1
        end = result['site'] + QTL_WINDOW if result['site'] + QTL_WINDOW < chom_size else chom_size

        if end <= start:
            raise Exception("ERROR: Window size for QTL #{} is negative ({:,} bp)".format(result['id'], end-start))

        # update the QTL record
        qtl = {'id': result['id'], 'chrom': result['chrom'], 'valid': 1, 'site': result['site'], 'start': start, 'end': end}

        dbc.save_record('qtls', qtl)

