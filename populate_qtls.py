#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from db_conn import *
from qtldb_api import *
from qtldb_lib import *

from utilities import *

import gzip
import shlex

# offset to use for the QTL window (+/- 50 Kb)
QTL_WINDOW = 50000

# sizes of each chrom in the given assemblies
CHROM_SIZE = {

    # UMD_3.1.1
    'cattle': {'1': 158337067, '2': 137060424, '3': 121430405, '4': 120829699, '5': 121191424, '6': 119458736,
               '7': 112638659, '8': 113384836, '9': 105708250, '10': 104305016, '11': 107310763, '12': 91163125,
               '13': 84240350, '14': 84648390, '15': 85296676, '16': 81724687, '17': 75158596, '18': 66004023,
               '19': 64057457, '20': 72042655, '21': 71599096, '22': 61435874, '23': 52530062, '24': 62714930,
               '25': 42904170, '26': 51681464, '27': 45407902, '28': 46312546, '29': 51505224, 'X': 148823899},

    # Sscrofa11.1
    'pig':    {'1': 274330532, '2': 151935994, '3': 132848913, '4': 130910915, '5': 104526007, '6': 170843587,
               '7': 121844099, '8': 138966237, '9': 139512083, '10': 69359453, '11': 79169978, '12': 61602749,
               '13': 208334590, '14': 141755446, '15': 140412725, '16': 79944280, '17': 63494081, '18': 55982971,
               'X': 125939595, 'Y': 43547828},

    # EquCab2.0
    'horse':  {'1': 185838109, '2': 120857687, '3': 119479920, '4': 108569075, '5': 99680356, '6': 84719076,
               '7': 98542428, '8': 94057673, '9': 83561422, '10': 83980604, '11': 61308211, '12': 33091231,
               '13': 42578167, '14': 93904894, '15': 91571448, '16': 87365405, '17': 80757907, '18': 82527541,
               '19': 59975221, '20': 64166202, '21': 57723302, '22': 49946797, '23': 55726280, '24': 46749900,
               '25': 39536964, '26': 41866177, '27': 39960074, '28': 46177339, '29': 33672925, '30': 30062385,
               '31': 24984650, 'X': 124114077},

    # CHIR_1.0
    'goat':   {'1': 155011307, '2': 135415751, '3': 116796116, '4': 115961478, '5': 111055201, '6': 114334461,
               '7': 106547263, '8': 111020524, '9': 90293942, '10': 99198151, '11': 105305070, '12': 82535142,
               '13': 80625018, '14': 92306894, '15': 78986926, '16': 77678508, '17': 71877645, '18': 61067880,
               '19': 62130014, '20': 71279863, '21': 66773250, '22': 57956300, '23': 49403180, '24': 61756751,
               '25': 41496684, '26': 50169583, '27': 44118849, '28': 43231948, '29': 48376377, 'X': 121952644},
}

ENSEMBL_DATA = {

    # see ftp://ftp.ensembl.org/pub/release-89/variation/gvf/sus_scrofa/
    #     ftp://ftp.ensembl.org/pub/release-89/gtf/sus_scrofa/
    'pig': {'gtf': 'data/ensembl/release-89/gtf/sus_scrofa/Sus_scrofa.Sscrofa10.2.89.gtf.gz',
            'gvf': 'data/ensembl/release-89/gvf/sus_scrofa/Sus_scrofa.gvf.gz'},
}

SWEEP_DATA = {

    # see https://www.nature.com/articles/ng.3394
    'pig': {'loci': 'data/sweep/EUD_Sweep_p001_FINAL_cutoff_MERGED10kb.bed',
            'snps': 'data/sweep/EUD_Sweep_p001_FINAL_cutoff.bed'}
}

SWEEP_NUM_SNPS = 3

SNP_CHIP_DATA = {

    # see http://bioinformatics.tecnoparco.org/SNPchimp
    'pig': 'data/SNPchimp/SNPchimp_pig.tsv.gz'
}

# TODO make global variable
VERBOSE = True


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

            qtl_id = dbc.save_record('qtls',qtl)

            num_loci += 1

            # get the all the SNPs from this locus
            snps = run_cmd(["printf '{locus}' | bedtools intersect -a {snps_file} -b stdin"
                           .format(locus=locus.strip(),snps_file=snps_file)], shell=True)

            for snp in snps.splitlines():
                # extract the SNP info
                chrom, start, end, cdf, p = snp.split()

                sweep_snp = {
                    'qtl_id': qtl_id,
                    'chrom':  chrom,
                    'site':   start,
                    'cdf':    cdf,
                    'p':      p
                }

                dbc.save_record('sweep_snps', sweep_snp)

    print("INFO: Loaded {} selective sweep loci for {}.".format(num_loci, species))


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

    # get all the QTL windows
    results = dbc.get_records_sql("""
        SELECT q.id, v.chrom, v.start AS site
          FROM qtls q
          JOIN ensembl_variants v
            ON v.rsnumber = q.peak                # only QTLs with a dbsnp peak
           AND v.type = 'SNV'                     # only single nucleotide variants
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

        # update the QTL record
        qtl = {'id': result['id'], 'chrom': result['chrom'], 'valid': 1, 'site': result['site'], 'start': start, 'end': end}

        dbc.save_record('qtls', qtl)

