#!/usr/bin/env python
# -*- coding: utf-8 -*-

from db_conn import *
from qtldb_api import *
from qtldb_lib import *

# maximum length of the QTL to process (100 Kb)
MAX_QTL_LENGTH = 100000

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

    print "INFO: Processing %s %s QTLs from '%s'" % (len(data['QTL_ID']), species, qtldb_file)

    # convert all the IDs to int
    qtl_ids = [int(qtl_id) for qtl_id in data['QTL_ID'] if qtl_id.isdigit()]

    # get all the QTLs already in the DB
    qtls = dbc.get_records('qtls', key='id')

    # find the new IDs in the list
    new_ids = list(set(qtl_ids) - set(qtls.keys()))

    if VERBOSE:
        print "INFO: Found %s new %s QTLs to add" % (len(new_ids), species)

    # rename these fields
    key_map = {
        'pubmedID':   'pubmed_id',
        'traitID':    'trait_id',
        'geneID':     'gene_id',
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
            trait['type'] = api.get_trait_type(species, trait['id'], trait['name'])

            dbc.save_record('traits', trait)

        # does the publication exist
        if not dbc.exists_record('pubmeds', {'id': record['pubmed_id']}):

            # setup the pubmed record
            pubmed = api.get_publication(species, record['pubmed_id'])

            if pubmed:
                pubmed['id'] = pubmed.pop('pubmed_ID')
                pubmed['year'] = re.search('\(([0-9]{4})\)', pubmed['authors']).group(1)
                pubmed['journal'] = pubmed['journal']['#text'][:-5]

                dbc.save_record('pubmeds', pubmed)
            else:
                # some records have a bogus pubmed_id
                record['pubmed_id'] = None

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

        dbc.save_record('qtls', qtl)

        added += 1

        if VERBOSE and added % 100 == 0:
            print "INFO: Added %5d new QTLs" % added

    if VERBOSE:
        print "INFO: Finished adding %s new %s QTLs" % (len(new_ids), species)


def compute_qtl_windows(species):

    # open a db connection
    dbc = db_conn()

    # get all the QTL windows
    results = dbc.get_records_sql("""
        SELECT DISTINCT q.id, d.chrom, d.site
          FROM qtls q
          JOIN dbsnp d on d.rsnumber = q.peak
         WHERE q.species = '{species}' 
           AND q.peak like 'rs%'                  # only QTLs with a dbsnp peak
           AND q.associationType = 'Association'  # only GWAS studies
           AND q.significance = 'Significant'     # only significant hits
           AND IFNULL(q.genomeLoc_end - q.genomeLoc_start, 0) <= {max}
           """.format(species=species, max=MAX_QTL_LENGTH), key=None)

    print "INFO: Computing window sizes for {:,} {} QTLs".format(len(results), species)

    # set offset to be half the max QTL length
    offset = MAX_QTL_LENGTH / 2

    for result in results:
        # get the size of the current chrom
        chom_size = CHROM_SIZE[species][result['chrom']]

        # calculate the bounded window size
        start = result['site'] - offset if result['site'] > offset else 1
        end = result['site'] + offset if result['site'] + offset < chom_size else chom_size

        # update the QTL record
        qtl = {'id': result['id'], 'valid': 1, 'start': start, 'end': end}

        dbc.save_record('qtls', qtl)


def cache_dbsnp_rsnumbers(species):

    # open a db connection
    dbc = db_conn()

    # copy all the rsnumbers we need into a new table (saves on costly lookups in `dbsnp_complete` which exceeds RAM)
    dbc.execute_sql("""
        INSERT INTO dbsnp
        SELECT DISTINCT dc.*
          FROM dbsnp_complete dc
          JOIN qtls q
            ON q.peak = dc.rsnumber
         WHERE q.species = '{species}'
           AND q.valid = 1
           """.format(species))

