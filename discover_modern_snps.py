#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import sys, socket
from collections import Counter

from pipeline_utils import *

import multiprocessing as mp

# show lots of debugging output
VERBOSE = False

if socket.gethostname() == 'macbookpro.local':
    # use test dataset
    PATH = '/Users/Evan/Dropbox/Code/alleletraj/tmpdata/fasta'
    CHROMS = ['10', '11']
    THREADS = 2

else:
    # use real dataset
    PATH = '/media/jbod/raid1-sdc1/laurent/full_run_results/Pig/modern/FASTA'
    # TODO replace with proper CHROM reference
    CHROMS = [str(c) for c in range(1, 19)] + ['X', 'Y']
    THREADS = 20

SPECIES = 'pig'

OUTGROUP = 'SVSV01U01_Sverrucosus_rh'  # Sus verrucosus / Javan warty pig

# the 81 European domestic pigs
MODERN_SAMPLES = ['AS01F01_AnglerSattleschwein_rh', 'AS01F09_Angler_Sattelsw_rh', 'BB01M47_Bunte_Bentheimer_rh',
                  'BK01F10_Berkshire_rh', 'BK01M20_Berkshire_rh', 'BS01F10_British_Saddle_rh',
                  'BS01F35_British_Saddleback_rh', 'CA01F14_Calabrese_rh', 'CM01F17_Chato_Murciano_rh',
                  'CM01F18_Chato_Murciano_rh', 'CS01F02_Cinta_Senese_rh', 'CT01F13_Cassertana_rh',
                  'CT01M12_Cassertana_rh', 'DU22M01_Duroc_rh', 'DU22M02_Duroc_rh', 'DU22M03_Duroc_rh',
                  'DU23M01_Duroc_rh', 'DU23M02_Duroc_rh', 'DU23M03_Duroc_rh', 'DU23M04_Duroc_rh',
                  'GO01F04_Gl_Old_Spots_rh', 'GO01F23_GloucesterOldSpot_rh', 'HA20U01_Hampshire_rh',
                  'HA20U02_Hampshire_rh', 'HA20U04_Hampshire_rh', 'HA20U06_Hampshire_rh', 'LB01F49_Large_Black_rh',
                  'LE01F25_Leicoma_rh', 'LR21M03_rh', 'LR24F01_rh', 'LR24F08_rh', 'LR24M17_rh', 'LR24M18_rh',
                  'LR24M19_rh', 'LR24M20_rh', 'LR24M21_rh', 'LR24M22_rh', 'LR30F02_rh', 'LR30F03_rh',
                  'LR30F04_Landrace_rh', 'LS01F04_Linderodsvin_rh', 'LW22F01_rh', 'LW22F02_rh', 'LW22F03_rh',
                  'LW22F04_rh', 'LW22F06_rh', 'LW22F07_rh', 'LW22F08_LargeWhite_rh', 'LW22F09_LargeWhite_rh',
                  'LW22M04_rh', 'LW36F01_rh', 'LW36F02_rh', 'LW36F03_rh', 'LW36F04_rh', 'LW36F05_rh', 'LW36F06_rh',
                  'LW37M01_rh', 'LW38MF02_rh', 'LW39M01_rh', 'LW39M02_rh', 'LW39M03_rh', 'LW39M04_rh', 'LW39M05_rh',
                  'LW39M07_rh', 'LW39M08_rh', 'MA01F18_Mangalica_rh', 'MA01F20_Mangalica_rh', 'MW01F29_Middle_White_rh',
                  'MW01F33_Middle_White_rh', 'NI01U07_Negro_Iberico_rh', 'NS01F05_Nera_Siciliana_rh', 'PI21F02_rh',
                  'PI21F06_rh', 'PI21F07_Pietrain_rh', 'PI21F08_Pietrain_rh', 'PI21F09_Pietrain_rh', 'PI21M17_rh',
                  'PI21M20_rh', 'RE01F51_Retinto_rh', 'TA01F19_Tamworth_rh', 'TA01M06_Tamworth_rh']


# possible values include:
fasta_map = {
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'S': ['C', 'G'],
    'W': ['A', 'T'],
}


def decode_fasta(pileup):
    """
    Return list of haploid genotype calls, filtering out N values
    """
    return sum([fasta_map.get(val, [val, val]) for val in pileup if val != 'N'], [])


def stream_fasta(fin):
    """
    Convert a file input into an iterable which returns a single character, filtering out newline characters
    """
    for character in itertools.chain.from_iterable(fin):
        if character != "\n":
            yield character


def process_chrom(chrom):
    """
    Extract all the SNPs from a given chromosome and estimate the allele frequency
    """
    site = 0
    num_snps = 0

    # open a private connection to the database
    dbc = db_conn()

    # get all the modern fasta files
    fasta_files = [PATH + '/' + sample + "/{}.fa".format(chrom) for sample in MODERN_SAMPLES]

    # get the outgroup fasta file
    outgroup_fasta = PATH + '/' + OUTGROUP + "/{}.fa".format(chrom)

    # add the outgroup to the front of the list
    fasta_files.insert(0, outgroup_fasta)

    print("STARTED: Parsing chr{} from {} fasta files.".format(chrom, len(fasta_files)))

    data = []

    for fasta in fasta_files:
        # load the fasta data
        fin = open(fasta, "rU")

        # discard header row
        fin.readline()

        # wrap the file in an iterator
        data.append(stream_fasta(fin))

        if VERBOSE:
            print("LOADED: {}".format(fasta))

    # zip the sequences together so we can iterate over them one site at a time
    for pileup in itertools.izip_longest(*data, fillvalue='N'):

        # increment the counter
        site += 1

        pileup = list(pileup)

        # get the outgroup allele
        outgroup_allele = pileup.pop(0)

        # decode the fasta format
        haploids = decode_fasta(pileup)

        # how many alleles are there at this site
        num_alleles = len(set(haploids))

        # we're looking for biallelic SNPs
        if num_alleles == 2:

            # count the haploid observations
            observations = Counter(haploids)

            # decode the outgroup allele
            ancestral = set(decode_fasta(outgroup_allele))

            if len(ancestral) != 1:
                if VERBOSE:
                    print("WARNING: Unknown ancestral allele chr{}:{} = {}".format(chrom, site, outgroup_allele),
                          file=sys.stderr)
                continue

            ancestral = ancestral.pop()
            alleles = observations.keys()

            if ancestral not in alleles:
                if VERBOSE:
                    print("WARNING: Pollyallelic site chr{}:{} = {}, ancestral {}".format(chrom, site, set(haploids),
                                                                                          ancestral), file=sys.stderr)
                continue

            # is this mutation a transition (A <-> G and C <-> T) or a transversion (everything else)
            type = 'ts' if set(alleles) == set('A', 'G') or set(alleles) == set('C', 'T') else 'tv'

            alleles.remove(ancestral)
            derived = alleles.pop()

            # calculate the minor allele frequency (maf)
            maf = observations[derived] / (observations[ancestral] + observations[derived])

            record = dict()
            record['species'] = SPECIES
            record['chrom'] = chrom
            record['site'] = site
            record['ancestral'] = ancestral
            record['ancestral_count'] = observations[ancestral]
            record['derived'] = derived
            record['derived_count'] = observations[derived]
            record['type'] = type
            record['maf'] = maf

            dbc.save_record('modern_snps', record)
            num_snps += 1

        elif num_alleles > 2:
            if VERBOSE:
                print("WARNING: Pollyallelic site chr{}:{} = {}".format(chrom, site + 1, set(haploids)),
                      file=sys.stderr)

    print("FINISHED: chr{} contained {} SNPs".format(chrom, num_snps))


def link_ensembl_variants():
    """
    Link modern SNPs to their Ensembl dbsnp variants
    """
    dbc = db_conn()

    print("INFO: Linking modern SNPs to their Ensembl dbsnp variants.")

    dbc.execute_sql("""
        UPDATE modern_snps ms
          JOIN ensembl_variants v
            ON ms.chrom = v.chrom
           AND ms.site = v.start
           SET ms.variant_id = v.id
         WHERE v.type = 'SNV'        
           AND CHAR_LENGTH(alt) = 1   
           AND v.ref IN (ms.derived, ms.ancestral)
           AND v.alt IN (ms.derived, ms.ancestral)""")


def link_ensembl_genes():
    """
    Link modern SNPs to their Ensembl genes
    """
    dbc = db_conn()

    print("INFO: Linking modern SNPs to their Ensembl genes.")

    dbc.execute_sql("""
        UPDATE modern_snps ms
          JOIN ensembl_genes eg
            ON eg.chrom = ms.chrom
           AND ms.site BETWEEN eg.start AND eg.end
           SET ms.gene_id = eg.id""")


def link_dbsnp_snpchip():
    """
    Link modern SNPs to their snpchip variants
    """
    dbc = db_conn()

    print("INFO: Linking modern SNPs to their snpchip variants.")

    dbc.execute_sql("""
        UPDATE modern_snps ms
          JOIN dbsnp_snpchip ds
            ON ds.chrom = ms.chrom
           AND ds.site = ms.site
           SET ms.snpchip_id = ds.id""")


def discover_modern_snps(species):

    if species != 'pig':
        # TODO make this work for all species not just pigs
        raise Exception('Not implemented yet for {}'.format(species))

    if THREADS > 1:
        # process the chromosomes in parallel
        pool = mp.Pool(THREADS)
        pool.map(process_chrom, CHROMS)
    else:
        for chrom in CHROMS:
            process_chrom(chrom)

    # link modern SNPs to their dbsnp, gene and snpchip records
    link_ensembl_variants()
    link_ensembl_genes()
    link_dbsnp_snpchip()
