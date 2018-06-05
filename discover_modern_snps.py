#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import sys, socket
from collections import Counter

from pipeline_utils import *

import multiprocessing as mp
import itertools
import traceback

SPECIES = 'pig'

OUTGROUP = 'SVSV01U01_Sverrucosus_rh'  # Sus verrucosus / Javan warty pig

# path to the folder containing the modern fasta files
FASTA_PATH = '/media/jbod/raid1-sdc1/laurent/full_run_results/Pig/modern/FASTA'

SAMPLES = {

    # # the 81 European domestic pigs
    # 'EUD': ['AS01F01_AnglerSattleschwein_rh', 'AS01F09_Angler_Sattelsw_rh', 'BB01M47_Bunte_Bentheimer_rh',
    #        'BK01F10_Berkshire_rh', 'BK01M20_Berkshire_rh', 'BS01F10_British_Saddle_rh', 'BS01F35_British_Saddleback_rh',
    #        'CA01F14_Calabrese_rh', 'CM01F17_Chato_Murciano_rh', 'CM01F18_Chato_Murciano_rh', 'CS01F02_Cinta_Senese_rh',
    #        'CT01F13_Cassertana_rh', 'CT01M12_Cassertana_rh', 'DU22M01_Duroc_rh', 'DU22M02_Duroc_rh', 'DU22M03_Duroc_rh',
    #        'DU23M01_Duroc_rh', 'DU23M02_Duroc_rh', 'DU23M03_Duroc_rh', 'DU23M04_Duroc_rh', 'GO01F04_Gl_Old_Spots_rh',
    #        'GO01F23_GloucesterOldSpot_rh', 'HA20U01_Hampshire_rh', 'HA20U02_Hampshire_rh', 'HA20U04_Hampshire_rh',
    #        'HA20U06_Hampshire_rh', 'LB01F49_Large_Black_rh', 'LE01F25_Leicoma_rh', 'LR21M03_rh', 'LR24F01_rh', 'LR24F08_rh',
    #        'LR24M17_rh', 'LR24M18_rh', 'LR24M19_rh', 'LR24M20_rh', 'LR24M21_rh', 'LR24M22_rh', 'LR30F02_rh', 'LR30F03_rh',
    #        'LR30F04_Landrace_rh', 'LS01F04_Linderodsvin_rh', 'LW22F01_rh', 'LW22F02_rh', 'LW22F03_rh', 'LW22F04_rh',
    #        'LW22F06_rh', 'LW22F07_rh', 'LW22F08_LargeWhite_rh', 'LW22F09_LargeWhite_rh', 'LW22M04_rh', 'LW36F01_rh',
    #        'LW36F02_rh', 'LW36F03_rh', 'LW36F04_rh', 'LW36F05_rh', 'LW36F06_rh', 'LW37M01_rh', 'LW38MF02_rh', 'LW39M01_rh',
    #        'LW39M02_rh', 'LW39M03_rh', 'LW39M04_rh', 'LW39M05_rh', 'LW39M07_rh', 'LW39M08_rh', 'MA01F18_Mangalica_rh',
    #        'MA01F20_Mangalica_rh', 'MW01F29_Middle_White_rh', 'MW01F33_Middle_White_rh', 'NI01U07_Negro_Iberico_rh',
    #        'NS01F05_Nera_Siciliana_rh', 'PI21F02_rh', 'PI21F06_rh', 'PI21F07_Pietrain_rh', 'PI21F08_Pietrain_rh',
    #        'PI21F09_Pietrain_rh', 'PI21M17_rh', 'PI21M20_rh', 'RE01F51_Retinto_rh', 'TA01F19_Tamworth_rh',
    #        'TA01M06_Tamworth_rh'],

    # the 22 Asian domestic pigs
    'ASD': ['JI01U08_Jinhua_rh', 'JI01U10_Jinhua_rh', 'JQ01U02_Jiangquhai_rh', 'JQ01U03_Jiangquahai_rh',
           'JQ01U08_Jiangquahai_rh', 'LSP01U16_LepingSpotted_rh', 'LSP01U18_LepingSpotted_rh', 'MS20M03_Meishan_rh',
           'MS20M05_Meishan_rh', 'MS20U10_Meishan_rh', 'MS20U11_Meishan_rh', 'MS20U13_Meishan_rh', 'MS21M01_Meishan_rh',
           'MS21M05_Meishan_rh', 'MS21M07_Meishan_rh', 'MS21M08_Meishan_rh', 'MS21M14_Meishan_rh',
           'WS01U03_WannanSpotted_rh', 'WS01U13_WannanSpotted_rh', 'XI01U03_rh', 'XI01U04_Xiang_rh', 'ZA01U02_Zang_rh'],

    # the 2 Sumatran scrofa (for ascertaining ancient alleles)
    'SUM': ['INDO22_Sumatra_rh', 'INDO33_Sumatra_rh'],
}

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


def process_chrom(args):
    """
    Extract all the SNPs from a given chromosome and estimate the allele frequency
    """

    try:

        # extract the nested tuple of arguments (an artifact of using izip to pass args to mp.Pool)
        (chrom, population) = args

        site = 0
        num_snps = 0

        # open a private connection to the database
        dbc = db_conn()

        # get all the modern fasta files
        fasta_files = [FASTA_PATH + '/' + sample + "/{}.fa".format(chrom) for sample in SAMPLES[population]]

        # get the outgroup fasta file
        outgroup_fasta = FASTA_PATH + '/' + OUTGROUP + "/{}.fa".format(chrom)

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

            print("LOADED: {}".format(fasta))

        # zip the sequences together so we can iterate over them one site at a time
        for pileup in itertools.izip_longest(*data, fillvalue='N'):

            # increment the counter
            site += 1

            # convert iterator into a list
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
                type = 'ts' if set(alleles) == set(['A', 'G']) or set(alleles) == set(['C', 'T']) else 'tv'

                alleles.remove(ancestral)
                derived = alleles.pop()

                # calculate the minor allele frequency (maf)
                maf = observations[derived] / (observations[ancestral] + observations[derived])

                record = dict()
                record['population'] = population
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

    except Exception:
        # Put all exception text into an exception and raise that
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))


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

    chroms = CHROM_SIZE[species].keys()

    if MULTI_THREADED:
        # chain together an iterator for the params
        params = itertools.chain.from_iterable(itertools.izip(chroms, itertools.repeat(pop)) for pop in SAMPLES)

        # process the chromosomes in parallel
        pool = mp.Pool(MAX_CPU_CORES)
        pool.map(process_chrom, params)

    else:
        for chrom in chroms:
            for pop in SAMPLES:
                # we only ascertain the modern SNPs in European domestic (EUD)
                process_chrom((chrom, pop))

    # link modern SNPs to their dbsnp, gene and snpchip records
    link_ensembl_variants()
    link_ensembl_genes()
    link_dbsnp_snpchip()
