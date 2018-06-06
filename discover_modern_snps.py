#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import sys, socket
from collections import Counter

from pipeline_utils import *

import multiprocessing as mp
import itertools
import traceback


def decode_fasta(pileup):
    """
    Return list of haploid genotype calls, filtering out N values
    """
    return sum([FASTA_MAP.get(val, [val, val]) for val in pileup if val != 'N'], [])


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
                type = 'ts' if set(alleles) == {'A', 'G'} or set(alleles) == {'C', 'T'} else 'tv'

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

    start = time()

    print("INFO: Linking modern SNPs to their Ensembl dbsnp variants...", end='')

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

    print("({}).".format(timedelta(seconds=time() - start)))


def link_ensembl_genes():
    """
    Link modern SNPs to their Ensembl genes
    """
    dbc = db_conn()

    start = time()

    print("INFO: Linking modern SNPs to their Ensembl genes...", end='')

    dbc.execute_sql("""
        UPDATE modern_snps ms
          JOIN ensembl_genes eg
            ON eg.chrom = ms.chrom
           AND ms.site BETWEEN eg.start AND eg.end
           SET ms.gene_id = eg.id""")

    print("({}).".format(timedelta(seconds=time() - start)))


def link_snpchip():
    """
    Link modern SNPs to their snpchip variants
    """
    dbc = db_conn()

    start = time()

    print("INFO: Linking modern SNPs to their snpchip variants...", end='')

    dbc.execute_sql("""
        UPDATE modern_snps ms
          JOIN snpchip sc
            ON sc.chrom = ms.chrom
           AND sc.site = ms.site
           SET ms.snpchip_id = sc.id""")

    print("({}).".format(timedelta(seconds=time() - start)))


def discover_modern_snps():

    if SPECIES != 'pig':
        # TODO make this work for all species not just pigs
        raise Exception('Not implemented yet for {}'.format(SPECIES))

    chroms = CHROM_SIZE[SPECIES].keys()

    if MULTI_THREADED:
        # chain together an iterator for the params
        params = itertools.chain.from_iterable(itertools.izip(chroms, itertools.repeat(pop)) for pop in SAMPLES)

        # process the chromosomes in parallel
        pool = mp.Pool(MAX_CPU_CORES)
        pool.map(process_chrom, params)

    else:
        for chrom in chroms:
            for pop in SAMPLES:
                # ascertain the modern SNPs separetly in each population
                process_chrom((chrom, pop))

    # link modern SNPs to their dbsnp, gene and snpchip records
    link_ensembl_variants()
    link_ensembl_genes()
    link_snpchip()
