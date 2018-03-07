#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob, itertools, sys, socket

from Bio import SeqIO
from collections import Counter

from db_conn import *

# use the Pathos library for improved multi-processing
import pathos.multiprocessing as mp

if socket.gethostname() == 'macbookpro.local':
    # use test dataset
    PATH = '/Users/Evan/Dropbox/Code/alleletraj/fasta'
    CHROMS = ['10', '11']
    THREADS = 1

else:
    # use real dataset
    PATH = '/media/jbod/raid1-sdc1/laurent/full_run_results/Pig/modern/FASTA'
    CHROMS = [str(chrom) for chrom in range(1, 19)] + ['X', 'Y']
    THREADS = 20


SPECIES = 'pig'
OUTGROUP = 'SVSV01U01_Sverrucosus_rh'
EXCLUDE = 'Sverrucosus' # Sus verrucosus / Javan warty pig

dbc = db_conn()

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


def process_chrom(chrom):
    """
    Extract all the SNPs from a given chromosome and estiamte the allele frequency
    """
    site = 0
    num_snps = 0

    # get all the fasta files (excluding the outgroup)
    fasta_files = [fasta for fasta in glob.glob(PATH + "/*/{}.fa".format(chrom)) if EXCLUDE not in fasta]

    # get the outgroup fasta file
    outgroup_fasta = PATH + '/' + OUTGROUP + "/{}.fa".format(chrom)

    # add the outgroup to the front of the list
    fasta_files.insert(0, outgroup_fasta)

    print "START: Parsing chr{} from {} fasta files.".format(chrom, len(fasta_files))

    data = {}

    for fasta in fasta_files:
        # load the fasta data
        data[fasta] = SeqIO.index(fasta, "fasta")

        print "LOADED: {}".format(fasta)

    # extract the sequence data from the current chrom in each fasta file
    seqs = [data[fasta][chrom].seq for fasta in fasta_files]

    # zip the sequences together so we can iterate over them one site at a time
    for pileup in itertools.izip(*seqs):

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
                print >> sys.stderr, "WARNING: Unknown ancestral allele chr{}:{} = {}".format(chrom, site, outgroup_allele)
                continue

            ancestral = ancestral.pop()
            alleles = observations.keys()

            if ancestral not in alleles:
                print >> sys.stderr, "WARNING: Pollyallelic site chr{}:{} = {}, ancestral {}".format(chrom, site, pileup, ancestral)
                continue

            alleles.remove(ancestral)
            derived = alleles.pop()

            record = dict()
            record['species'] = SPECIES
            record['chrom'] = chrom
            record['site'] = site
            record['ancestral'] = ancestral
            record['ancestral_count'] = observations[ancestral]
            record['derived'] = derived
            record['derived_count'] = observations[derived]

            dbc.save_record('modern_snps', record)
            num_snps += 1

        elif num_alleles > 2:
            print >> sys.stderr, "WARNING: Pollyallelic site chr{}:{} = {}".format(chrom, site + 1, pileup)

    print "FINISHED: chr{} contained {} SNPs".format(chrom, num_snps)


if THREADS > 1:
    # process the chromosomes in parallel
    pool = mp.ProcessingPool(THREADS)
    results = pool.map(process_chrom, CHROMS)
else:
    for chrom in CHROMS:
        process_chrom(chrom)
