#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob, itertools, sys

from Bio import SeqIO
from collections import Counter

from db_conn import *

PATH = '/media/jbod/raid1-sdc1/laurent/full_run_results/Pig/modern/FASTA'
# PATH = '/Users/Evan/Dropbox/Code/alleletraj/fasta'

SPECIES = 'pig'

OUTGROUP = 'SVSV01U01_Sverrucosus_rh.fa'

EXCLUDE = 'Sverrucosus' # Sus verrucosus / Javan warty pig

dbc = db_conn()

# get all the fasta files (excluding the outgroup)
fasta_files = [fasta for fasta in glob.glob(PATH + "/*.fa") if EXCLUDE not in fasta]

data = {}
chroms = set()

for fasta in fasta_files:
    # load the fasta data
    data[fasta] = SeqIO.index(fasta, "fasta")

    # extract the list of chromosomes from the fasta file
    chroms |= set(data[fasta].keys())

    print "Loaded: {}".format(fasta)

# load the outgroup data
outgroup = SeqIO.index(PATH + '/' + OUTGROUP, "fasta")

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
    Return list of haploid genotype calls
    """
    if 'N' in pileup:
        pileup.remove('N')

    return sum([fasta_map.get(val, [val, val]) for val in pileup], [])


for chrom in sorted(chroms):

    print "Processing chr{}".format(chrom)

    # initialise the position counter
    pos = 0

    # extract the current chrom from each fasta file
    seqs = [data[fasta][chrom].seq for fasta in fasta_files]

    # zip the sequences together so we can iterate over them one site at a time
    for pileup in itertools.izip(*seqs):

        # increment the counter
        pos += 1

        # decode the fasta format
        haploids = decode_fasta(list(pileup))

        # how many alleles are there at this site
        num_alleles = len(set(haploids))

        # we're looking for biallelic SNPs
        if num_alleles == 2:

            # count the haploid observations
            observations = Counter(haploids)

            # get the ancestral allele
            ancestral = set(decode_fasta(list(outgroup[chrom].seq[pos - 1])))

            if len(ancestral) != 1:
                print >> sys.stderr, "WARNING: Unknown ancestral allele chr{}:{} = {}".format(chrom, pos, outgroup[chrom].seq[pos - 1])
                continue

            ancestral = ancestral.pop()
            alleles = observations.keys()

            if ancestral not in alleles:
                print >> sys.stderr, "WARNING: Pollyallelic site chr{}:{} = {}, ancestral {}".format(chrom, pos, pileup, ancestral)
                continue

            alleles.remove(ancestral)
            derived = alleles.pop()

            record = dict()
            record['species'] = SPECIES
            record['chrom'] = chrom
            record['site'] = pos
            record['ancestral'] = ancestral
            record['ancestral_count'] = observations[ancestral]
            record['derived'] = derived
            record['derived_count'] = observations[derived]

            dbc.save_record('modern_snps', record)

            print '.',

        elif num_alleles > 2:
            print >> sys.stderr, "WARNING: Pollyallelic site chr{}:{} = {}".format(chrom, pos + 1, pileup)
