#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import itertools
import multiprocessing as mp
import sys
import traceback
from collections import Counter
import pysam as ps

from pipeline_utils import *

# TODO convert to luigi pipeline
# TODO update to use pipeline_snp_call
# NB. polarized  should we drop sites that are hom-REF in outgroup and hom-ALT in ingroup?


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


def process_fasta_files(args):
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
        fasta_files = [FASTA_PATH + '/' + sample + "/{}.fa".format(chrom) for sample in SAMPLES[SPECIES][population]]

        # get the outgroup fasta file
        outgroup_fasta = FASTA_PATH + '/' + OUTGROUP + "/{}.fa".format(chrom)

        # add the outgroup to the front of the list
        fasta_files.insert(0, outgroup_fasta)

        print("STARTED: Parsing {} chr{} from {} fasta files.".format(population, chrom, len(fasta_files)))

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
                        print("WARNING: Pollyallelic site chr{}:{} = {}, ancestral {}"
                              .format(chrom, site, set(haploids), ancestral), file=sys.stderr)
                    continue

                # is this mutation a transition (A <-> G and C <-> T) or a transversion (everything else)
                type = 'ts' if set(alleles) == {'A', 'G'} or set(alleles) == {'C', 'T'} else 'tv'

                alleles.remove(ancestral)
                derived = alleles.pop()

                # calculate the derived allele frequency (daf)
                daf = float(observations[derived]) / (observations[ancestral] + observations[derived])

                record = dict()
                record['population'] = population
                record['chrom'] = chrom
                record['site'] = site
                record['ancestral'] = ancestral
                record['ancestral_count'] = observations[ancestral]
                record['derived'] = derived
                record['derived_count'] = observations[derived]
                record['type'] = type
                record['daf'] = daf

                dbc.save_record('modern_snps', record)
                num_snps += 1

            elif num_alleles > 2:
                if VERBOSE:
                    print("WARNING: Polyallelic site chr{}:{} = {}".format(chrom, site + 1, set(haploids)),
                          file=sys.stderr)

        print("FINISHED: {} chr{} contained {:,} SNPs".format(population, chrom, num_snps))

    except Exception:
        # Put all exception text into an exception and raise that
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))


def process_vcf_files(args):
    """
    Extract all the SNPs from a given chromosome and estimate the allele frequency
    """

    try:
        # extract the nested tuple of arguments (an artifact of using izip to pass args to mp.Pool)
        (chrom, population) = args

        num_snps = 0

        # open a private connection to the database
        dbc = db_conn()

        vcf_file = "vcf/horse_{}_chr{}.vcf".format(population, chrom)

        print("STARTED: Parsing {} chr{} in {}.".format(population, chrom, vcf_file))

        # parse the VCF with pysam
        for rec in ps.VariantFile(vcf_file).fetch():

            # skip low quality sites
            if int(rec.qual) < MIN_GENO_QUAL:
                continue

            # skip polyallelic sites
            if len(rec.alleles) != 2:
                continue

            # get the outgroup genotype
            out_geno = set(rec.samples[OUTGROUP]['GT'])

            # skip heterozygous sites in the outgroup (as they can't be polarized)
            if len(out_geno) != 1 or out_geno == {None}:
                continue

            # resolve the ancestral and derived alleles from the outgroup genotype
            ancestral = rec.alleles[out_geno.pop()]
            derived = [base for base in rec.alleles if base != ancestral].pop()

            # lets collate all the haploid observations for the two alleles
            haploids = []

            # resolve the genotypes of all the samples
            for sample in SAMPLES[SPECIES][population]:

                # get the genotype call for this site
                geno = rec.samples[sample]['GT']

                # decode the GT notation into allele calls (e.g. 0/0, 0/1, 1/1)
                alleles = [rec.alleles[idx] for idx in geno if idx is not None]

                # add them to the list
                haploids += alleles

            # count the haploid observations
            observations = Counter(haploids)

            if len(observations.keys()) != 2:
                # skip non-variant sites
                continue

            # is this mutation a transition (A <-> G and C <-> T) or a transversion (everything else)
            type = 'ts' if set(rec.alleles) == {'A', 'G'} or set(rec.alleles) == {'C', 'T'} else 'tv'

            # calculate the derived allele frequency (daf)
            daf = float(observations[derived]) / (observations[ancestral] + observations[derived])

            record = dict()
            record['population'] = population
            record['chrom'] = chrom
            record['site'] = rec.pos
            record['ancestral'] = ancestral
            record['ancestral_count'] = observations[ancestral]
            record['derived'] = derived
            record['derived_count'] = observations[derived]
            record['type'] = type
            record['daf'] = daf

            dbc.save_record('modern_snps', record)
            num_snps += 1

        print("FINISHED: {} chr{} contained {:,} SNPs".format(population, chrom, num_snps))

    except Exception:

        # Put all exception text into an exception and raise that
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))


def link_ensembl_genes():
    """
    Link modern SNPs to their Ensembl genes
    """
    dbc = db_conn()

    start = time()

    print("INFO: Linking modern SNPs to their Ensembl genes...", end='')

    # chunk the queries by chrom (to avoid temp tables)
    chroms = CHROM_SIZE[SPECIES].keys()

    for chrom in chroms:
        dbc.execute_sql("""
            UPDATE modern_snps ms
              JOIN ensembl_genes eg
                ON eg.chrom = ms.chrom
               AND ms.site BETWEEN eg.start AND eg.end
               SET ms.gene_id = eg.id
             WHERE ms.chrom = '{chrom}'""".format(chrom=chrom))

    print("({}).".format(timedelta(seconds=time() - start)))


def link_ensembl_variants():
    """
    Link modern SNPs to their Ensembl dbsnp variants
    """
    dbc = db_conn()

    start = time()

    print("INFO: Linking modern SNPs to their Ensembl dbsnp variants...", end='')

    # chunk the queries by chrom (to avoid temp tables)
    chroms = CHROM_SIZE[SPECIES].keys()

    for chrom in chroms:
        dbc.execute_sql("""
            UPDATE modern_snps ms
              JOIN ensembl_variants v
                ON ms.chrom = v.chrom
               AND ms.site = v.start
               SET ms.variant_id = v.id
             WHERE ms.chrom = '{chrom}'
               and v.type = 'SNV'        
               AND CHAR_LENGTH(alt) = 1   
               AND v.ref IN (ms.derived, ms.ancestral)
               AND v.alt IN (ms.derived, ms.ancestral)""".format(chrom=chrom))

    print("({}).".format(timedelta(seconds=time() - start)))


def link_snpchip():
    """
    Link modern SNPs to their snpchip variants
    """
    dbc = db_conn()

    start = time()

    print("INFO: Linking modern SNPs to their snpchip variants...", end='')

    # chunk the queries by chrom (to avoid temp tables)
    chroms = CHROM_SIZE[SPECIES].keys()

    for chrom in chroms:
        dbc.execute_sql("""
            UPDATE modern_snps ms
              JOIN snpchip sc
                ON sc.chrom = ms.chrom
               AND sc.site = ms.site
               SET ms.snpchip_id = sc.id
             WHERE ms.chrom = '{chrom}'""".format(chrom=chrom))

    print("({}).".format(timedelta(seconds=time() - start)))


def discover_modern_snps():
    """
    Ascertain SNPs in modern whole genome data.
    """

    if SPECIES == 'pig':
        # ascertain modern pig SNPs from whole genome FASTA files
        func = process_fasta_files

    elif SPECIES == 'horse':
        # ascertain modern horse SNPs from VCF files
        func = process_vcf_files

    else:
        raise Exception('Not implemented yet for {}'.format(SPECIES))

    chroms = CHROM_SIZE[SPECIES].keys()

    if MULTI_THREADED:
        # chain together an iterator for the params
        params = itertools.chain.from_iterable(itertools.izip(chroms, itertools.repeat(pop)) for pop in SAMPLES[SPECIES])

        # process the chromosomes in parallel
        pool = mp.Pool(CPU_CORES_MAX)
        pool.map(func, params)

    else:
        for chrom in chroms:
            for pop in SAMPLES[SPECIES]:
                # ascertain the modern SNPs separately in each population
                func((chrom, pop))


def links_modern_snps():
    """
    Link modern SNPs to their dbsnp, gene and snpchip records
    """
    link_ensembl_genes()
    link_ensembl_variants()
    link_snpchip()
