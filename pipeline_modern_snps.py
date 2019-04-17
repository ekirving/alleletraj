#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import itertools
import luigi
import pysam
import sys

# import my custom modules
from pipeline_consts import *
from pipeline_utils import PipelineTask, db_conn

from collections import Counter
from datetime import timedelta
from time import time


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


class ProcessFASTAs(PipelineTask):
    """
    Extract all the SNPs from a given chromosome and estimate the allele frequency
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    # def requires(self):
    #     return ExtractSNPsVCF(self.species, self.population)

    # def output(self):
    #     return luigi.LocalTarget('sfs/{}/dadi/{}.sfs'.format(self.basename, self.population))

    def run(self):

        site = 0
        num_snps = 0

        # open a private connection to the database
        dbc = db_conn(self.species)

        # get all the modern fasta files
        fasta_files = [FASTA_PATH + '/' + sample + "/{}.fa".format(self.chrom)
                       for sample in SAMPLES[self.species][self.population]]

        # get the outgroup fasta file
        outgroup_fasta = FASTA_PATH + '/' + OUT_GROUP[self.species] + "/{}.fa".format(self.chrom)

        # add the outgroup to the front of the list
        fasta_files.insert(0, outgroup_fasta)

        print("STARTED: Parsing {} chr{} from {} fasta files.".format(self.population, self.chrom, len(fasta_files)))

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
                        print("WARNING: Unknown ancestral allele chr{}:{} = {}"
                              .format(self.chrom, site, outgroup_allele), file=sys.stderr)
                    continue

                ancestral = ancestral.pop()
                alleles = observations.keys()

                if ancestral not in alleles:
                    if VERBOSE:
                        print("WARNING: Pollyallelic site chr{}:{} = {}, ancestral {}"
                              .format(self.chrom, site, set(haploids), ancestral), file=sys.stderr)
                    continue

                # is this mutation a transition (A <-> G and C <-> T) or a transversion (everything else)
                snp_type = 'ts' if set(alleles) == {'A', 'G'} or set(alleles) == {'C', 'T'} else 'tv'

                alleles.remove(ancestral)
                derived = alleles.pop()

                # calculate the derived allele frequency (daf)
                daf = float(observations[derived]) / (observations[ancestral] + observations[derived])

                record = dict()
                record['population'] = self.population
                record['chrom'] = self.chrom
                record['site'] = site
                record['ancestral'] = ancestral
                record['ancestral_count'] = observations[ancestral]
                record['derived'] = derived
                record['derived_count'] = observations[derived]
                record['type'] = snp_type
                record['daf'] = daf

                dbc.save_record('modern_snps', record)
                num_snps += 1

            elif num_alleles > 2:
                if VERBOSE:
                    print("WARNING: Polyallelic site chr{}:{} = {}".format(self.chrom, site + 1, set(haploids)),
                          file=sys.stderr)

        print("FINISHED: {} chr{} contained {:,} SNPs".format(self.population, self.chrom, num_snps))


class ProcessVCFs(PipelineTask):
    """
    Extract all the SNPs from a given chromosome and estimate the allele frequency
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    # def requires(self):
    #     return ExtractSNPsVCF(self.species, self.population)

    # def output(self):
    #     return luigi.LocalTarget('sfs/{}/dadi/{}.sfs'.format(self.basename, self.population))

    def run(self):
        num_snps = 0

        # open a private connection to the database
        dbc = db_conn(self.species)

        vcf_file = "vcf/horse_{}_chr{}.vcf".format(self.population, self.chrom)

        print("STARTED: Parsing {} chr{} in {}.".format(self.population, self.chrom, vcf_file))

        # parse the VCF with pysam
        for rec in pysam.VariantFile(vcf_file).fetch():

            # skip low quality sites
            if int(rec.qual) < MIN_GENO_QUAL:
                continue

            # skip polyallelic sites
            if len(rec.alleles) != 2:
                continue

            # get the outgroup genotype
            out_geno = set(rec.samples[OUT_GROUP[self.species]]['GT'])

            # skip heterozygous sites in the outgroup (as they can't be polarized)
            if len(out_geno) != 1 or out_geno == {None}:
                continue

            # resolve the ancestral and derived alleles from the outgroup genotype
            ancestral = rec.alleles[out_geno.pop()]
            derived = [base for base in rec.alleles if base != ancestral].pop()

            # lets collate all the haploid observations for the two alleles
            haploids = []

            # resolve the genotypes of all the samples
            for sample in SAMPLES[self.species][self.population]:

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
            snp_type = 'ts' if set(rec.alleles) == {'A', 'G'} or set(rec.alleles) == {'C', 'T'} else 'tv'

            # calculate the derived allele frequency (daf)
            daf = float(observations[derived]) / (observations[ancestral] + observations[derived])

            record = dict()
            record['population'] = self.population
            record['chrom'] = self.chrom
            record['site'] = rec.pos
            record['ancestral'] = ancestral
            record['ancestral_count'] = observations[ancestral]
            record['derived'] = derived
            record['derived_count'] = observations[derived]
            record['type'] = snp_type
            record['daf'] = daf

            dbc.save_record('modern_snps', record)
            num_snps += 1

        print("FINISHED: {} chr{} contained {:,} SNPs".format(self.population, self.chrom, num_snps))


class DiscoverSNPs(luigi.WrapperTask):
    """
    Ascertain SNPs in modern whole genome data.
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        if self.species == 'pig':
            # special case of FASTA data for pig ascertainment
            yield ProcessFASTAs(self.species, self.population, self.chrom)
        else:
            # parse the VCF files
            yield ProcessVCFs(self.species, self.population, self.chrom)


class LinkEnsemblGenes(PipelineTask):
    """
    Link modern SNPs to their Ensembl genes
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        return DiscoverSNPs(self.species, self.population, self.chrom)

    # def output(self):
    #     return luigi.LocalTarget('sfs/{}/dadi/{}.sfs'.format(self.basename, self.population))

    def run(self):
        dbc = db_conn(self.species)

        start = time()

        print("INFO: Linking modern SNPs to their Ensembl genes...", end='')

        dbc.execute_sql("""
            UPDATE modern_snps ms
              JOIN ensembl_genes eg
                ON eg.chrom = ms.chrom
               AND ms.site BETWEEN eg.start AND eg.end
               SET ms.gene_id = eg.id
             WHERE ms.population = '{pop}'
               AND ms.chrom = '{chrom}'""".format(pop=self.population, chrom=self.chrom))

        print("({}).".format(timedelta(seconds=time() - start)))


class LinkEnsemblVariants(PipelineTask):
    """
    Link modern SNPs to their Ensembl dbsnp variants
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        return DiscoverSNPs(self.species, self.population, self.chrom)

    # def output(self):
    #     return luigi.LocalTarget('sfs/{}/dadi/{}.sfs'.format(self.basename, self.population))

    def run(self):
        dbc = db_conn(self.species)

        start = time()

        print("INFO: Linking modern SNPs to their Ensembl dbsnp variants...", end='')

        dbc.execute_sql("""
            UPDATE modern_snps ms
              JOIN ensembl_variants v
                ON ms.chrom = v.chrom
               AND ms.site = v.start
               SET ms.variant_id = v.id
             WHERE ms.population = '{pop}'
               AND ms.chrom = '{chrom}'
               AND v.type = 'SNV'        
               AND CHAR_LENGTH(alt) = 1   
               AND v.ref IN (ms.derived, ms.ancestral)
               AND v.alt IN (ms.derived, ms.ancestral)""".format(pop=self.population, chrom=self.chrom))

        print("({}).".format(timedelta(seconds=time() - start)))


class LinkSNPChip(PipelineTask):
    """
    Link modern SNPs to their snpchip variants
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        return DiscoverSNPs(self.species, self.population, self.chrom)

    # def output(self):
    #     return luigi.LocalTarget('sfs/{}/dadi/{}.sfs'.format(self.basename, self.population))

    def run(self):
        dbc = db_conn(self.species)

        start = time()

        print("INFO: Linking modern SNPs to their snpchip variants...", end='')

        dbc.execute_sql("""
            UPDATE modern_snps ms
              JOIN snpchip sc
                ON sc.chrom = ms.chrom
               AND sc.site = ms.site
               SET ms.snpchip_id = sc.id
             WHERE ms.population = '{pop}'
               AND ms.chrom = '{chrom}'""".format(pop=self.population, chrom=self.chrom))

        print("({}).".format(timedelta(seconds=time() - start)))


class DiscoverModernSNPs(luigi.WrapperTask):
    """
    Populate the modern_snps table and link records to genes, dbsnp and snpchip records.
    """

    def requires(self):
        species = 'horse'

        for pop in SAMPLES[species]:
            for chrom in CHROM_SIZE[species]:
                yield LinkEnsemblGenes(species, pop, chrom)
                yield LinkEnsemblVariants(species, pop, chrom)
                yield LinkSNPChip(species, pop, chrom)
