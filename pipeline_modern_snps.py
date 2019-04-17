#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import itertools
import luigi
import pysam
import sys

from collections import Counter

# import my custom modules
from pipeline_consts import SAMPLES, OUT_GROUP, CHROM_SIZE
from pipeline_snp_call import BiallelicSNPsVCF
from pipeline_utils import PipelineTask, db_conn

# path to the folder containing the modern fasta files
FASTA_PATH = '/media/jbod/raid1-sdc1/laurent/full_run_results/Pig/modern/FASTA'

# extended FASTA codes for biallelic sites
FASTA_MAP = {
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

    # TODO
    # def requires(self):
    #     return ExtractSNPsVCF(self.species, self.population)

    # TODO
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
                    print("WARNING: Unknown ancestral allele chr{}:{} = {}"
                          .format(self.chrom, site, outgroup_allele), file=sys.stderr)
                    continue

                ancestral = ancestral.pop()
                alleles = observations.keys()

                if ancestral not in alleles:
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
                print("WARNING: Polyallelic site chr{}:{} = {}".format(self.chrom, site + 1, set(haploids)),
                      file=sys.stderr)

        print("FINISHED: {} chr{} contained {:,} SNPs".format(self.population, self.chrom, num_snps))


class AlleleFrequencyFromVCF(PipelineTask):
    """
    Estimate the allele frequency of all the SNPs in a given chromosome VCF.

    :type species: str
    :type population: str
    :type chrom: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        return BiallelicSNPsVCF(self.species, self.population)

    def output(self):
        return luigi.LocalTarget('db/{}-modern_snps.log'.format(self.basename))

    def run(self):
        dbc = db_conn(self.species)

        # count the number of SNPs added
        num_snps = 0

        # parse the VCF with pysam
        for rec in pysam.VariantFile(self.input().path).fetch():

            # the VCF has already been polarised, so the REF/ALT are the ancestral/derived
            ancestral = rec.ref
            derived = rec.alts[0]

            # lets collate all the haploid observations for the two alleles
            haploids = []

            # resolve the genotypes of all the samples
            for sample in SAMPLES[self.species][self.population]:
                # get the alleles, but skip any missing genotypes
                alleles = [a for a in rec.samples[sample].alleles if a is not None]

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

            record = {
                'population': self.population,
                'chrom': self.chrom,
                'site': rec.pos,
                'ancestral': ancestral,
                'ancestral_count': observations[ancestral],
                'derived': derived,
                'derived_count': observations[derived],
                'type': snp_type,
                'daf': daf,
            }

            dbc.save_record('modern_snps', record)
            num_snps += 1

        with self.output().open('w') as fout:
            fout.write("Added {:,} SNPs".format(num_snps))


class ProcessSNPs(luigi.WrapperTask):
    """
    Ascertain modern SNPs in whole-genome data.

    :type species: str
    :type population: str
    :type chrom: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        if self.species == 'pig':
            # special case of FASTA data for pig ascertainment
            yield ProcessFASTAs(self.species, self.population, self.chrom)
        else:
            # parse the VCF files for all other species
            yield AlleleFrequencyFromVCF(self.species, self.population, self.chrom)


class LinkEnsemblGenes(PipelineTask):
    """
    Link modern SNPs to their Ensembl genes

    :type species: str
    :type population: str
    :type chrom: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        return ProcessSNPs(self.species, self.population, self.chrom)

    def output(self):
        return luigi.LocalTarget('db/{}-ensembl_genes.log'.format(self.basename))

    def run(self):
        dbc = db_conn(self.species)

        exec_time = dbc.execute_sql("""
            UPDATE modern_snps ms
              JOIN ensembl_genes eg
                ON eg.chrom = ms.chrom
               AND ms.site BETWEEN eg.start AND eg.end
               SET ms.gene_id = eg.id
             WHERE ms.population = '{pop}'
               AND ms.chrom = '{chrom}'""".format(pop=self.population, chrom=self.chrom))

        with self.output().open('w') as fout:
            fout.write('Execution took {}'.format(exec_time))


class LinkEnsemblVariants(PipelineTask):
    """
    Link modern SNPs to their Ensembl dbsnp variants

    :type species: str
    :type population: str
    :type chrom: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        return ProcessSNPs(self.species, self.population, self.chrom)

    def output(self):
        return luigi.LocalTarget('db/{}-ensembl_vars.log'.format(self.basename))

    def run(self):
        dbc = db_conn(self.species)

        exec_time = dbc.execute_sql("""
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

        with self.output().open('w') as fout:
            fout.write('Execution took {}'.format(exec_time))


class LinkSNPChip(PipelineTask):
    """
    Link modern SNPs to their SNPChip variants

    :type species: str
    :type population: str
    :type chrom: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        return ProcessSNPs(self.species, self.population, self.chrom)

    def output(self):
        return luigi.LocalTarget('db/{}-snpchip.log'.format(self.basename))

    def run(self):
        dbc = db_conn(self.species)

        exec_time = dbc.execute_sql("""
            UPDATE modern_snps ms
              JOIN snpchip sc
                ON sc.chrom = ms.chrom
               AND sc.site = ms.site
               SET ms.snpchip_id = sc.id
             WHERE ms.population = '{pop}'
               AND ms.chrom = '{chrom}'""".format(pop=self.population, chrom=self.chrom))

        with self.output().open('w') as fout:
            fout.write('Execution took {}'.format(exec_time))


class DiscoverModernSNPs(luigi.WrapperTask):
    """
    Populate the modern_snps table and link records to genes, dbsnp and snpchip records.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):

        # process all the populations in chromosome chunks
        for pop in SAMPLES[self.species]:
            for chrom in CHROM_SIZE[self.species]:
                yield LinkEnsemblGenes(self.species, pop, chrom)
                yield LinkEnsemblVariants(self.species, pop, chrom)
                yield LinkSNPChip(self.species, pop, chrom)
