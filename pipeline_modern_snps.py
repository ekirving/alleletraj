#!/usr/bin/env python
# -*- coding: utf-8 -*-

import itertools
import luigi
import pysam

from collections import Counter

# import my custom modules
from pipeline_consts import SAMPLES, OUT_GROUP, CHROM_SIZE  # TODO make these into PipelineTask properties
from pipeline_snp_call import BiallelicSNPsVCF
from pipeline_utils import PipelineTask, PipelineExternalTask, PipelineWrapperTask
from db_conn import db_conn

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


def mutation_type(alleles):
    """
    Is this mutation a transition (A <-> G and C <-> T) or a transversion (everything else)
    """
    return 'ts' if set(alleles) == {'A', 'G'} or set(alleles) == {'C', 'T'} else 'tv'


class ExternalFASTA(PipelineExternalTask):
    """
    External task dependency for a whole-genome FASTA file.

    N.B. These have been created outside the workflow of this pipeline.

    :type sample: str
    """
    sample = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget("{}/{}/{}.fa".format(FASTA_PATH, self.sample, self.chrom))


class AlleleFrequencyFromFASTA(PipelineTask):
    """
    Extract all the SNPs from a given chromosome and estimate the allele frequency
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    @property
    def samples(self):
        """
        Include the outgroup in the sample list, as we need it for polarization
        """
        return [OUT_GROUP[self.species]] + SAMPLES[self.species][self.population]

    def requires(self):
        for sample in self.samples:
            return ExternalFASTA(sample)

    def output(self):
        return luigi.LocalTarget('db/{}-modern_snps.log'.format(self.basename))

    def run(self):
        # get all the input fasta files
        fasta_files = [fasta_file.path for fasta_file in self.input()]

        # open a private connection to the database
        dbc = db_conn(self.species)

        with self.output().open('w') as fout:
            fout.write("STARTED: Parsing {} fasta files.".format(len(fasta_files)))

            site = 0
            num_snps = 0
            data = []

            for fasta in fasta_files:
                # load the fasta data
                fin = open(fasta, "rU")

                # discard header row
                fin.readline()

                # wrap the file in an iterator
                data.append(stream_fasta(fin))

                fout.write("LOADED: {}".format(fasta))

            # zip the sequences together so we can iterate over them one site at a time
            for pileup in itertools.izip_longest(*data, fillvalue='N'):

                # increment the counter
                site += 1

                # convert iterator into a list
                pileup = list(pileup)

                # get the outgroup alleles
                out_allele = pileup.pop(0)

                # decode the fasta format
                haploids = decode_fasta(pileup)

                # how many alleles are there at this site
                num_alleles = len(set(haploids))

                # we only want biallelic SNPs
                if num_alleles > 2:
                    fout.write("WARNING: Polyallelic site chr{}:{} = {}".format(self.chrom, site, set(haploids)))
                    continue

                # count the haploid observations
                observations = Counter(haploids)

                # decode the outgroup allele
                ancestral = set(decode_fasta(out_allele))

                # we cannot handle variable sites in the outgroup
                if len(ancestral) != 1:
                    fout.write("WARNING: Unknown ancestral allele chr{}:{} = {}".format(self.chrom, site, out_allele))
                    continue

                ancestral = ancestral.pop()
                alleles = observations.keys()

                if ancestral not in alleles:
                    fout.write("WARNING: Polyallelic site chr{}:{} = {}, ancestral {}".format(self.chrom, site,
                                                                                              set(haploids), ancestral))
                    continue

                # is this mutation a transition or a transversion
                snp_type = mutation_type(alleles)

                alleles.remove(ancestral)
                derived = alleles.pop()

                # calculate the derived allele frequency (DAF)
                daf = float(observations[derived]) / (observations[ancestral] + observations[derived])

                record = {
                    'population': self.population,
                    'chrom': self.chrom,
                    'site': site,
                    'ancestral': ancestral,
                    'ancestral_count': observations[ancestral],
                    'derived': derived,
                    'derived_count': observations[derived],
                    'type': snp_type,
                    'daf': daf,
                }

                dbc.save_record('modern_snps', record)
                num_snps += 1

            fout.write("FINISHED: Added {:,} SNPs".format(num_snps))


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
        # open a private connection to the database
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
                haploids += [alleles for alleles in rec.samples[sample].alleles if alleles is not None]

            # count the haploid observations
            observations = Counter(haploids)
            alleles = observations.keys()

            if len(alleles) != 2:
                # skip non-variant sites
                continue

            # is this mutation a transition or a transversion
            snp_type = mutation_type(alleles)

            # calculate the derived allele frequency (DAF)
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


class ProcessSNPs(PipelineWrapperTask):
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
            yield AlleleFrequencyFromFASTA(self.species, self.population, self.chrom)
        else:
            # parse the VCF files for all other species
            yield AlleleFrequencyFromVCF(self.species, self.population, self.chrom)


class ModernSNPsPipeline(PipelineWrapperTask):
    """
    Populate the modern_snps table.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):

        # process SNPs for all populations and all chromosomes
        for pop in SAMPLES[self.species]:
            for chrom in CHROM_SIZE[self.species]:
                yield ProcessSNPs(self.species, pop, chrom)


if __name__ == '__main__':
    luigi.run()
