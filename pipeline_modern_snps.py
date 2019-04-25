#!/usr/bin/env python
# -*- coding: utf-8 -*-

import itertools
import luigi
import pysam

from collections import Counter

# import my custom modules
from pipeline_database import CreateDatabase
from pipeline_snp_call import BiallelicSNPsVCF
from pipeline_utils import PipelineTask, PipelineExternalTask, PipelineWrapperTask

# TODO move into spreadsheet
FASTA_PATH = '/media/jbod2/raid1-sdc1/laurent/full_run_results/Pig/modern/FASTA'

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
    External task dependency for a chromosome FASTA file.

    N.B. These have been created outside the workflow of this pipeline.

    :type sample: str
    :type chrom: str
    """
    sample = luigi.Parameter()
    chrom = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget("{}/{}/{}.fa".format(FASTA_PATH, self.sample, self.chrom))


class ModernSNPsFromFASTA(PipelineTask):
    """
    Ascertain moderns SNPS, and estimate their allele frequency, from the given chromosome FASTA files.

    :type species: str
    :type chrom: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        yield CreateDatabase(self.species)

        # include the outgroup in the sample list, as we need it for polarization
        for sample in self.all_samples:
            yield ExternalFASTA(sample)

    def output(self):
        return luigi.LocalTarget('db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        # get all the input fasta files
        fasta_files = [fasta_file.path for fasta_file in self.input()[1:]]

        # open a private connection to the database
        dbc = self.db_conn()

        with self.output().open('w') as fout:
            fout.write("STARTED: Parsing {} fasta files.\n".format(len(fasta_files)))

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

                fout.write("LOADED: {}\n".format(fasta))

            # zip the sequences together so we can iterate over them one site at a time
            for pileup in itertools.izip_longest(*data, fillvalue='N'):

                # increment the counter
                site += 1

                # convert iterator into a list
                pileup = list(pileup)

                # decode the outgroup allele(s)
                out_alleles = decode_fasta(pileup.pop(0))

                # we cannot handle variable sites in the outgroup
                if len(set(out_alleles)) != 1:
                    fout.write("WARNING: Cannot polarize site chr{}:{} = {}\n".format(self.chrom, site, out_alleles))
                    continue

                # get the ancestral allele
                ancestral = out_alleles.pop()

                # TODO how to group by population
                # decode all the sample genotypes
                haploids = decode_fasta(pileup)

                all_alleles = set(out_alleles + haploids)

                # how many alleles are there at this site
                num_alleles = len(all_alleles)

                # we only want biallelic SNPs
                if num_alleles == 1:
                    # skip non-variant sites
                    continue
                elif num_alleles > 2:
                    fout.write("WARNING: Polyallelic site chr{}:{} = {}\n".format(self.chrom, site, all_alleles))
                    continue
                elif len(set(haploids)) == 1:
                    # exclude sites exclusively hom-ALT, as these are likely mispolarised
                    fout.write("WARNING: Mispolarized site chr{}:{} = {}\n".format(self.chrom, site, all_alleles))
                    continue

                # count the haploid observations
                observations = Counter(haploids)
                alleles = observations.keys()

                # is this mutation a transition or a transversion
                snp_type = mutation_type(alleles)

                alleles.remove(ancestral)
                derived = alleles.pop()

                # calculate the derived allele frequency (DAF)
                daf = float(observations[derived]) / (observations[ancestral] + observations[derived])

                modern_snp = {
                    'chrom': self.chrom,
                    'site': site,
                    'ancestral': ancestral,
                    'derived': derived,
                    'type': snp_type,
                }

                modsnp_id = dbc.save_record('modern_snps', modern_snp)
                num_snps += 1

                modern_snp_daf = {
                    'modsnp_id': modsnp_id,
                    'ancestral_count': observations[ancestral],
                    'derived_count': observations[derived],
                    'daf': daf,
                }

                dbc.save_record('modern_snp_daf', modern_snp_daf)

            fout.write("FINISHED: Added {:,} SNPs\n".format(num_snps))


class ModernSNPsFromVCF(PipelineTask):
    """
    Ascertain moderns SNPS, and estimate their allele frequency, from the given chromosome VCF file.

    :type species: str
    :type chrom: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        yield CreateDatabase(self.species)
        yield BiallelicSNPsVCF(self.species, self.chrom)

    def output(self):
        return luigi.LocalTarget('db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        # unpack the inputs
        _, vcf_file = self.input()

        # open a private connection to the database
        dbc = self.db_conn()

        # count the number of SNPs added
        num_snps = 0

        # parse the VCF with pysam
        for rec in pysam.VariantFile(vcf_file.path).fetch():

            # the VCF has already been polarised, so the REF/ALT are the ancestral/derived
            ancestral = rec.ref
            derived = rec.alts[0]

            # is this mutation a transition or a transversion
            snp_type = mutation_type([ancestral, derived])

            modern_snp = {
                'chrom': self.chrom,
                'site': rec.pos,
                'ancestral': ancestral,
                'derived': derived,
                'type': snp_type,
            }

            modsnp_id = dbc.save_record('modern_snps', modern_snp)

            num_snps += 1

            # resolve the genotypes of the samples in one population at a time
            for pop in self.populations:

                # lets collate all the haploid observations for the two alleles
                haploids = []

                for sample in self.populations[pop]:
                    # get the alleles, but skip any missing genotypes
                    haploids += [alleles for alleles in rec.samples[sample].alleles if alleles is not None]

                # count the haploid observations
                observations = Counter(haploids)

                # calculate the derived allele frequency (DAF)
                daf = float(observations[derived]) / (observations[ancestral] + observations[derived])

                modern_snp_daf = {
                    'modsnp_id': modsnp_id,
                    'ancestral_count': observations[ancestral],
                    'derived_count': observations[derived],
                    'daf': daf,
                }

                dbc.save_record('modern_snp_daf', modern_snp_daf)

        with self.output().open('w') as fout:
            fout.write("Added {:,} SNPs".format(num_snps))


class LoadModernSNPs(PipelineWrapperTask):
    """
    Ascertain modern SNPs in whole-genome data.

    :type species: str
    :type chrom: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        if self.species == 'pig':
            # special case of FASTA data for pig ascertainment
            yield ModernSNPsFromFASTA(self.species, self.chrom)
        else:
            # parse the VCF files for all other species
            yield ModernSNPsFromVCF(self.species, self.chrom)


class ModernSNPsPipeline(PipelineWrapperTask):
    """
    Populate the modern_snps table.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):

        # process SNPs for all chromosomes
        for chrom in self.chromosomes:
            yield LoadModernSNPs(self.species, chrom)


if __name__ == '__main__':
    luigi.run()
