#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
from collections import Counter

# third party modules
import luigi
import pysam

# local modules
from alleletraj import utils
from alleletraj.db.load import CreateDatabase
from vcf import BiallelicSNPsVCF


def mutation_type(alleles):
    """
    Is this mutation a transition (A <-> G and C <-> T) or a transversion (everything else)
    """
    return 'ts' if set(alleles) == {'A', 'G'} or set(alleles) == {'C', 'T'} else 'tv'


class LoadModernSNPs(utils.PipelineTask):
    """
    Ascertain moderns SNPs and estimate their allele frequency from a chromosome VCF file.

    :type species: str
    :type chrom: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        yield CreateDatabase(self.species)
        yield BiallelicSNPsVCF(self.species, self.chrom)

    def output(self):
        return luigi.LocalTarget('data/db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        # unpack the inputs
        _, vcf_file = self.input()

        # open a private connection to the db
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
            for pop in self.modern_pops:

                # lets collate all the haploid observations for the two alleles
                haploids = []

                for sample in self.modern_pops[pop]:
                    # get the alleles, but skip any missing genotypes
                    haploids += [alleles for alleles in rec.samples[sample].alleles if alleles is not None]

                if len(haploids) == 0:
                    # skip sites with no genotypes for this population
                    continue

                # count the haploid observations
                observations = Counter(haploids)

                # calculate the derived allele frequency (DAF)
                daf = float(observations[derived]) / (observations[ancestral] + observations[derived])

                modern_snp_daf = {
                    'modsnp_id': modsnp_id,
                    'population': pop,
                    'ancestral_count': observations[ancestral],
                    'derived_count': observations[derived],
                    'daf': daf,
                }

                dbc.save_record('modern_snp_daf', modern_snp_daf)

        with self.output().open('w') as fout:
            fout.write('Added {:,} SNPs'.format(num_snps))


class ModernSNPsPipeline(utils.PipelineWrapperTask):
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
