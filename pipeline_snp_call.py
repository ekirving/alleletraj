#!/usr/bin/env python
# -*- coding: utf-8 -*-

import luigi
import numpy

# import my custom modules
from pipeline_consts import *
from pipeline_utils import PipelineTask, run_cmd

# VCF parser
from pysam import VariantFile


class BAMfile(luigi.ExternalTask):
    """
    External task dependency for the aligned BAM file.

    N.B. These have been created outside the workflow of this pipeline.

    :type species: str
    :type sample: str
    """
    species = luigi.Parameter()
    sample = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(BAM_FILES[self.species][self.sample])


class BCFToolsCall(PipelineTask):
    """
    Make genotype calls using the bcftools mpileup workflow

    :type species: str
    :type population: str
    :type chrom: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        for sample in SAMPLES[self.species][self.population]:
            yield BAMfile(self.species, sample)

    def output(self):
        return luigi.LocalTarget('vcf/{}.vcf.gz'.format(self.basename))

    def run(self):

        # bcftools needs the sex specified in a separate file
        sex_file = 'data/{}_{}.sex'.format(self.species, self.population)

        with open(sex_file, 'w') as fout:
            for sample in SAMPLES[self.species][self.population]:
                fout.write('{}\t{}\n'.format(sample, SAMPLE_SEX[self.species][sample]))

        with self.output().temporary_path() as vcf_out:
            params = {
                'ref': REF_FILE[self.species],
                'chr': self.chrom,
                'bam': ' '.join([bam.path for bam in self.input()]),
                'pld': 'data/{}.ploidy'.format(self.species),
                'sex': sex_file,
                'vcf': vcf_out
            }

            cmd = "bcftools mpileup --fasta-ref {ref} --regions {chr} --output-type u {bam} | " \
                  "bcftools call  --multiallelic-caller --ploidy-file {pld} --samples-file {sex} " \
                  " --output-type z --output {vcf}".format(**params)

            run_cmd([cmd], shell=True)


class QuantilesOfCoverageVCF(PipelineTask):
    """
    Calculate the upper and lower quantiles of the depth of coverage for a VCF.

    :type species: str
    :type population: str
    :type chrom: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        return BCFToolsCall(self.species, self.population, self.chrom)

    def output(self):
        return luigi.LocalTarget('vcf/{}.DoC'.format(self.basename))

    def run(self):
        depth = []

        # iterate over the VCF and extract the depth of coverage at each site
        for rec in VariantFile(self.input().path).fetch():
            try:
                depth.append(rec.info['DP'])
            except KeyError:
                pass

        # calculate the quantiles
        quants = numpy.quantile(depth, [QUANTILE_LOW, QUANTILE_HIGH])

        # save them to disk
        with self.output().open('w') as fout:
            fout.write('{} {}'.format(int(quants[0]), int(quants[1])))


class FilterVCF(PipelineTask):
    """
    Remove low quality sites and those in the upper and lower quantiles of the depth of coverage distribution.

    :type species: str
    :type population: str
    :type chrom: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        yield BCFToolsCall(self.species, self.population, self.chrom)
        yield QuantilesOfCoverageVCF(self.species, self.population, self.chrom)

    def output(self):
        return luigi.LocalTarget('vcf/{}-filtered.vcf.gz'.format(self.basename))

    def run(self):

        # unpack the input params
        vcf_input, quant_file = self.input()

        # get the quantiles
        qlow, qhigh = numpy.loadtxt(quant_file.path)

        with self.output().temporary_path() as vcf_out:
            run_cmd(['bcftools',
                     'filter',
                     '--exclude', 'QUAL<{} | DP<{} | DP>{}'.format(MIN_GENO_QUAL, int(qlow), int(qhigh)),
                     '--output-type', 'z',
                     '--output', vcf_out,
                     vcf_input.path])


class ConcatFilteredVCFs(PipelineTask):
    """
    Concatenate the chromosome level VCFs into a single file.

    :type species: str
    :type population: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()

    def requires(self):
        for chrom in CHROM_SIZE[self.species]:
            yield FilterVCF(self.species, self.population, 'chr{}'.format(chrom))

    def output(self):
        return luigi.LocalTarget('vcf/{}-chrAll-filtered.vcf.gz'.format(self.basename))

    def run(self):

        # unpack the input params
        vcf_files = [vcf.path for vcf in self.input()]

        with self.output().temporary_path() as vcf_out:
            run_cmd(['bcftools',
                     'concat',
                     '--output-type', 'z',
                     '--output', vcf_out,
                     ] + vcf_files)


class PolarizeVCF(PipelineTask):
    """
    Switch the REF allele for the ancestral allele, based on an outgroup present in the VCF.

    Drops any sites which are heterozygous or uncallable in the outgroup.

    :type species: str
    :type population: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()

    def requires(self):
        return ConcatFilteredVCFs(self.species, self.population)

    def output(self):
        return luigi.LocalTarget('vcf/{}-chrAll-filtered-polar.vcf.gz'.format(self.basename))

    def run(self):

        # open both VCF files
        vcf_in = VariantFile(self.input().path)
        vcf_out = VariantFile(self.output().path, 'w', header=vcf_in.header)

        # iterate over the VCF and determine the ancestral allele
        for rec in vcf_in.fetch():

            # get the outgroup alleles
            out_alleles = rec.samples[OUTGROUP].alleles

            # skip sites that are either heterozygous or missing in the outgroup
            if len(set(out_alleles)) != 1 or out_alleles[0] is None:
                continue

            # get the ancestral allele
            anc = out_alleles[0]

            # do we need to polarize this site
            if rec.ref != anc:

                # get all the alleles at this site, minus the ancestral
                alt = set(tuple(rec.ref) + tuple(rec.alts))
                alt.remove(anc)

                # VCFs store the GT as an allele index, so we have to update the indices
                alleles = list(anc) + list(alt)
                indices = dict(zip(alleles, range(0, len(alleles))))

                for sample in rec.samples:
                    rec.samples[sample].allele_indices = [indices.get(gt, None) for gt in rec.samples[sample].alleles]

                # polarize the REF/ALT alleles
                rec.ref = anc
                rec.alts = alt

            vcf_out.write(rec)


class SubsetSNPsVCF(PipelineTask):
    """
    Extract all the SNPs from the VCF.

    :type species: str
    :type population: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()

    def requires(self):
        return PolarizeVCF(self.species, self.population)

    def output(self):
        return luigi.LocalTarget('vcf/{}-chrAll-filtered-polar-SNPs.vcf.gz'.format(self.basename))

    def run(self):

        with self.output().temporary_path() as vcf_out:
            run_cmd(['bcftools',
                     'view',
                     '--types', 'snps',
                     '--exclude', 'INFO/INDEL=1',
                     '--output-type', 'z',
                     '--output', vcf_out,
                     self.input().path])


class BCFtoolsCallSNPs(luigi.WrapperTask):
    """
    Call SNPs using the bcftools mpileup | call workflow.
    """

    def requires(self):
        return SubsetSNPsVCF('horse', 'DOM2')


if __name__ == '__main__':
    luigi.run()
