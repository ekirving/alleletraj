#!/usr/bin/env python
# -*- coding: utf-8 -*-

import luigi
import numpy

# import my custom modules
# TODO make these into PipelineTask properties
from pipeline_consts import BAM_FILES, SAMPLE_SEX, MIN_GENO_QUAL
from pipeline_utils import PipelineTask, PipelineExternalTask, PipelineWrapperTask, run_cmd

# VCF parser
from pysam import VariantFile

# quantiles for filtering VCF files
QUANTILE_LOW = 0.05
QUANTILE_HIGH = 0.95


class ExternalFASTA(PipelineExternalTask):
    """
    External task dependency for a reference assembly FASTA file.

    N.B. These have been downloaded outside the workflow of this pipeline.

    :type species: str
    """
    species = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget('fasta/{}.{}.dna.toplevel.fa'.format(self.binomial, self.assembly))


class ExternalBAM(PipelineExternalTask):
    """
    External task dependency for an aligned BAM file.

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
        yield ExternalFASTA(self.species)

        # include the outgroup in the sample list, as we need it for polarization
        for sample in [self.outgroup] + self.samples:
            yield ExternalBAM(self.species, sample)

    def output(self):
        return luigi.LocalTarget('vcf/{}.vcf.gz'.format(self.basename))

    def run(self):
        # unpack the input params
        ref_file, bam_files = self.input()[0], self.input()[1:]

        # bcftools needs the sex specified in a separate file
        sex_file = 'vcf/{}_{}.sex'.format(self.species, self.population)

        with open(sex_file, 'w') as fout:
            for sample in [self.outgroup] + self.samples:
                fout.write('{}\t{}\n'.format(sample, SAMPLE_SEX[self.species][sample]))

        with self.output().temporary_path() as vcf_out:
            params = {
                'ref': ref_file.path,
                'chr': self.chrom,
                'bam': ' '.join([bam.path for bam in bam_files]),
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
    Calculate the lower (5%) and upper (95%) quantiles of the depth of coverage for a VCF.

    The DoC distribution is calculated exclusively on sites that pass the genotype quality filter.

    :type species: str
    :type population: str
    :type chrom: str
    :type qual: int
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()
    qual = luigi.IntParameter()

    def requires(self):
        return BCFToolsCall(self.species, self.population, self.chrom)

    def output(self):
        return luigi.LocalTarget('vcf/{}.quant'.format(self.basename))

    def run(self):
        depth = []

        # iterate over the VCF and extract the depth of coverage at each site
        for rec in VariantFile(self.input().path).fetch():
            try:
                # only count depth at sites which pass the genotype quality filter
                if int(rec.qual) >= self.qual:
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
    Remove low quality sites and sites in the upper and lower quantiles of the depth of coverage distribution.

    Also, normalise indels and merge multiallelic sites.

    :type species: str
    :type population: str
    :type chrom: str
    :type qual: int
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()
    qual = luigi.IntParameter()

    def requires(self):
        yield ExternalFASTA(self.species)
        yield BCFToolsCall(self.species, self.population, self.chrom)
        yield QuantilesOfCoverageVCF(self.species, self.population, self.chrom, self.qual)

    def output(self):
        return luigi.LocalTarget('vcf/{}-quant.vcf.gz'.format(self.basename))

    def run(self):

        # unpack the input params
        ref_file, vcf_input, quant_file = self.input()

        # get the quantiles
        qlow, qhigh = numpy.loadtxt(quant_file.path)

        with self.output().temporary_path() as vcf_out:
            params = {
                'qual':  MIN_GENO_QUAL,
                'qlow':  int(qlow),
                'qhigh': int(qhigh),
                'vcf':   vcf_input.path,
                'ref':   ref_file.path,
                'out':   vcf_out,
            }

            cmd = "bcftools filter --exclude 'QUAL<{qual} | DP<{qlow} | DP>{qhigh}' --output-type u {vcf} | " \
                  "bcftools norm --fasta-ref {ref} --multiallelics +any --output-type z --output {out}".format(**params)

            run_cmd([cmd], shell=True)


class PolarizeVCF(PipelineTask):
    """
    Switch the REF allele for the ancestral allele, based on an outgroup present in the VCF.

    Drops any sites which are heterozygous or uncallable in the outgroup.

    :type species: str
    :type population: str
    :type chrom: str
    :type qual: int
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()
    qual = luigi.IntParameter(default=MIN_GENO_QUAL)

    def requires(self):
        return FilterVCF(self.species, self.population, self.chrom, self.qual)

    def output(self):
        return luigi.LocalTarget('vcf/{}-quant-polar.vcf.gz'.format(self.basename))

    def run(self):

        with self.output().temporary_path() as tmp_out:
            # open both VCF files
            vcf_in = VariantFile(self.input().path)
            vcf_out = VariantFile(tmp_out, 'wz', header=vcf_in.header)  # `z` tells pysam to bgzip output

            # iterate over the VCF and determine the ancestral allele
            for rec in vcf_in.fetch():

                # get the outgroup alleles
                out_alleles = rec.samples[self.outgroup].alleles

                # skip sites in the outgroup that are heterozygous or missing
                if len(set(out_alleles)) != 1 or out_alleles[0] is None:
                    continue

                # get the ancestral allele
                anc = out_alleles[0]

                # do we need to polarize this site
                if rec.ref != anc:

                    # get all the alleles at this site, minus the ancestral
                    alt = set([rec.ref] + list(rec.alts))
                    alt.remove(anc)

                    # VCFs store the GT as an allele index, so we have to update the indices
                    alleles = [anc] + list(alt)  # NOTE important distinction between [] and list()
                    indices = dict(zip(alleles, range(0, len(alleles))))

                    for sample in rec.samples:
                        rec.samples[sample].allele_indices = [indices.get(gt, None)
                                                              for gt in rec.samples[sample].alleles]

                    # polarize the REF/ALT alleles
                    rec.ref = anc
                    rec.alts = alt

                vcf_out.write(rec)


class BiallelicSNPsVCF(PipelineTask):
    """
    Extract all the biallelic SNPs from the filtered and polarised VCF.

    :type species: str
    :type population: str
    :type chrom: str
    :type qual: int
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()
    qual = luigi.IntParameter(default=MIN_GENO_QUAL)

    def requires(self):
        return PolarizeVCF(self.species, self.population, self.chrom, self.qual)

    def output(self):
        return luigi.LocalTarget('vcf/{}-quant-polar-SNPs.vcf.gz'.format(self.basename))

    def run(self):
        with self.output().temporary_path() as vcf_out:
            run_cmd(['bcftools',
                     'view',
                     '--types', 'snps',            # only keep SNPs
                     '--min-alleles', 2,           # which are biallelic
                     '--max-alleles', 2,
                     '--exclude', 'INFO/INDEL=1',  # exclude sites marked as INDELs in INFO tag
                     '--samples', '^' + self.outgroup,  # drop the outgroup
                     '--min-ac', '1:nref',         # exclude sites exclusively hom-ALT, as these are likely mispolarised
                     '--output-type', 'z',
                     '--output-file', vcf_out,
                     self.input().path])

        # index the vcf
        run_cmd(['bcftools', 'index', '--tbi', self.output().path])


class WholeGenomeSNPsVCF(PipelineTask):
    """
    Concatenate the all chromosome VCFs into a single file.

    :type species: str
    :type population: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()

    def requires(self):
        for chrom in self.chromosomes:
            # TODO handle chr prefixes better
            yield BiallelicSNPsVCF(self.species, self.population, 'chr{}'.format(chrom))

    def output(self):
        return luigi.LocalTarget('vcf/{}-chrAll-filtered-polar-SNPs.vcf.gz'.format(self.basename))

    def run(self):

        # get all the vcf files to concatenate
        vcf_files = [vcf.path for vcf in self.input()]

        with self.output().temporary_path() as tmp_out:
            run_cmd(['bcftools',
                     'concat',
                     '--output-type', 'z',
                     '--output', tmp_out
                     ] + vcf_files)

        # index the vcf
        run_cmd(['bcftools', 'index', '--tbi', self.output().path])


class SNPCallPipeline(PipelineWrapperTask):
    """
    Call SNPs using the bcftools `mpileup | call` workflow.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        for pop in self.populations:
            yield WholeGenomeSNPsVCF(self.species, pop)


if __name__ == '__main__':
    luigi.run()
