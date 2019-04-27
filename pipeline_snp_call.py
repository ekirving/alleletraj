#!/usr/bin/env python
# -*- coding: utf-8 -*-

import luigi
import numpy

# VCF parser
from pysam import VariantFile

# import my custom modules
from pipeline_alignment import ReferenceFASTA, AlignedBAM
from pipeline_consts import SAMPLE_SEX
from pipeline_utils import PipelineTask, PipelineExternalTask, PipelineWrapperTask, run_cmd

# the minimum phred scaled genotype quality (30 = 99.9%)
MIN_GENO_QUAL = 30

# quantiles for filtering VCF files
QUANTILE_LOW = 0.05
QUANTILE_HIGH = 0.95


class ReferencePloidy(PipelineExternalTask):
    """
    Make a ploidy-file file for bcftools which specifies the sex based ploidy of chromosomes in an assembly.

    See --ploidy-file in https://samtools.github.io/bcftools/bcftools.html#call

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return ReferenceFASTA(self.species)

    def output(self):
        return luigi.LocalTarget('ensembl/{}.{}.dna.toplevel.ploidy'.format(self.binomial, self.assembly))

    def run(self):
        # get the reference index
        _, fai_file = self.input()

        # build the ploidy file
        with self.output().open('w') as fout:

            # by iterating over the reference index
            with fai_file.open('r') as fin:

                for line in fin:
                    # get the chrom name and size
                    contig, size, _, _, _ = line.split()

                    # strip any annoying `chr` prefix
                    nochr = contig.replace('chr', '')

                    if nochr == 'X':
                        for sex, ploidy in [('M', '1'), ('F', '2')]:
                            fout.write('\t'.join([contig, '1', size, sex, ploidy]) + '\n')
                    elif nochr == 'Y':
                        for sex, ploidy in [('M', '1'), ('F', '0')]:
                            fout.write('\t'.join([contig, '1', size, sex, ploidy]) + '\n')
                    elif nochr == 'MT':
                        for sex, ploidy in [('M', '1'), ('F', '1')]:
                            fout.write('\t'.join([contig, '1', size, sex, ploidy]) + '\n')

                # all other chroms are diploid
                for sex in ['M', 'F']:
                    fout.write('\t'.join(['*', '*', '*', sex, '2']) + '\n')


class BCFToolsCall(PipelineTask):
    """
    Make genotype calls using the bcftools mpileup workflow

    :type species: str
    :type chrom: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        yield ReferenceFASTA(self.species)
        yield ReferencePloidy(self.species)

        for sample in self.all_samples:
            yield AlignedBAM(self.species, sample)

    def output(self):
        return luigi.LocalTarget('vcf/{}.vcf.gz'.format(self.basename))

    def run(self):
        # unpack the input params
        (ref_file, _), pld_file, bam_files = self.input()[0], self.input()[1], self.input()[2:]

        # bcftools needs the sex specified in a separate file
        sex_file = 'vcf/{}-modern.sex'.format(self.species)
        with open(sex_file, 'w') as fout:
            for sample in self.all_samples:
                fout.write('{}\t{}\n'.format(sample, SAMPLE_SEX[self.species][sample]))

        # sample names in the BAM file(s) may not be consistent, so override the @SM code
        rgs_file = 'vcf/{}-modern.rgs'.format(self.species)
        with open(rgs_file, 'w') as fout:
            for idx, sample in enumerate(self.all_samples):
                fout.write('*\t{}\t{}\n'.format(bam_files[idx].path, sample))

        with self.output().temporary_path() as vcf_out:
            params = {
                'ref': ref_file.path,
                'chr': 'chr{}'.format(self.chrom) if self.species == 'horse' else self.chrom,
                'rgs': rgs_file,
                'bam': ' '.join([bam.path for bam in bam_files]),
                'pld': pld_file.path,
                'sex': sex_file,
                'vcf': vcf_out
            }

            cmd = "bcftools mpileup --fasta-ref {ref} --regions {chr} --read-groups {rgs} --output-type u {bam} | " \
                  "bcftools call  --multiallelic-caller --ploidy-file {pld} --samples-file {sex} --output-type z " \
                  " --output {vcf}".format(**params)

            run_cmd([cmd], shell=True)


class QuantilesOfCoverageVCF(PipelineTask):
    """
    Calculate the lower (5%) and upper (95%) quantiles of the depth of coverage for a VCF.

    The DoC distribution is calculated exclusively on sites that pass the genotype quality filter.

    :type species: str
    :type chrom: str
    :type qual: int
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()
    qual = luigi.IntParameter()

    def requires(self):
        return BCFToolsCall(self.species, self.chrom)

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
    :type chrom: str
    :type qual: int
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()
    qual = luigi.IntParameter()

    def requires(self):
        yield ReferenceFASTA(self.species)
        yield BCFToolsCall(self.species, self.chrom)
        yield QuantilesOfCoverageVCF(self.species, self.chrom, self.qual)

    def output(self):
        return luigi.LocalTarget('vcf/{}-quant.vcf.gz'.format(self.basename))

    def run(self):
        # unpack the input params
        (ref_file, _), vcf_input, quant_file = self.input()

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
    :type chrom: str
    :type qual: int
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()
    qual = luigi.IntParameter(default=MIN_GENO_QUAL)

    def requires(self):
        return FilterVCF(self.species, self.chrom, self.qual)

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
    :type chrom: str
    :type qual: int
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()
    qual = luigi.IntParameter(default=MIN_GENO_QUAL)

    def requires(self):
        return PolarizeVCF(self.species, self.chrom, self.qual)

    def output(self):
        return luigi.LocalTarget('vcf/{}-quant-polar-SNPs.vcf.gz'.format(self.basename))

    def run(self):
        with self.output().temporary_path() as vcf_out:
            run_cmd(['bcftools',
                     'view',
                     '--samples', '^' + self.outgroup,  # drop the outgroup
                     '--types', 'snps',                 # only keep SNPs
                     '--min-alleles', 2,                # which are biallelic
                     '--max-alleles', 2,
                     '--min-ac', '1:minor',             # and they must be variable, after dropping the outgroup
                     '--exclude', 'INFO/INDEL=1',       # exclude sites marked as INDELs in INFO tag
                     '--output-type', 'z',
                     '--output-file', vcf_out,
                     self.input().path])

        # index the vcf
        run_cmd(['bcftools', 'index', '--tbi', self.output().path])


class WholeGenomeSNPsVCF(PipelineTask):
    """
    Concatenate the all chromosome VCFs into a single file.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        for chrom in self.chromosomes:
            yield BiallelicSNPsVCF(self.species, chrom)

    def output(self):
        return luigi.LocalTarget('vcf/{}-chrAll-filtered-polar-SNPs.vcf.gz'.format(self.species))

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
        yield WholeGenomeSNPsVCF(self.species)


if __name__ == '__main__':
    luigi.run()
