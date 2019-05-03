#!/usr/bin/env python
# -*- coding: utf-8 -*-

import luigi
import numpy

# VCF parser
from pysam import VariantFile

# import my custom modules
from alleletraj.modern.alignment import ReferenceFASTA, AlignedBAM
from alleletraj import utils

# the minimum phred scaled genotype quality (30 = 99.9%)
MIN_GENO_QUAL = 30

# quantiles for filtering VCF files
QUANTILE_LOW = 0.05
QUANTILE_HIGH = 0.95


class ReferencePloidy(utils.PipelineExternalTask):
    """
    Make a ploidy-file file for bcftools which specifies the sex based ploidy of chromosomes in an assembly.

    See --ploidy-file in https://samtools.github.io/bcftools/bcftools.html#call

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return ReferenceFASTA(self.species)

    def output(self):
        return luigi.LocalTarget('data/ensembl/{}.{}.dna.toplevel.ploidy'.format(self.binomial, self.assembly))

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


class BCFToolsSamplesFile(utils.PipelineTask):
    """
    Make a samples sex file and a samples readgroup file for bcftools.

    The sex file specifies the sex of each samples, to be used with the ploidy file, and the readgroup file ensures that
    BAM files containing multiple libraries with different accession codes are treated as one sample.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        for pop, sample in self.all_samples:
            yield AlignedBAM(self.species, pop, sample)

    def output(self):
        return [luigi.LocalTarget('data/vcf/{}.{}'.format(self.species, ext)) for ext in ['sex', 'rgs']]

    def run(self):
        # unpack the params
        bam_files = [bam_file for bam_file, bai_file in self.input()]
        sex_file, rgs_file = self.output()

        # bcftools needs the sex specified in a separate file
        with sex_file.open('w') as fout:
            for pop, sample in self.all_samples:
                sex = self.all_populations[pop][sample]['sex'][0].upper()  # use single letter uppercase
                fout.write('{}\t{}\n'.format(sample, sex))

        # sample names in the BAM file(s) may not be consistent, so override the @SM code
        with rgs_file.open('w') as fout:
            for idx, (pop, sample) in enumerate(self.all_samples):
                fout.write('*\t{}\t{}\n'.format(bam_files[idx].path, sample))


class BCFToolsCall(utils.PipelineTask):
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
        yield BCFToolsSamplesFile(self.species)

        for pop, sample in self.all_samples:
            yield AlignedBAM(self.species, pop, sample)

    def output(self):
        return luigi.LocalTarget('data/vcf/{}.vcf.gz'.format(self.basename))

    def run(self):
        # unpack the input params
        ref_file, _ = self.input()[0]
        pld_file = self.input()[1]
        sex_file, rgs_file = self.input()[2]
        bam_files = [bam_file for bam_file, bai_file in self.input()[3:]]

        with self.output().temporary_path() as vcf_out:
            params = {
                'ref': ref_file.path,
                'chr': self.chrom,
                'rgs': rgs_file.path,
                'bam': ' '.join([bam.path for bam in bam_files]),
                'pld': pld_file.path,
                'sex': sex_file.path,
                'vcf': vcf_out
            }

            cmd = "bcftools mpileup --fasta-ref {ref} --regions {chr} --read-groups {rgs} --output-type u {bam} | " \
                  "bcftools call  --multiallelic-caller --ploidy-file {pld} --samples-file {sex} --output-type z " \
                  " --output {vcf}".format(**params)

            utils.run_cmd([cmd], shell=True)


class QuantilesOfCoverageVCF(utils.PipelineTask):
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
        return luigi.LocalTarget('data/vcf/{}.quant'.format(self.basename))

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


class FilterVCF(utils.PipelineTask):
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
        return luigi.LocalTarget('data/vcf/{}-quant.vcf.gz'.format(self.basename))

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

            utils.run_cmd([cmd], shell=True)


class PolarizeVCF(utils.PipelineTask):
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
        return luigi.LocalTarget('data/vcf/{}-quant-polar.vcf.gz'.format(self.basename))

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


class BiallelicSNPsVCF(utils.PipelineTask):
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
        return luigi.LocalTarget('data/vcf/{}-quant-polar-SNPs.vcf.gz'.format(self.basename))

    def run(self):
        with self.output().temporary_path() as vcf_out:
            utils.run_cmd(['bcftools',
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
        utils.run_cmd(['bcftools', 'index', '--tbi', self.output().path])


class WholeGenomeSNPsVCF(utils.PipelineTask):
    """
    Concatenate the all chromosome VCFs into a single file.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        for chrom in self.chromosomes:
            yield BiallelicSNPsVCF(self.species, chrom)

    def output(self):
        return luigi.LocalTarget('data/vcf/{}-chrAll-filtered-polar-SNPs.vcf.gz'.format(self.species))

    def run(self):

        # get all the vcf files to concatenate
        vcf_files = [vcf.path for vcf in self.input()]

        with self.output().temporary_path() as tmp_out:
            utils.run_cmd(['bcftools',
                           'concat',
                           '--output-type', 'z',
                           '--output', tmp_out
                           ] + vcf_files)

        # index the vcf
        utils.run_cmd(['bcftools', 'index', '--tbi', self.output().path])


class SNPCallPipeline(utils.PipelineWrapperTask):
    """
    Call SNPs using the bcftools `mpileup | call` workflow.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        yield WholeGenomeSNPsVCF(self.species)


if __name__ == '__main__':
    luigi.run()
