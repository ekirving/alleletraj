#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import my custom modules
from pipeline_utils import *

# VCF parser
from pysam import VariantFile


class ModernVCF(luigi.ExternalTask):
    """
    Make sure the VCF files exist

    :type species: str
    :type population: str
    :type chrom: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget("vcf/{}-{}-{}.vcf".format(self.species, self.population, self.chrom))


class DepthOfCoverage(PipelineTask):
    """
    Extract the depth of coverage for every site in a VCF.

    :type species: str
    :type population: str
    :type chrom: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        # TODO link back to task that made the VCF
        return ModernVCF(self.species, self.population, self.chrom)

    def output(self):
        return luigi.LocalTarget("vcf/{}.depth".format(self.basename))

    def run(self):
        with self.output().open('w') as fout:

            # iterate over the VCF and extract the depth of coverage at each site
            for rec in VariantFile(self.input().path).fetch():
                try:
                    fout.write('{}\n'.format(rec.info['DP']))
                except KeyError:
                    pass


class QuantilesOfCoverage(PipelineTask):
    """
    Calculate the quantiles of the depth of coverage for a VCF.

    :type species: str
    :type population: str
    :type chrom: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        return DepthOfCoverage(self.species, self.population, self.chrom)

    def output(self):
        return luigi.LocalTarget("vcf/{}.quant".format(self.basename))

    def run(self):
        depth = np.loadtxt(self.input().path)
        quant = np.quantile(depth, [QUANTILE_LOW, QUANTILE_HIGH])
        np.savetxt(self.output().path, quant)


class FilterQuantiles(PipelineTask):
    """
    Remove sites in the upper and lower quantiles of the depth of coverage distribution.

    :type species: str
    :type population: str
    :type chrom: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        yield ModernVCF(self.species, self.population, self.chrom)
        yield QuantilesOfCoverage(self.species, self.population, self.chrom)

    def output(self):
        return luigi.LocalTarget("vcf/{}.filtered.vcf".format(self.basename))

    def run(self):

        # unpack the input params
        vcf_input, quant_file = self.input()

        # get the quantiles
        qlow, qhigh = np.loadtxt(quant_file.path)

        with self.output().open('w') as vcf_out:

            run_cmd(["bcftools",
                     "filter",
                     "--exclude", "DP<{} | DP>{}".format(int(qlow), int(qhigh)),
                     vcf_input.path], stdout=vcf_out)


class MergeVCF(PipelineTask):
    """
    Merge the chromosome level VCFs into a single file.

    :type species: str
    :type population: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()

    def requires(self):
        # TODO replace with CHROM[self.species]
        for chrom in ['chr1', 'chr2']:
            yield FilterQuantiles(self.species, self.population, chrom)

    def output(self):
        return luigi.LocalTarget("vcf/{}.chrAll.filtered.vcf".format(self.basename))

    def run(self):

        # unpack the input params
        vcf_files = [vcf.path for vcf in self.input()]

        with self.output().open('w') as vcf_out:
            run_cmd(["bcftools", "concat"] + vcf_files, stdout=vcf_out)


class EasySFS(PipelineTask):
    """
    Calculate the Site Frequency Spectrum.

    :type species: str
    :type population: str
    :type chrom: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        yield FilterQuantiles(self.species, self.population, self.chrom)

    def output(self):
        return luigi.LocalTarget("dadi/{}.sfs".format(self.basename))

    def run(self):

        # unpack the input params
        vcf_input, quant_file = self.input()

        # get the quantiles
        qlow, qhigh = np.loadtxt(quant_file.path)

        with self.output().open('w') as out:

            run_cmd(["bcftools",
                     "filter",
                     "--exclude", "DP<{} | DP>{}".format(int(qlow), int(qhigh)),
                     vcf_input.path], stdout=out)



class PipelineDadi(luigi.WrapperTask):
    """
    Blah...
    """

    # print('\n\n\n-----------------\n\n\n')

    def requires(self):
        yield MergeVCF('horse', 'DOM')


if __name__ == '__main__':
    luigi.run()
