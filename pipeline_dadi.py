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


class QuantilesOfCoverageVCF(PipelineTask):
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
        # TODO link back to task that made the VCF
        return ModernVCF(self.species, self.population, self.chrom)

    def output(self):
        return luigi.LocalTarget("vcf/{}.quant".format(self.basename))

    def run(self):
        depth = []

        # iterate over the VCF and extract the depth of coverage at each site
        for rec in VariantFile(self.input().path).fetch():
            try:
                depth.append(rec.info['DP'])
            except KeyError:
                pass

        # calculate the quantiles
        quants = np.quantile(depth, [QUANTILE_LOW, QUANTILE_HIGH])

        with self.output().open('w') as fout:
            fout.write('{} {}'.format(int(quants[0]), int(quants[1])))


class FilterQuantilesVCF(PipelineTask):
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
        yield QuantilesOfCoverageVCF(self.species, self.population, self.chrom)

    def output(self):
        return luigi.LocalTarget("vcf/{}.quant.vcf".format(self.basename))

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


class MergeFilteredVCFs(PipelineTask):
    """
    Merge the chromosome level VCFs into a single file.

    :type species: str
    :type population: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()

    def requires(self):
        # TODO replace with CHROM[self.species]
        for chrom in ['1', '2']:
            yield FilterQuantilesVCF(self.species, self.population, 'chr{}'.format(chrom))

    def output(self):
        return luigi.LocalTarget("vcf/{}.chrAll.quant.vcf".format(self.basename))

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

    def requires(self):
        yield MergeFilteredVCFs(self.species, self.population)

    def output(self):
        return luigi.LocalTarget("dadi/{}.sfs".format(self.basename))

    def run(self):
        pass


class PipelineDadi(luigi.WrapperTask):
    """
    Run the dadi pipeline
    """

    # print('\n\n\n-----------------\n\n\n')

    def requires(self):
        yield EasySFS('horse', 'DOM')


if __name__ == '__main__':
    luigi.run()
