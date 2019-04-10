#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import my custom modules
from pipeline_utils import *

# VCF parser
from pysam import VariantFile


class BAMfile(luigi.ExternalTask):
    """
    External task dependency for the aligned BAM file.

    N.B. These have been created outside the workflow of this pipeline.

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
            yield BAMfile(sample)

    def output(self):
        # TODO make vcf.gz and propagate
        return luigi.LocalTarget('vcf/{}.vcf'.format(self.basename))

    def run(self):

        params = {
            'ref': REF_FILE[self.species],
            'chr': self.chrom,
            'bams': ' '.join([bam.path for bam in self.input()]),
            'ploidy':  './data/{}.ploidy'.format(self.species),
            'samples': './data/horses_DOM.samples',
            'vcfout':  './vcf/horse_DOM_chr{}.vcf'
        }

        cmd = "bcftools mpileup --fasta-ref {ref} --regions chr{chr} {bams}" \
              " | bcftools call --multiallelic-caller --ploidy-file {ploidy} --samples-file {samples} --output-type v" \
              " | bcftools view --exclude INFO/INDEL=1 --output-file {vcfout}".format(**params)

        run_cmd([cmd], shell=True)


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
        return BCFToolsCall(self.species, self.population, self.chrom)

    def output(self):
        return luigi.LocalTarget('vcf/{}.quant'.format(self.basename))

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
        return luigi.LocalTarget('vcf/{}.quant.vcf'.format(self.basename))

    def run(self):

        # unpack the input params
        vcf_input, quant_file = self.input()

        # get the quantiles
        qlow, qhigh = np.loadtxt(quant_file.path)

        with self.output().open('w') as vcf_out:

            run_cmd(['bcftools',
                     'filter',
                     '--exclude', 'DP<{} | DP>{}'.format(int(qlow), int(qhigh)),
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
        for chrom in ['10']:
            yield FilterQuantilesVCF(self.species, self.population, 'chr{}'.format(chrom))

    def output(self):
        return luigi.LocalTarget('vcf/{}-chrAll-quant.vcf'.format(self.basename))

    def run(self):

        # unpack the input params
        vcf_files = [vcf.path for vcf in self.input()]

        with self.output().open('w') as vcf_out:
            run_cmd(['bcftools', 'concat'] + vcf_files, stdout=vcf_out)


class PolarizeVCF(PipelineTask):
    """
    Switch the REF/ALT to match the ancestral/derived allele, based on an outgroup in the VCF.

    :type species: str
    :type population: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()

    def requires(self):
        return MergeFilteredVCFs(self.species, self.population)

    def output(self):
        return luigi.LocalTarget('vcf/{}-chrAll-quant-polar.vcf'.format(self.basename))

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
                alt = set((rec.ref,) + rec.alts)
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


class PipelineDadi(luigi.WrapperTask):
    """
    Run the dadi pipeline
    """

    # print('\n\n\n-----------------\n\n\n')

    def requires(self):
        yield PolarizeVCF('horse', 'DOM2')


if __name__ == '__main__':
    luigi.run()
