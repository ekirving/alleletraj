#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.bam import SampleBAM
from alleletraj.modern.vcf import ReferencePloidy, WholeGenomeSNPsVCF, MIN_GENO_QUAL
from alleletraj.ref import ReferenceFASTA

# minimum depth of coverage to call diploid genotypes
MIN_GENO_DEPTH = 10


class BCFToolsTargetsFile(utils.DatabaseTask):
    """
    Make a targets file to contstrain the callable sites and alleles in the ancint samples.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        yield ReferenceFASTA(self.species)
        yield WholeGenomeSNPsVCF(self.species)

    def output(self):
        return [luigi.LocalTarget('data/vcf/{}-chrAll-quant-polar-SNPs.{}'.format(self.species, ext))
                for ext in ['tsv.gz', 'tsv.gz.tbi']]

    def run(self):
        # unpack the params
        (ref_file, _), vcf_file = self.input()
        tsv_file, _ = self.output()

        params = {
            'ref': ref_file.path,
            'vcf': vcf_file.path,
            'tsv': tsv_file.path
        }

        # use `bcftools norm` to unpolarize the SNP VCF
        cmd = "bcftools norm --check-ref s --fasta-ref {ref} --do-not-normalize {vcf} | " \
              "bcftools query -f'%CHROM\\t%POS\\t%REF,%ALT\\n' | " \
              "bgzip -c > {tsv} && tabix -s1 -b2 -e2 {tsv}".format(**params)

        utils.run_cmd([cmd], shell=True)


class BCFToolsCallAncient(utils.DatabaseTask):
    """
    Make genotype calls using the bcftools mpileup workflow, but contstrain the callable sites and alleles.

    Also, normalise indels and merge multiallelic sites, then recalculate AN and AC tags as `norm` doesn't do this.

    :type species: str
    :type population: str
    :type sample: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    sample = luigi.Parameter()

    def requires(self):
        yield ReferenceFASTA(self.species)
        yield ReferencePloidy(self.species)
        yield BCFToolsTargetsFile(self.species)
        yield SampleBAM(self.species, self.population, self.sample)

    def output(self):
        return [luigi.LocalTarget('data/vcf/{}.vcf.{}'.format(self.basename, ext)) for ext in ['gz', 'log']]

    def run(self):
        # unpack the params
        (ref_file, _), pld_file, (tsv_file, _), (bam_file, _) = self.input()
        vcf_file, log_file = self.output()

        sex_file = luigi.LocalTarget(is_tmp=True)
        rgs_file = luigi.LocalTarget(is_tmp=True)

        # bcftools needs the sex specified in a separate file
        with sex_file.open('w') as fout:
            fout.write('{}\t{}\n'.format(self.sample, self.sample_data['sex']))

        # sample names in the BAM file(s) may not be consistent, so override the @SM code
        with rgs_file.open('w') as fout:
            fout.write('*\t{}\t{}\n'.format(bam_file.path, self.sample))

        with vcf_file.temporary_path() as vcf_path, log_file.open('w') as fout:

            params = {
                'ref': ref_file.path,
                'tsv': tsv_file.path,
                'rgs': rgs_file.path,
                'bam': bam_file.path,
                'pld': pld_file.path,
                'sex': sex_file.path,
                'vcf': vcf_path
            }

            cmd = "bcftools mpileup --fasta-ref {ref} --targets-file {tsv} --read-groups {rgs} -O u {bam} | " \
                  "bcftools call --multiallelic-caller --ploidy-file {pld} --samples-file {sex} --targets-file {tsv} " \
                  " --constrain alleles --output-type u | " \
                  "bcftools norm --fasta-ref {ref} --multiallelics +any --output-type u | " \
                  "bcftools +fill-tags --output-type z --output {vcf} -- -t AN,AC ".format(**params)

            utils.run_cmd([cmd], shell=True, stderr=fout)


class BCFToolsFilterAncientVCF(utils.PipelineTask):
    """
    Filter out sites with low quality or low depth, and all INDELS and other junk.

    Retains only high quality biallelic SNPs and homozygous reference calls.

    :type species: str
    :type population: str
    :type sample: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    sample = luigi.Parameter()

    def requires(self):
        yield ReferenceFASTA(self.species)
        yield BCFToolsCallAncient(self.species, self.population, self.sample)

    def output(self):
        return [luigi.LocalTarget('data/vcf/{}-filter.vcf.{}'.format(self.basename, ext)) for ext in ['gz', 'gz.tbi']]

    def run(self):
        # unpack the input params
        (ref_file, _), (vcf_input, _) = self.input()
        vcf_file, _ = self.output()

        with vcf_file.temporary_path() as vcf_out:

            params = {
                'qual': MIN_GENO_QUAL,
                'dp':   MIN_GENO_DEPTH,
                'vcf':  vcf_input.path,
                'ref':  ref_file.path,
                'out':  vcf_out,
            }

            cmd = "bcftools filter --exclude 'QUAL<{qual} | DP<{dp}' --output-type u {vcf} | " \
                  "bcftools view --exclude-types indels,mnps,bnd,other -O z --output-file {out}".format(**params)

            utils.run_cmd([cmd], shell=True)

        # index the vcf
        utils.run_cmd(['bcftools', 'index', '--tbi', vcf_file.path])
