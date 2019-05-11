#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
import glob
import os
import random
import shutil

# third party modules
import struct

import luigi

# local modules
from alleletraj import utils
from alleletraj.ancient.call import CallAncientGenotypes
from alleletraj.modern.vcf import WholeGenomeSNPsVCF

# minimum genotype call rate
PLINK_MIN_GENO = 0.9


def plink_extract_sps(input_prefix, bim_file, output_prefix):
    """
    Extracts the positions given by the bim file.

    :param input_prefix:
    :param bim_file:
    :param output_prefix:
    :return:
    """

    # copy the FAM file
    shutil.copyfile(input_prefix + '.fam', output_prefix + '.fam')

    # extract the list of sites from the BIM file
    sites = set()

    # extract the list of sites from the SNP array using the supplied BIM file
    with open(bim_file, 'r') as infile:
        for line in infile:
            # get the chromosome and position of each site
            chrom, _, _, pos, _, _ = line.split()

            # add the locus to the set
            sites.add((chrom, pos))

    # open the output files for writing
    with open(output_prefix + '.bed', 'wb') as bed_fout, open(output_prefix + '.bim', 'w') as bim_fout:

        # initialise the new BED file (see https://www.cog-genomics.org/plink2/formats#bed)
        bed_fout.write(struct.pack('b', 0x6c))
        bed_fout.write(struct.pack('b', 0x1b))
        bed_fout.write(struct.pack('b', 0x01))

        # and the input files for reading
        with open(input_prefix + '.bim', 'r') as bim_fin:

            for line in bim_fin:
                # get the chromosome and position of each site
                chrom, _, _, pos, _, _ = line.split()

                # skip any sites not found in the SNP array
                if (chrom, pos) not in sites:
                    continue

                # otherwise add this locus to the output files
                bim_fout.write(line)

                # TODO read the input BED rather than assume a genotype of b11
                bed_fout.write(struct.pack('b', 3))


class PlinkTask(utils.PipelineTask):
    """
    Luigi class for running plink tasks

    see https://www.cog-genomics.org/plink/1.9/index
    """

    @property
    def chrset(self):
        """
        Tell plink know which chromosomes to expect.
        """
        chrset = max([chrom for chrom in self.chromosomes if chrom.isdigit()])
        for chrom in ['X', 'Y', 'MT']:
            if chrom not in self.chromosomes:
                chrset += ' no-{}'.format(chrom.lower())

        return chrset


class PlinkVCFtoBED(PlinkTask):
    """
    Convert a VCF to binary plink format (i.e. BED, BIM, FAM)

    :type species: str
    """
    species = luigi.Parameter()

    resources = {'cpu-cores': 1, 'ram-gb': 80}

    def requires(self):
        return WholeGenomeSNPsVCF(self.species)

    def output(self):
        return [luigi.LocalTarget('data/plink/{}.{}'.format(self.basename, ext)) for ext in
                ['bed', 'bim', 'fam', 'log']]

    def run(self):
        # unpack the params
        vcf_file = self.input()
        bed_file, _, _, _ = self.output()

        # convert GB into MB
        mem_mb = self.resources['ram-gb'] * 1000

        cmd = ['plink',
               '--chr-set', self.chrset,
               '--make-bed',
               '--double-id',
               '--biallelic-only', 'strict', 'list',
               '--set-missing-var-ids', '@-#',
               '--vcf-require-gt',
               '--vcf-half-call', 'missing',
               '--memory', mem_mb,
               '--vcf', vcf_file.path,
               '--out', utils.trim_ext(bed_file.path)]

        utils.run_cmd(cmd)

        # TODO set sex of modern samples (--update-sex) https://www.cog-genomics.org/plink/1.9/data#update_indiv


class PlinkExtractSNPs(PlinkTask):
    """
    Extract only those SNPs which appear in the modern dataset.

    :type species: str
    :type population: str
    :type sample: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    sample = luigi.Parameter()

    def requires(self):
        yield PlinkVCFtoBED(self.species)
        yield CallAncientGenotypes(self.population, self.sample)

    def output(self):
        return [luigi.LocalTarget('data/plink/{}.{}'.format(self.basename, ext)) for ext in
                ['bed', 'bim', 'fam', 'log']]

    def run(self):
        # unpack the inputs/outputs
        (_, ref_bim, _, _), (bed_input, _, _, _) = self.input()
        bed_output, _, _, log_file = self.output()

        # random call the ancient genotypes
        plink_extract_sps(utils.trim_ext(bed_input.path), ref_bim.path, utils.trim_ext(bed_output.path))

        # write a log file to show this task finished
        with log_file.open('w') as log_fout:
            log_fout.write('Done!')


class PlinkMergeBeds(PlinkTask):
    """
    Merge all the ancient samples into into the modern dataset

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        yield PlinkVCFtoBED(self.species)

        for population, sample in self.all_samples:
            yield PlinkExtractSNPs(self.species, population, sample)

    def output(self):
        return [luigi.LocalTarget('data/plink/{}.merged.{}'.format(self.basename, ext)) for ext in
                ['bed', 'bim', 'fam', 'log']]

    def run(self):
        bed_files = []

        # generate a unique suffix for temporary files
        suffix = 'luigi-tmp-{:010}'.format(random.randrange(0, 1e10))

        # compose the list of files to merge
        for inputs in self.input():
            bed_files.append(utils.trim_ext(inputs[0].path) + '.' + suffix)

            for plink_file in inputs:
                # make a copy of the BED/BIM/FAM files because we'll need to filter them
                shutil.copyfile(plink_file.path, utils.insert_suffix(plink_file.path, suffix))

        merge_list = 'data/plink/{}.merged.list'.format(self.basename)

        # plink requires the name of the first bed file to be given in the command
        bed_first = bed_files[0]

        # and the rest of the bed files to be listed in a merge file
        with open(merge_list, 'w') as fout:
            fout.write('\n'.join(bed_files[1:]))

        # compose the merge command, because we are going to need it twice
        merge = ['plink',
                 '--chr-set', self.chrset,
                 '--make-bed',
                 '--bfile', bed_first,
                 '--merge-list', merge_list,
                 '--out', 'data/plink/{}.merged'.format(self.basename)]

        try:
            # attempt the merge
            utils.run_cmd(merge)

        except Exception as e:
            missnp_file = 'data/plink/{}.merged-merge.missnp'.format(self.basename)

            # handle merge errors
            if os.path.isfile(missnp_file):

                # filter all the BED files, using the missnp file created by the failed merge
                for bed_file in bed_files:
                    utils.run_cmd(['plink',
                                   '--chr-set', self.chrset,
                                   '--make-bed',
                                   '--exclude', missnp_file,
                                   '--bfile', bed_file,
                                   '--out', bed_file])

                # reattempt the merge
                utils.run_cmd(merge)

            else:
                raise Exception(e)

        # tidy up all the temporary files
        for tmp in glob.glob('data/plink/*{}*'.format(suffix)):
            os.remove(tmp)


class PlinkIndepPairwise(PlinkTask):
    """
    Produce a list of SNPs with high discriminating power.

    Filter out sites with low minor allele frequency and high linkage disequilibrium.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return PlinkMergeBeds(self.species)

    def output(self):
        return [luigi.LocalTarget('data/plink/{}.prune.{}'.format(self.basename, ext)) for ext in ['in', 'out', 'log']]

    def run(self):
        # unpack the inputs/outputs
        bed_input, _, _, _ = self.input()
        _, _, log_file = self.output()

        # calculate the prune list (prune.in / prune.out)
        cmd = ['plink',
               '--chr-set', self.chrset,
               '--indep-pairwise', 50, 10, 0.5,  # accept R^2 coefficient of up to 0.5
               '--bfile', utils.trim_ext(bed_input.path),
               '--out',   utils.trim_ext(bed_input.path)]

        log = utils.run_cmd(cmd)

        # write the log file
        with log_file.open('w') as fout:
            fout.write(log)


class PlinkPruneBed(PlinkTask):
    """
    Prune the merged group BED file using the prune.in list from PlinkIndepPairwise()

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        yield PlinkMergeBeds(self.species)
        yield PlinkIndepPairwise(self.species)

    def output(self):
        return [luigi.LocalTarget('data/plink/{}.pruned.{}'.format(self.basename, ext)) for ext
                in ['bed', 'bim', 'fam', 'log']]

    def run(self):
        # unpack the inputs/outputs
        (bed_input, _, _, _), (prune_in, _, _) = self.input()
        bed_output, _, _, _ = self.output()

        # apply the prune list
        utils.run_cmd(['plink',
                       '--chr-set', self.chrset,
                       '--make-bed',
                       '--extract', prune_in.path,
                       '--bfile',   utils.trim_ext(bed_input.path),
                       '--out',     utils.trim_ext(bed_output.path)])


class PlinkExtractPop(PlinkTask):
    """
    Extract a single population from a larger BED file

    :type species: str
    :type population: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()

    def requires(self):
        return PlinkMergeBeds(self.species)

    def output(self):
        return [luigi.LocalTarget('data/plink/{}.{}'.format(self.basename, ext)) for ext in
                ['bed', 'bim', 'fam', 'log']]

    def run(self):
        # unpack the inputs/outputs
        bed_input, _, _, _ = self.input()
        bed_output, _, _, _ = self.output()

        pop_list = 'data/plink/{}.poplist'.format(self.basename)

        # write the population IDs to disk so we can filter out any unwanted "families"
        with open(pop_list, 'w') as fout:
            fout.write(self.population)

        # apply the prune list
        utils.run_cmd(['plink',
                       '--chr-set', self.chrset,
                       '--make-bed',
                       '--geno',     '0.999',  # drop sites with no coverage, or rather, where less than 0.1% is missing
                       '--keep-fam', pop_list,
                       '--bfile',    utils.trim_ext(bed_input.path),
                       '--out',      utils.trim_ext(bed_output.path)])


class PlinkHighGeno(PlinkTask):
    """
    Drop all sites with a genotyping call rate below the given threshold.

    :type species: str
    :type geno: int
    """
    species = luigi.Parameter()
    geno = luigi.IntParameter(default=PLINK_MIN_GENO)

    def requires(self):
        return PlinkMergeBeds(self.species)

    def output(self):
        return [luigi.LocalTarget('data/plink/{}.{}'.format(self.basename, ext)) for ext in
                ['bed', 'bim', 'fam', 'log']]

    def run(self):
        # unpack the inputs/outputs
        bed_input, _, _, _ = self.input()
        bed_output, _, _, _ = self.output()

        # plink requires the genotyping rate to be expressed as the missing threshold (i.e. 90% = 0.1)
        plink_geno = 1 - self.geno

        utils.run_cmd(['plink',
                       '--chr-set', self.chrset,
                       '--make-bed',
                       '--geno',  plink_geno,
                       '--bfile', utils.trim_ext(bed_input.path),
                       '--out',   utils.trim_ext(bed_output.path)])


if __name__ == '__main__':
    luigi.run()
