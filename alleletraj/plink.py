#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
import os

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.modern.vcf import BiallelicSNPsVCF

# the species flag for plink telling it how many chromosomes to expect
PLINK_TAXA = '--dog'

# column delimiter for PED files
PLINK_COL_DELIM = ' '

# status codes for PED files
PLINK_UNKNOWN = 0
PLINK_MISSING_PHENO = -9
PLINK_MISSING_GENO = 0

# min genotyping rate
PLINK_MIN_GENO = 90


class PlinkVCFtoBED(utils.PipelineTask):
    """
    Convert a VCF to binary plink format (i.e. BED, BIM, FAM)

    :type species: str
    :type chrom: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    resources = {'cpu-cores': 1, 'ram-gb': 80}  # TODO make this smaller now we're splitting on chrom

    def requires(self):
        return BiallelicSNPsVCF(self.species, self.chrom)

    def output(self):
        return [luigi.LocalTarget('data/plink/{}.{}'.format(self.basename, ext)) for ext in
                ['bed', 'bim', 'fam', 'log']]

    def run(self):
        # unpack the inputs/outputs
        vcf_file = self.input()
        bed_file, bim_file, _, log_file = self.output()

        # convert GB into MB
        mem_mb = self.resources['ram-gb'] * 1000

        # setup the base query
        cmd = ['plink',
               PLINK_TAXA,
               '--make-bed',
               '--double-id',
               '--allow-extra-chr',
               '--biallelic-only', 'strict', 'list',
               '--set-missing-var-ids', '@-#',
               '--vcf-require-gt',
               '--vcf-half-call', 'missing',
               '--memory', mem_mb,
               '--vcf', vcf_file.path,
               '--out', utils.trim_ext(bed_file.path)]

        if self.dataset == 'NHGRI_722g_hq':
            # drop all filtered sites
            cmd.append('--vcf-filter')

        # convert the VCF
        utils.run_cmd(cmd)

        # if the VCF had named variants (e.g. rsnumbers) then they will clash with our chrom-pos notation
        utils.run_cmd(['mv {bim} {bim}.bak'.format(bim=bim_file.path)], shell=True)
        utils.run_cmd(["awk '{{$2=$1\"-\"$4; print $0}}' {bim}.bak > {bim}".format(bim=bim_file.path)], shell=True)
        os.remove('{bim}.bak'.format(bim=bim_file.path))


class PlinkExtractSNPs(utils.PipelineTask):
    """
    Extract only those SNPs which appear in the reference dataset.

    :type dataset: str
    :type population: str
    :type sample: str
    """
    dataset = luigi.Parameter()
    population = luigi.Parameter()
    sample = luigi.Parameter()

    def requires(self):
        yield PlinkVCFtoBED(self.dataset)
        yield CallAncientGenotypes(self.population, self.sample)

    def output(self):
        return [luigi.LocalTarget('plink/{}.{}'.format(self.basename, ext)) for ext in ['bed', 'bim', 'fam', 'log']]

    def run(self):

        # unpack the inputs/outputs
        (_, ref_bim, _, _), (bed_input, _, _, _) = self.input()
        bed_output, _, _, log_file = self.output()

        # random call the ancient genotypes
        plink_extract_sps(utils.trim_ext(bed_input.path), ref_bim.path, utils.trim_ext(bed_output.path))

        # write a log file to show this task finished
        with log_file.open('w') as log_fout:
            log_fout.write('Done!')


class PlinkMergeBeds(utils.PipelineTask):
    """
    Merge multiple BED files into one

    :type dataset: str
    """
    dataset = luigi.Parameter()

    def requires(self):

        if self.dataset in ['NHGRI_722g', 'NHGRI_722g_hq']:
            # convert the VCF of modern data
            yield PlinkVCFtoBED(self.dataset)
        else:
            # get the modern reference data
            yield ModernReferencePanel(self.dataset)

        # and merge it with all the ancient samples
        for population in ANCIENT_POPS:
            for sample in ANCIENT_POPS[population]:
                yield PlinkExtractSNPs(self.dataset, population, sample)

    def output(self):
        return [luigi.LocalTarget('plink/{0}.merged.{1}'.format(self.dataset, ext)) for ext
                in ['bed', 'bim', 'fam', 'log']]

    def run(self):

        # generate a unique suffix for temporary files
        suffix = 'tmp' + str(random.getrandbits(100))

        bed_files = []

        # compose the list of files to merge
        for inputs in self.input():
            bed_files.append(utils.trim_ext(inputs[0].path) + '.' + suffix)

            for plink_file in inputs:
                # make a copy of the BED/BIM/FAM files because we'll need to filter them
                shutil.copyfile(plink_file.path, insert_suffix(plink_file.path, suffix))

        # use the first input as the named BED file for the merge command
        bed_snparray = bed_files[0]

        merge_list = 'plink/{}.merged.list'.format(self.dataset)

        # save the merge-list
        with open(merge_list, 'w') as fout:
            fout.write('\n'.join(bed_files[1:]))

        # compose the merge command, because we are going to need it twice
        merge = ['plink',
                 PLINK_TAXA,
                 '--make-bed',
                 '--allow-extra-chr',  # prevent errors with X, MT and unplaced scaffolds
                 '--bfile', bed_snparray,
                 '--merge-list', merge_list,
                 '--out', 'plink/{}.merged'.format(self.dataset)]

        try:
            # attempt the merge
            utils.run_cmd(merge)

        except Exception as e:
            missnp_file = 'plink/{}.merged-merge.missnp'.format(self.dataset)

            # handle merge errors
            if os.path.isfile(missnp_file):

                # filter all the BED files, using the missnp file created by the failed merge
                for bed_file in bed_files:
                    utils.run_cmd(['plink',
                                   PLINK_TAXA,
                                   '--make-bed',
                                   '--allow-extra-chr',  # prevent errors with X, MT and unplaced scaffolds
                                   '--exclude', missnp_file,
                                   '--bfile', bed_file,
                                   '--out', bed_file])

                # reattempt the merge
                utils.run_cmd(merge)

            else:
                raise Exception(e)

        # tidy up all the temporary intermediate files
        for tmp in glob.glob('./plink/*{}*'.format(suffix)):
            os.remove(tmp)
        for tmp in glob.glob('./data/*{}*'.format(suffix)):
            os.remove(tmp)


class MesoMergeBeds(luigi.WrapperTask):
    """
    Merge all the meso-dogs datasets
    """

    def requires(self):
        for dataset in DATASETS:
            yield PlinkMergeBeds(dataset)

            # TODO rename populations
            # TODO merge sex


if __name__ == '__main__':
    luigi.run()
