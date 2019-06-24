#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
import glob
import os
import random
import shutil

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.ancient.angsd import CallAncientGenotypes
from alleletraj.const import GROUP_BY_SMPL
from alleletraj.modern.vcf import WholeGenomeSNPsVCF

# PLINK fam file columns (see https://www.cog-genomics.org/plink2/formats#fam)
PLINK_COL_FID = 0
PLINK_COL_IID = 1
PLINK_COL_FATHER = 2
PLINK_COL_MOTHER = 3
PLINK_COL_SEX = 4
PLINK_COL_PHENO = 5

# sex codes for PLINK
PLINK_SEX_UNKNOWN = 0
PLINK_SEX_MALE = 1
PLINK_SEX_FEMALE = 2

# minimum genotyping call rate for samples (--mind) and SNPs (--geno)
PLINK_MIN_MIND = 50
PLINK_MIN_GENO = 90


def plink_sex_code(sex):
    """
    Return the plink sex code for the given string.
    """
    if isinstance(sex, basestring):
        if sex.upper().startswith('M'):
            return PLINK_SEX_MALE
        elif sex.upper().startswith('F'):
            return PLINK_SEX_FEMALE

    return PLINK_SEX_UNKNOWN


class PlinkTask(utils.PipelineTask):
    """
    Luigi class for running plink tasks.

    see https://www.cog-genomics.org/plink/1.9/index
    """

    @property
    def chrset(self):
        """
        Tell plink which chromosomes to expect
        """
        chrset = str(max(map(int, self.autosomes)))
        for chrom in ['X', 'Y', 'MT']:
            if chrom not in self.chromosomes:
                chrset += ' no-{}'.format(chrom.lower())

        return chrset


class PlinkVCFtoBED(PlinkTask, utils.DatabaseTask):
    """
    Convert a VCF to binary PLINK binary format (i.e. bed, bim, fam).

    :type species: str
    """
    species = luigi.Parameter()

    # resources = {'cpu-cores': 1, 'ram-gb': 80}  # TODO what do we actually need here?

    def requires(self):
        return WholeGenomeSNPsVCF(self.species)

    def output(self):
        return [luigi.LocalTarget('data/plink/{}-modern.{}'.format(self.basename, ext)) for ext in
                ['bed', 'bim', 'fam', 'log']]

    def run(self):
        # unpack the params
        vcf_file = self.input()
        bed_file, _, fam_file, _ = self.output()

        # convert GB into MB
        # mem_mb = self.resources['ram-gb'] * 1000

        cmd = ['plink',
               '--chr-set', self.chrset,
               '--make-bed',
               '--double-id',
               '--biallelic-only', 'strict', 'list',
               '--set-missing-var-ids', '@_#',  # NOTE use `CHR_POS` to match how angsd encodes variant IDs
               '--vcf-require-gt',
               '--vcf-half-call', 'missing',
               # '--memory', mem_mb,
               '--vcf', vcf_file.path,
               '--out', utils.trim_ext(bed_file.path)]

        utils.run_cmd(cmd)

        # load the fam file
        with fam_file.open('r') as fin:
            fam_data = fin.readlines()

        # get all the sample records, indexed by sample name
        samples = self.list_samples(modern=True, outgroup=True)
        samples = dict([(sample, samples[(pop, sample)]) for pop, sample in samples])

        # set the correct population and sex codes
        with fam_file.open('w') as fout:
            for line in fam_data:
                fam = line.split()

                # get the sample code
                sample = fam[PLINK_COL_IID]

                # set the population and sex codes
                fam[PLINK_COL_FID] = samples[sample]['population']
                fam[PLINK_COL_SEX] = plink_sex_code(samples[sample]['sex'])

                fout.write(' '.join(map(str, fam)) + '\n')

        # remove the .nosex file
        os.remove('{}.nosex'.format(utils.trim_ext(fam_file.path)))


class HaploToPlink(PlinkTask, utils.DatabaseTask):
    """
    Convert an angsd haplo.gz file to PLINK transposed text format (tped).

    :type species: str
    :type chrom: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        return CallAncientGenotypes(self.species, self.chrom)

    def output(self):
        return [luigi.LocalTarget('data/plink/{}-ancient.{}'.format(self.basename, ext)) for ext in
                ['tped', 'tfam', 'log']]

    def run(self):
        # unpack the params
        hap_file, _, _, _ = self.input()
        tped_file, tfam_file, log_file = self.output()

        log = utils.run_cmd(['haploToPlink', hap_file.path, utils.trim_ext(tfam_file.path)])

        # angst converts all the names to ind[0-9]+, so we need to restore the real population and sample codes
        with tfam_file.open('w') as fout:
            samples = self.list_samples(ancient=True)
            for pop, sample in samples:
                fam = [0 for _ in range(6)]
                fam[PLINK_COL_FID] = pop
                fam[PLINK_COL_IID] = sample
                fam[PLINK_COL_SEX] = plink_sex_code(samples[(pop, sample)]['sex'])
                fam[PLINK_COL_PHENO] = -9

                fout.write(' '.join(map(str, fam)) + '\n')

        with log_file.open('w') as fout:
            fout.write(log)


class PlinkTpedToBed(PlinkTask):
    """
    Convert PLINK tped format to binary (bed), and only keep the sites that intersect with the modern SNPs.

    :type species: str
    :type chrom: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        yield PlinkVCFtoBED(self.species)
        yield HaploToPlink(self.species, self.chrom)

    def output(self):
        return [luigi.LocalTarget('data/plink/{}-ancient.{}'.format(self.basename, ext)) for ext in
                ['bed', 'bim', 'fam', 'log']]

    def run(self):
        # unpack the params
        (_, bim_file, _, _), (tped_file, tfam_file, _) = self.input()
        _, _, _, log_file = self.output()

        # NOTE angst encodes missing genotypes as N when 0 is the default in plink
        cmd = ['plink',
               '--chr-set', self.chrset,
               '--make-bed',
               '--extract', bim_file.path,  # only keep the sites that intersect with the modern SNPs
               '--missing-genotype', 'N',
               '--output-missing-genotype', '0',
               '--tped', tped_file.path,
               '--tfam', tfam_file.path,
               '--out',  utils.trim_ext(log_file.path)]

        utils.run_cmd(cmd)


class PlinkMergeBeds(PlinkTask):
    """
    Merge all the ancient samples into into the modern dataset.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        yield PlinkVCFtoBED(self.species)
        for chrom in self.chromosomes:
            yield PlinkTpedToBed(self.species, chrom)

    def output(self):
        return [luigi.LocalTarget('data/plink/{}-merged.{}'.format(self.basename, ext)) for ext in
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

        merge_list = 'data/plink/{}-merged-{}.list'.format(self.basename, suffix)

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
                 '--out', 'data/plink/{}-merged'.format(self.basename)]

        try:
            # attempt the merge
            utils.run_cmd(merge)

        except Exception as e:
            # TODO rather than dumping all polyallelic sites, we could do the merge iteratively per ancient sample
            missnp_file = 'data/plink/{}-merged-merge.missnp'.format(self.basename)

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
    Produce a list of SNPs in LD with each other to remove from downstream analyses.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return PlinkMergeBeds(self.species)

    def output(self):
        return [luigi.LocalTarget('data/plink/{}-merged-indep.{}'.format(self.basename, ext)) for ext in
                ['prune.in', 'prune.out', 'log']]

    def run(self):
        # unpack the inputs/outputs
        bed_input, _, _, _ = self.input()
        _, _, log_file = self.output()

        # calculate the prune list (prune.in / prune.out)
        cmd = ['plink',
               '--chr-set', self.chrset,
               '--indep-pairwise', 50, 10, 0.3,  # TODO ADMIXTURE manual (chapter 2.3) suggests R^2 of 0.1
               '--bfile', utils.trim_ext(bed_input.path),
               '--out',   utils.trim_ext(log_file.path)]

        utils.run_cmd(cmd)


class PlinkPruneBed(PlinkTask):
    """
    Prune the merged group BED file using the prune.in list from PlinkIndepPairwise.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        yield PlinkMergeBeds(self.species)
        yield PlinkIndepPairwise(self.species)

    def output(self):
        return [luigi.LocalTarget('data/plink/{}-merged-pruned.{}'.format(self.basename, ext)) for ext
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
    Extract a specific population from a larger BED file.

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
                       '--geno',     '0.9999',  # drop sites with no coverage (i.e. less than 0.01% missing)
                       '--keep-fam', pop_list,
                       '--bfile',    utils.trim_ext(bed_input.path),
                       '--out',      utils.trim_ext(bed_output.path)])


class PlinkHighGeno(PlinkTask):
    """
    Drop all sites with a genotyping call rate below the given threshold.

    :type species: str
    :type mind: int
    :type geno: int
    """
    species = luigi.Parameter()
    mind = luigi.IntParameter(default=PLINK_MIN_MIND)
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
        plink_mind = 1 - (self.mind / 100.0)
        plink_geno = 1 - (self.geno / 100.0)

        utils.run_cmd(['plink',
                       '--chr-set', self.chrset,
                       '--make-bed',
                       '--mind',  plink_mind,
                       '--geno',  plink_geno,
                       '--bfile', utils.trim_ext(bed_input.path),
                       '--out',   utils.trim_ext(bed_output.path)])


class PlinkDistMatrix(PlinkTask):
    """
    Compute a distance matrix using IBS.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return PlinkHighGeno(self.species)

    def output(self):
        return [luigi.LocalTarget('data/njtree/{}.{}'.format(self.basename, ext)) for ext in ['data', 'mdist']]

    def run(self):
        # unpack the inputs/outputs
        bed_file, _, fam_file, _ = self.input()
        data_file, mdist_file = self.output()

        # TODO what about bootstrapping?
        # make the distance matrix
        utils.run_cmd(['plink',
                       '--chr-set', self.chrset,
                       '--distance', 'square', '1-ibs',
                       '--bfile',    utils.trim_ext(bed_file.path),
                       '--out',      utils.trim_ext(mdist_file.path)])

        # use awk to extract the sample names
        awk = "awk '{print $2}' " + fam_file.path

        # transpose them into a row
        head = utils.run_cmd([awk + ' | xargs'], shell=True)

        # use awk to extract the sample and population names
        awk = "awk '{print $1\"\\t\"$2}' " + fam_file.path

        # add the samples names as a column to the mdist data
        data = utils.run_cmd([awk + ' | paste - {}'.format(mdist_file.path)], shell=True)

        # save the labeled file
        with data_file.open('w') as fout:
            fout.write('Population\tSample\t' + head)
            fout.write(data)


class PlinkBedToFreq(PlinkTask):
    """
    Convert a BED file into a minor allele frequency report, needed for input into Treemix.

    :type species: str
    :type groupby: str
    """
    species = luigi.Parameter()
    groupby = luigi.Parameter()

    def requires(self):
        return PlinkHighGeno(self.species)

    def output(self):
        return [luigi.LocalTarget('data/plink/{}.{}'.format(self.basename, ext)) for ext in ['frq.strat.gz', 'log']]

    def run(self):
        # unpack the inputs/outputs
        bed_path, bim_path, fam_path, _ = [in_file.path for in_file in self.input()]
        _, log_file = self.output()

        # to group by samples, we need to reassign them to their own families
        if self.groupby == GROUP_BY_SMPL:
            # replace family with sample code
            fam = utils.run_cmd(["awk '{$1=$2}$0' " + fam_path.path], shell=True)

            # make a new fam file
            fam_path = utils.insert_suffix(fam_path, GROUP_BY_SMPL)
            with open(fam_path, 'w') as fout:
                fout.write(fam)

        utils.run_cmd(['plink',
                       '--chr-set', self.chrset,
                       '--freq', 'gz',  # make a gzipped MAF report
                       '--family',      # group by population
                       '--bed', bed_path,
                       '--bim', bim_path,
                       '--fam', fam_path,
                       '--out', utils.trim_ext(log_file.path)])


if __name__ == '__main__':
    luigi.run()
