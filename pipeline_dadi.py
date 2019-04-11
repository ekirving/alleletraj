#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import my custom modules
from pipeline_utils import *

from pipeline_snp_call import PolarizeVCF, SubsetSNPsVCF


class CountCallableSites(PipelineTask):
    """
    Count the number of callable sites, as dadi needs this number to estimate the ancestral population size from theta.

    :type species: str
    :type population: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()

    def requires(self):
        return PolarizeVCF(self.species, self.population)

    def output(self):
        return luigi.LocalTarget('sfs/{}.L'.format(self.basename))

    def run(self):

        # count all unique sites
        size = run_cmd(["bcftools query -f '%CHROM %POS\\n' {} | uniq | wc -l".format(self.input().path)], shell=True)

        with self.output().open('w') as fout:
            fout.write(size)


class EasySFS(PipelineTask):
    """
    Calculate the Site Frequency Spectrum.

    :type species: str
    :type population: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()

    def requires(self):
        return SubsetSNPsVCF(self.species, self.population)

    def output(self):
        return luigi.LocalTarget('sfs/{}/dadi/{}.sfs'.format(self.basename, self.population))

    def run(self):

        # make a sample/population file
        pop_file = 'sfs/{}.pops'.format(self.basename)
        num_samples = 0

        with open(pop_file, 'w') as fout:
            for sample in SAMPLES[self.species][self.population]:
                if sample not in SFS_EXCLUSIONS[self.species]:
                    fout.write('{}\t{}\n'.format(sample, self.population))
                    num_samples += 1

        params = {
            'vcf':  self.input().path,
            'pops': pop_file,
            'out':  self.basename,
            'proj': num_samples * 2  # don't project down
        }

        # pass 'yes' into easySFS to get past the interactive prompt
        cmd = "echo 'yes' | easySFS.py -a -i {vcf} -p {pops} -o sfs/{out} --proj {proj} --unfolded".format(**params)

        run_cmd([cmd], shell=True)


class DadiModelDemography(luigi.WrapperTask):
    """
    Find the best fitting of 5 sequential epoch models (i.e. 1 epoch, 2 epoch, etc.).
    """

    def requires(self):
        return EasySFS('horse', 'DOM2')


if __name__ == '__main__':
    luigi.run()
