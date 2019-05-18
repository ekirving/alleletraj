#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.const import GROUP_BY_POPS
from alleletraj.plink import PlinkMergeBeds


class ConvertfBedToEigenstrat(utils.PipelineTask):
    """
    Convert a BED file into Eigenstrat format, for use by admixtools.

    :type species: str
    :type groupby: str
    """
    species = luigi.Parameter()
    groupby = luigi.Parameter(default=GROUP_BY_POPS)

    def requires(self):
        return PlinkMergeBeds(self.species)

    def output(self):
        return [luigi.LocalTarget('data/eigenstrat/{}.{}'.format(self.basename, ext)) for ext
                in ['geno', 'snp', 'ind', 'par', 'log']]

    def run(self):
        # unpack the inputs/outputs
        bed_input, bim_input, fam_input, _ = self.input()
        geno_file, snp_file, ind_file, par_file, log_file = self.output()

        # N.B. admixtools requires an invalid version of the .fam file, where the phenotype column has been replaced
        # with the population name. If you don't do this then none of the populations in poplistname will be found.
        # Also, if you happen to be using qtmode:YES (because otherwise a value of -9 throws errors!) then it will fail
        # to find the populations, even if you've make the invalid version of the .fam file that poplistname wants.
        if self.groupby == GROUP_BY_POPS:
            fam = utils.run_cmd(["awk '$6=$1' " + fam_input.path], shell=True)

        else:  # self.groupby == GROUP_BY_SMPL
            # also... we want to pretend that all our samples are populations!
            fam = utils.run_cmd(["awk '$6=$2' " + fam_input.path], shell=True)

        # save the fam file
        fam_file = utils.insert_suffix(fam_input.path, self.groupby)
        with open(fam_file, 'w') as fout:
            fout.write(fam)

        # compose the config settings for convertf
        config = [
            'genotypename:    {}'.format(bed_input.path),
            'snpname:         {}'.format(bim_input.path),
            'indivname:       {}'.format(fam_file),
            'outputformat:    EIGENSTRAT',
            'genotypeoutname: {}'.format(geno_file.path),
            'snpoutname:      {}'.format(snp_file.path),
            'indivoutname:    {}'.format(ind_file.path),
            'familynames:     NO',
            'pordercheck:     NO'
        ]

        # convertf needs the params to be defined in a .par file
        with par_file.open('w') as fout:
            fout.write('\n'.join(config))

        # run convertf
        log = utils.run_cmd(['convertf', '-p', par_file.path])

        # save the log file
        with log_file.open('w') as fout:
            fout.write(log)
