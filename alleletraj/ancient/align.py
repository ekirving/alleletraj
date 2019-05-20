#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.ancient.trim import AdapterRemoval
from alleletraj.const import CPU_CORES_LOW, CPU_CORES_MED
from alleletraj.ref import ReferenceFASTA, BwaIndexBWTSW

# BWA-ALN settings
BWA_ALN_SEED = 1024  # long seed effectively disables seeding
BWA_ALN_EDIT = 0.03  # small edit distance improves mapping when reads have errors


class BwaAln(utils.PipelineTask):
    """
    Align the fastq file to the reference genome using the BWA-ALN algorithm.

    Step 1: Find the SA coordinates of the input reads.

    Sets aDNA specific thresholds in the alignment to improve mapping sensitivity.

    https://www.mdpi.com/2073-4425/9/3/157

    :type species: str
    :type sample: str
    :type accession: str
    :type paired: bool
    """
    species = luigi.Parameter()
    sample = luigi.Parameter()
    accession = luigi.Parameter()
    paired = luigi.BoolParameter()

    resources = {'cpu-cores': CPU_CORES_MED, 'ram-gb': 8}  # TODO check ram usage

    def requires(self):
        yield AdapterRemoval(self.accession, self.paired)
        yield ReferenceFASTA(self.species)
        yield BwaIndexBWTSW(self.species)

    def output(self):
        return luigi.LocalTarget('data/bam/{}.sai'.format(self.accession))

    def run(self):
        # unpack the params
        (fastq_file, _, _), (ref_file, _), _ = self.input()
        sai_out = self.output()

        params = {
            'seed': BWA_ALN_SEED,  # disable seeding by setting very long length
            'edit': BWA_ALN_EDIT,  # maximum edit distance
            'threads': self.resources['cpu-cores'],
            'reference': ref_file.path,
            'fastq': fastq_file.path,
        }

        with sai_out.temporary_path() as sai_path:
            # get the temporary path for the sai file
            params['sai'] = sai_path

            # align using the ALN algorithm
            cmd = "bwa aln -l {seed} -n {edit} -t {threads} {reference} {fastq} > {sai}".format(**params)

            # perform the alignment
            utils.run_cmd([cmd], shell=True)


class BwaSamSe(utils.PipelineTask):
    """
    Align the fastq file to the reference genome using the BWA-ALN algorithm.

    Step 2: Generate alignments in the SAM format, and convert to sorted BAM.

    Becase we collapse paired-end reads we use SAM-SE, not SAM-PE, to produce the alignment.

    :type species: str
    :type sample: str
    :type accession: str
    :type paired: bool
    """
    species = luigi.Parameter()
    sample = luigi.Parameter()
    accession = luigi.Parameter()
    paired = luigi.BoolParameter()

    resources = {'cpu-cores': CPU_CORES_LOW, 'ram-gb': 8}  # TODO check ram usage

    def requires(self):
        yield AdapterRemoval(self.accession, self.paired)
        yield BwaAln(self.species, self.sample, self.accession, self.paired)
        yield ReferenceFASTA(self.species)

    def output(self):
        return [luigi.LocalTarget('data/bam/{}.sort.{}'.format(self.accession, ext)) for ext in
                ['bam', 'bam.bai','log']]

    def run(self):
        # unpack the params
        (fastq_file, _, _), sai_file, (ref_file, _) = self.input()
        bam_out, _, log_file = self.output()

        params = {
            'readgroup': r'@RG\tID:{sample}\tSM:{sample}'.format(sample=self.sample),
            'reference': ref_file.path,
            'sai':       sai_file.path,
            'fastq':     fastq_file.path,
            'threads':   self.resources['cpu-cores'],
        }

        with bam_out.temporary_path() as bam_path, log_file.open('w') as fout:
            # get the temporary path for the bam file
            params['bam'] = bam_path

            # align using the ALN algorithm, and convert to a sorted BAM file
            cmd = "bwa samse -r '{readgroup}' {reference} {sai} {fastq} " \
                  " | samtools sort -@ {threads} -O bam -o {bam} -".format(**params)

            # perform the alignment
            utils.run_cmd([cmd], shell=True, stderr=fout)

        # index the BAM file
        utils.run_cmd(['samtools', 'index', '-b', bam_out.path])


if __name__ == '__main__':
    luigi.run()
