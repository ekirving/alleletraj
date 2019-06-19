#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.const import CPU_CORES_MED
from alleletraj.modern.trim import TrimGalore
from alleletraj.ref import ReferenceFASTA, BwaIndexBWTSW


class BwaMem(utils.PipelineTask):
    """
    Align the fastq file(s) to the reference genome using the BWA-MEM algorithm.

    http://www.htslib.org/workflow/#mapping_to_variant

    :type species: str
    :type sample: str
    :type accession: str
    :type paired: bool
    """
    species = luigi.Parameter()
    sample = luigi.Parameter()
    accession = luigi.Parameter()
    paired = luigi.BoolParameter()

    resources = {'cpu-cores': CPU_CORES_MED, 'ram-gb': 16}

    def requires(self):
        yield TrimGalore(self.accession, self.paired)
        yield ReferenceFASTA(self.species)
        yield BwaIndexBWTSW(self.species)

    def output(self):
        return [luigi.LocalTarget('data/bam/{}.sort.{}'.format(self.accession, ext)) for ext in
                ['bam', 'bam.bai', 'log']]

    def run(self):
        # unpack the input params
        fastq_input, (ref_file, _), _ = self.input()
        bam_out, _, log_file = self.output()

        params = {
            # CPU threads to use
            'threads': self.resources['cpu-cores'],

            # compose the read group metadata
            'readgroup': r'@RG\tID:{sample}\tSM:{sample}'.format(sample=self.sample),

            # reference genome
            'reference': ref_file.path,

            # get the fastq file(s) to align
            'fastq': ' '.join(fastq.path for fastq in fastq_input),
        }

        with bam_out.temporary_path() as bam_path, log_file.open('w') as fout:
            # get the temporary path for the bam file
            params['bam'] = bam_path

            # align using bwa-mem
            cmd = "bwa mem -t {threads} -R '{readgroup}' {reference} {fastq} | "

            if self.paired:
                # fix paired end flags
                cmd += "samtools fixmate - - | "

            # drop unaligned reads and convert to a sorted BAM file
            cmd += "samtools view -h -F 0x0004 | " \
                   "samtools sort -@ {threads} -O bam -o {bam}".format(**params)

            utils.run_cmd([cmd], shell=True, stderr=fout)

        # index the BAM file
        utils.run_cmd(['samtools', 'index', '-b', bam_out.path])


if __name__ == '__main__':
    luigi.run()
