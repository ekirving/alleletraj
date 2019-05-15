#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.const import CPU_CORES_LOW, CPU_CORES_MED
from alleletraj.ref import ReferenceFASTA, BwaIndexBWTSW
from alleletraj.sra import SraToolsFastqDump

# hard filters for AdapterRemoval
TRIM_MIN_BASEQ = 20
TRIM_MIN_LENGTH = 25

# BWA-ALN settings
BWA_ALN_SEED = 1024  # long seed effectively disables seeding
BWA_ALN_EDIT = 0.03  # small edit distance improves mapping when reads have errors


class AdapterRemoval(utils.PipelineTask):
    """
    Trim adapters, low-quality bases from the 5' and 3' ends, drop any resulting reads below a min-length threshold,
    and collapse overlapping read pairs.

    Becase we collapse paired-end reads, this task always produces a single fastq file.

    :type accession: str
    :type paired: bool
    """
    accession = luigi.Parameter()
    paired = luigi.BoolParameter()

    def requires(self):
        return SraToolsFastqDump(self.accession, self.paired)

    def output(self):
        yield luigi.LocalTarget('data/fastq/{}-trim.fq.gz'.format(self.accession))
        yield luigi.LocalTarget('data/fastq/{}-trim-discard.fq.gz'.format(self.accession))
        yield luigi.LocalTarget('data/fastq/{}-trim.log'.format(self.accession))

    def run(self):
        fastq_files = self.input()
        fastq_keep, fastq_drop, log_file = self.output()

        cmd = [
            'AdapterRemoval',
            '--discarded',  fastq_drop.path,
            '--settings',   log_file.path,
            '--trimns',                       # trim ambiguous bases (N) at both 5'/3' termini
            '--trimqualities',                # trim bases at both 5'/3' termini with quality scores...
            '--minquality', TRIM_MIN_BASEQ,   # below this value
            '--minlength',  TRIM_MIN_LENGTH,  # discard reads shorter than this length following trimming
            '--gzip'
        ]

        if self.paired:
            cmd += [
                '--file1', fastq_files[0].path,
                '--file2', fastq_files[1].path,
                '--collapse',
                '--outputcollapsed', fastq_keep.path
            ]
        else:
            # TODO still need to test paired end mode works correctly
            cmd += [
                '--basename', utils.trim_ext(log_file.path),
                '--file1',    fastq_files[0].path,
                '--output1',  fastq_keep.path
            ]

        # perform the trimming
        utils.run_cmd(cmd)


class BwaAln(utils.PipelineTask):
    """
    Align the fastq file to the reference genome using the BWA-ALN algorithm.

    Step 1: Find the SA coordinates of the input reads.

            Sets aDNA specific thresholds in the alignment to improve mapping sensitivity.

            See https://www.mdpi.com/2073-4425/9/3/157

    :type species: str
    :type sample: str
    :type accession: str
    :type paired: bool
    """
    species = luigi.Parameter()
    sample = luigi.Parameter()
    accession = luigi.Parameter()
    paired = luigi.BoolParameter(default=True)

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
    paired = luigi.BoolParameter(default=True)  # TODO this should come from the CSV

    resources = {'cpu-cores': CPU_CORES_LOW, 'ram-gb': 8}  # TODO check ram usage

    def requires(self):
        yield AdapterRemoval(self.accession, self.paired)
        yield BwaAln(self.species, self.sample, self.accession)
        yield ReferenceFASTA(self.species)

    def output(self):
        return [luigi.LocalTarget('data/bam/{}.sort.{}'.format(self.accession, ext)) for ext in ['bam', 'bam.bai']]

    def run(self):
        # unpack the params
        (fastq_file, _, _), sai_file, (ref_file, _) = self.input()
        bam_out, _ = self.output()

        params = {
            'readgroup': r'@RG\tID:{sample}\tSM:{sample}'.format(sample=self.sample),
            'reference': ref_file.path,
            'sai':       sai_file.path,
            'fastq':     fastq_file.path,
            'threads':   self.resources['cpu-cores'],
        }

        with bam_out.temporary_path() as bam_path:
            # get the temporary path for the bam file
            params['bam'] = bam_path

            # align using the ALN algorithm, and convert to a sorted BAM file
            cmd = "bwa samse -r '{readgroup}' {reference} {sai} {fastq} " \
                  " | samtools sort -@ {threads} -O bam -o {bam} -".format(**params)

            # perform the alignment
            utils.run_cmd([cmd], shell=True)

        # index the BAM file
        utils.run_cmd(['samtools', 'index', '-b', bam_out.path])


class FilterUniqueSAMCons(utils.PipelineTask):
    """
    Remove PCR duplicates, so we don't overestimate coverage.

    FilterUniqueSAMCons calls consensus on bases in duplicate reads, rather than simply keeping the read with best mapq.

    See https://bioinf.eva.mpg.de/fastqProcessing/

    :type species: str
    :type sample: str
    :type accession: str
    """
    species = luigi.Parameter()
    sample = luigi.Parameter()
    accession = luigi.Parameter()

    resources = {'cpu-cores': 1, 'ram-gb': 8}  # TODO check ram usage

    def requires(self):
        return BwaSamSe(self.species, self.sample, self.accession)

    def output(self):
        return [luigi.LocalTarget('data/bam/{}.sort.rmdup.{}'.format(self.accession, ext)) for ext in
                ['bam', 'bam.bai']]

    def run(self):
        # unpack the params
        bam_in, _ = self.input()
        bam_out, _ = self.output()

        with bam_out.temporary_path() as bam_path:
            # filter duplicates
            cmd = 'samtools view -h {} | FilterUniqueSAMCons.py | samtools view -b -o {}'.format(bam_in.path, bam_path)

            utils.run_cmd([cmd], shell=True)

        # index the BAM file
        utils.run_cmd(['samtools', 'index', '-b', bam_out.path])


class AncientAlignPipeline(utils.PipelineWrapperTask):
    """
    Test the ancient alignment pipeline.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        for accession in ['ERR2528430', 'ERR2528431', 'ERR2528432']:
            yield BwaSamSe('horse', 'Blagotin2', accession, paired=False)


if __name__ == '__main__':
    luigi.run()
