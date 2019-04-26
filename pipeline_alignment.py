#!/usr/bin/env python
# -*- coding: utf-8 -*-

import luigi

# import my custom modules
from pipeline_consts import CPU_CORES_MED, BAM_FILES
from pipeline_utils import PipelineTask, PipelineExternalTask, PipelineWrapperTask, run_cmd, trim_ext

# hard filters for TrimGalore!
TRIM_MIN_BASEQ = 20
TRIM_MIN_LENGTH = 25

# location of software tools
GATK = "/usr/local/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar"
PICARD = "/usr/local/picard-tools-2.5.0/picard.jar"


class SraToolsFastqDump(PipelineTask):
    """
    Fetches the paired-end and single end FASTQ files for a given accession code, using SRA Tools

    :type accession: str
    :type paired: bool
    """
    accession = luigi.Parameter()
    paired = luigi.BoolParameter()

    def output(self):
        if self.paired:
            yield luigi.LocalTarget('fastq/{}_1.fastq.gz'.format(self.accession))
            yield luigi.LocalTarget('fastq/{}_2.fastq.gz'.format(self.accession))
        else:
            yield [luigi.LocalTarget('fastq/{}.fq.gz'.format(self.accession))]

        yield luigi.LocalTarget('fastq/{}.log'.format(self.accession))

    def run(self):
        log_file = list(self.output()).pop()

        # use the NCBI SRA toolkit to fetch the fastq files
        log = run_cmd(['fastq-dump',
                       '--gzip',               # output gzipped files
                       '--split-3',            # split into two paired end fastq files + one unpaired fastq
                       '--outdir', './fastq',  # output directory
                       self.accession])

        with log_file.open('w') as fout:
            fout.write(log)


class TrimGalore(PipelineTask):
    """
    Trim adapters, low-quality bases from the 3' end, and drop any resulting reads below a min-length threshold.

    :type accession: str
    :type paired: bool
    """
    accession = luigi.Parameter()
    paired = luigi.BoolParameter()

    def requires(self):
        return SraToolsFastqDump(self.accession)

    def output(self):
        if self.paired:
            for pair in [1, 2]:
                yield luigi.LocalTarget('fastq/{}_{}_val_{}.fq.gz'.format(self.accession, pair, pair))
                yield luigi.LocalTarget('fastq/{}_{}_val_{}_fastqc.zip'.format(self.accession, pair, pair))
        else:
            yield luigi.LocalTarget('fastq/{}_trimmed.fq.gz'.format(self.accession))
            yield luigi.LocalTarget('fastq/{}_fastqc.zip'.format(self.accession))

    def run(self):
        fastq_files = self.input()[:-1]

        cmd = ['trim_galore',
               '--quality', TRIM_MIN_BASEQ,   # trim low-quality ends from reads in addition to adapter removal
               '--length',  TRIM_MIN_LENGTH,  # discard reads that became shorter than length because of trimming
               '--fastqc',                    # run FastQC once trimming is complete
               '--output_dir', 'fastq']

        if self.paired:
            cmd.append('--paired')

        # add the fastq files
        for fastq in fastq_files:
            cmd.append(fastq.path)

        # perform the trimming
        run_cmd(cmd)


class ReferenceFASTA(PipelineTask):
    """
    Get the reference genome and index it.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        # avoid circular dependency
        from pipeline_ensembl import DownloadEnsemblData
        return DownloadEnsemblData(self.species, 'fasta', bgzip=True)

    def output(self):
        return [luigi.LocalTarget('ensembl/{}.{}.dna.toplevel.{}'.format(self.binomial, self.assembly, ext)) for ext in
                ['fa.gz', 'fa.gz.fai']]

    def run(self):
        # get the downloaded reference assembly
        ref_file = self.input()

        # build an index
        run_cmd(['samtools', 'faidx', ref_file.path])


class BwaIndexBWTSW(PipelineTask):
    """
    Builds the BWA index for the reference genome, needed for performing alignments

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return ReferenceFASTA(self.species)

    def output(self):
        ref_file, _ = self.input()
        return [luigi.LocalTarget('{}.{}'.format(ref_file.path, ext)) for ext in ['amb', 'ann', 'bwt', 'pac', 'sa']]

    def run(self):
        ref_file, _ = self.input()

        run_cmd(['bwa',
                 'index',        # index needed for bwa alignment
                 '-a', 'bwtsw',  # algorithm suitable for mammals
                 ref_file.path])


class BwaMem(PipelineTask):
    """
    Align a fastq file to the reference genome

    :type species: str
    :type accession: str
    :type paired: bool
    """
    species = luigi.Parameter()
    accession = luigi.Parameter()
    paired = luigi.BoolParameter(default=True)

    resources = {'cpu-cores': CPU_CORES_MED, 'ram-gb': 8}

    def requires(self):
        yield TrimGalore(self.accession, self.paired)
        yield ReferenceFASTA(self.species)
        yield BwaIndexBWTSW(self.species)

    def output(self):
        return [luigi.LocalTarget('bam/{}.sort.{}'.format(self.accession, ext)) for ext in ['bam', 'bam.bai']]

    def run(self):
        # unpack the input params
        fastq_input, (ref_file, _), _ = self.input()
        bam_out, _ = self.output()

        params = {
            # CPU threads to use
            'threads': self.resources['cpu-cores'],

            # compose the read group metadata
            'readgroup': r'@RG\tID:{sample}\tSM:{sample}'.format(sample=self.accession),

            # reference genome
            'reference': ref_file.path,

            # get the fastq file(s) to align
            'fastq': ' '.join(fastq.path for fastq in fastq_input if '.fq.gz' in fastq.path),
        }

        with bam_out.temporary_path() as bam_path:

            # get the temporary path for the bam file
            params['bam'] = bam_path

            # align using the mem algorithm, and convert to a sorted BAM file
            cmd = "bwa mem -t {threads} -R '{readgroup}' {reference} {fastq} " \
                  " | samtools sort -@ {threads} -O bam -o {bam} -".format(**params)

            # perform the alignment
            run_cmd([cmd], shell=True)

        # index the BAM file
        run_cmd(['samtools', 'index', '-b', bam_out.path])


class PicardMarkDuplicates(PipelineTask):
    """
    Remove PCR duplicates, so we don't overestimate coverage

    :type species: str
    :type accession: str
    """
    species = luigi.Parameter()
    accession = luigi.Parameter()

    resources = {'cpu-cores': 1, 'ram-gb': 8}

    def requires(self):
        return BwaMem(self.species, self.accession)

    def output(self):
        return [luigi.LocalTarget('bam/{}.sort.rmdup.{}'.format(self.accession, ext)) for ext in
                ['bam', 'bam.bai', 'log']]

    def run(self):
        # unpack the params
        bam_in, _ = self.input()
        bam_out, _, log_file = self.output()

        with bam_out.temporary_path() as bam_path:

            run_cmd(['java', self.java_mem,
                     '-jar', PICARD,
                     'MarkDuplicates',
                     'INPUT=' + bam_in.path,
                     'OUTPUT=' + bam_path,
                     'METRICS_FILE=' + log_file.path,
                     'REMOVE_DUPLICATES=true',
                     'QUIET=true'])

        # index the BAM file
        run_cmd(['samtools', 'index', '-b', bam_out.path])


class PicardSequenceDictionary(PipelineTask):
    """
    Unzip the reference genome and create the sequence dictionary, because GATK is stupid and cannot handle bgzip.

    :type species: str
    """
    species = luigi.Parameter()

    resources = {'cpu-cores': 1, 'ram-gb': 8}

    def requires(self):
        return ReferenceFASTA(self.species)

    def output(self):
        return [luigi.LocalTarget('ensembl/{}.{}.dna.toplevel.{}'.format(self.binomial, self.assembly, ext)) for ext in
                ['fa', 'fa.fai', 'dict']]

    def run(self):
        ref_in, _ = self.input()
        ref_out, _, dict_file = self.output()

        # unzip the reference genome
        run_cmd(['gunzip', '--keep', '--force', ref_in.path])

        # build a regular index
        run_cmd(['samtools', 'faidx', ref_out.path])

        # create the sequence dictionary
        with dict_file.temporary_path() as dict_path:
            run_cmd(['java', self.java_mem,
                     '-jar', PICARD,
                     'CreateSequenceDictionary',
                     'R=' + ref_out.path,
                     'O=' + dict_path])


class GATKRealignerTargetCreator(PipelineTask):
    """
    Use the GATK RealignerTargetCreator to find suspicious indels to target for local realignment

    :type species: str
    :type accession: str
    """
    species = luigi.Parameter()
    accession = luigi.Parameter()

    resources = {'cpu-cores': 1, 'ram-gb': 4}

    def requires(self):
        yield PicardMarkDuplicates(self.species, self.accession)
        yield PicardSequenceDictionary(self.species)

    def output(self):
        return [luigi.LocalTarget('bam/{}.sort.rmdup.realign.{}'.format(self.accession, ext)) for ext in
                ['intervals', 'intervals.log']]

    def run(self):
        # unpack the inputs/outputs
        (bam_in, _, _), (ref_file, _, _) = self.input()
        itv_file, log_file = self.output()

        with itv_file.temporary_path() as itv_path, open(log_file.path, 'w') as log_fout:

            run_cmd(['java', self.java_mem,
                     '-jar', GATK,
                     '--analysis_type', 'RealignerTargetCreator',
                     '--reference_sequence', ref_file.path,
                     '--input_file', bam_in.path,
                     '--out', itv_path],
                    stdout=log_fout)


class GATKIndelRealigner(PipelineTask):
    """
    Use the GATK IndelRealigner to fix perform local realignment of reads around indels

    :type species: str
    :type accession: str
    """
    species = luigi.Parameter()
    accession = luigi.Parameter()

    resources = {'cpu-cores': 1, 'ram-gb': 4}

    def requires(self):
        yield PicardMarkDuplicates(self.species, self.accession)
        yield GATKRealignerTargetCreator(self.species, self.accession)
        yield PicardSequenceDictionary(self.species)

    def output(self):
        return [luigi.LocalTarget('bam/{}.sort.rmdup.realign.{}'.format(self.accession, ext)) for ext in
                ['bam', 'bai', 'log']]

    def run(self):
        # unpack the inputs/outputs
        (bam_in, _, _), (itv_file, _), (ref_file, _, _) = self.input()
        bam_out, _, log_file = self.output()

        with bam_out.temporary_path() as bam_path, open(log_file.path, 'w') as log_fout:

            run_cmd(['java', self.java_mem,
                     '-jar', GATK,
                     '--analysis_type', 'IndelRealigner',
                     '--reference_sequence', ref_file.path,
                     '--input_file', bam_in.path,
                     '--targetIntervals', itv_file.path,
                     '--out', bam_path],
                    stdout=log_fout)


class SAMToolsMerge(PipelineTask):
    """
    Merge multiple libraries into a single BAM file

    :type species: str
    :type sample: str
    """
    species = luigi.Parameter()
    sample = luigi.Parameter()

    resources = {'cpu-cores': 1, 'ram-gb': 8}

    def requires(self):
        for accession in self.accessions:
            yield GATKIndelRealigner(self.species, accession)

    def output(self):
        return [luigi.LocalTarget('bam/{}.sort.rmdup.realign.{}'.format(self.sample, ext)) for ext in ['bam', 'bai']]

    def run(self):
        bam_inputs = [bam_file.path for bam_file, bai_file, log_file in self.input()]
        bam_out, _ = self.output()

        with bam_out.temporary_path() as bam_path:

            run_cmd(['samtools',
                     'merge',
                     bam_path
                     ] + bam_inputs)

        # index the BAM file
        run_cmd(['samtools', 'index', '-b', bam_out.path])


class ExternalBAM(PipelineExternalTask):
    """
    External task dependency for an aligned BAM file.

    N.B. These have been created outside the workflow of this pipeline.

    :type species: str
    :type sample: str
    """
    species = luigi.Parameter()
    sample = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(BAM_FILES[self.species][self.sample])


class AlignedBAM(PipelineWrapperTask):
    """
    Align raw reads to a reference genome, deduplicate, realign indels, and merge multiple libraries.

    :type species: str
    :type sample: str
    """
    species = luigi.Parameter()
    sample = luigi.Parameter()

    def requires(self):
        if self.species == 'pig':
            # use the BAM file we aligned ourselves
            yield SAMToolsMerge(self.species, self.sample)
        else:
            # use the BAM file provide by collaborators
            yield ExternalBAM(self.species, self.sample)


class AlignmentPipeline(PipelineWrapperTask):
    """
    Get BAM files for all samples

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        for samples in self.samples:
            yield AlignedBAM(self.species, samples)


if __name__ == '__main__':
    luigi.run()
