#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.ancient.dedupe import FilterUniqueSAMCons
from alleletraj.ancient.rescale import MapDamageRescale
from alleletraj.gatk import GATKIndelRealigner
from alleletraj.modern.dedupe import PicardMarkDuplicates

# ignore these inconsequential errors when running ValidateSamFile
VALIDATE_BAM_IGNORE = [
    'MISSING_PLATFORM_VALUE',         # The read group is missing its PL (platform unit) field
    'MATE_NOT_FOUND',                 # Read is marked as paired, but its pair was not found
    'INVALID_FLAG_FIRST_OF_PAIR',     # First of pair flag set for unpaired read
    'INVALID_FLAG_MATE_UNMAPPED',     # Mate unmapped flag is incorrectly set
    'INVALID_FLAG_SECOND_OF_PAIR',    # Second of pair flag set for unpaired read
    'INVALID_MATE_REF_INDEX',         # Mate reference index (MRNM) set for unpaired read
    'MISMATCH_MATE_ALIGNMENT_START',  # Mate alignment does not match alignment start of mate
    'MISMATCH_MATE_CIGAR_STRING',     # The mate cigar tag does not match its mate's cigar string
]


class ExternalSampleBAM(utils.PipelineExternalTask):
    """
    External task to return a fully processed BAM file for a given sample.

    N.B. These have been created outside the workflow of this pipeline.

    :type path: str
    """
    path = luigi.Parameter()

    def output(self):
        yield luigi.LocalTarget(self.path)
        yield luigi.LocalTarget(self.path + '.bai')


class ValidateBamFile(utils.DatabaseTask):
    """
    Validate an external BAM file using Picard.

    This tool reports on the validity of a BAM file relative to the SAM format specification.

    https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/picard_sam_ValidateSamFile.php
    https://software.broadinstitute.org/gatk/documentation/article.php?id=7571

    :type species: str
    :type population: str
    :type sample: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    sample = luigi.Parameter()

    resources = {'cpu-cores': 1, 'ram-gb': 8}

    # do not retry after failure, as this just chews CPU cycles
    retry_count = 0

    def requires(self):
        return ExternalSampleBAM(self.sample_data['path'])

    def output(self):
        # included the external bam and bai files in the output
        return list(self.input()) + [luigi.LocalTarget('data/bam/{}.{}'.format(self.sample, ext)) for ext in
                                     ['log', 'err']]

    def run(self):
        bam_file, _ = self.input()
        _, _, log_file, err_file = self.output()

        # ValidateSamFile is super pedantic and there are many error types we don't really care about
        ignore = ['IGNORE=' + error for error in VALIDATE_BAM_IGNORE]

        # validate the BAM file
        with log_file.temporary_path() as log_path:
            utils.run_cmd(['java',
                           self.java_mem,
                           self.java_gc_threads,
                           '-jar', 'jar/picard.jar',
                           'ValidateSamFile',
                           'MODE=SUMMARY',
                           'INPUT=' + bam_file.path,
                           'OUTPUT=' + err_file.path
                           ] + ignore, stderr=open(log_path, 'w'))


class AccessionBAM(utils.DatabaseTask):
    """
    Wrapper taks to return a processed BAM file for a given accession code.

    Pipeline forks depending on ancient/modern status of sample.

    :type species: str
    :type sample: str
    :type accession: str
    """
    species = luigi.Parameter()
    sample = luigi.Parameter()
    accession = luigi.Parameter()

    def requires(self):
        if self.sample_data['ancient']:
            return MapDamageRescale(self.species, self.sample, self.accession)
        else:
            return GATKIndelRealigner(self.species, self.sample, self.accession)

    def output(self):
        return self.input()


class DeduplicatedBAM(utils.DatabaseTask):
    """
    Wrapper task to return a deduplicated BAM file.

    We use different trimming, alignment and deduplication methods depending on ancient/modern status of sample.

    When `accession` is not set this returns a deduplicated BAM file containing all the sample accessions.

    :type species: str
    :type sample: str
    :type accession: str
    """
    species = luigi.Parameter()
    sample = luigi.Parameter()
    accession = luigi.OptionalParameter(default=None)

    def complete(self):
        # this task is complete if our input task is complete (necessary to force ValidateBamFile to run)
        return self.requires().complete()

    def requires(self):
        if self.sample_data['ancient']:
            # use the aDNA pipeline
            return FilterUniqueSAMCons(self.species, self.sample, self.accession)
        else:
            # use the modern pipeline
            return PicardMarkDuplicates(self.species, self.sample, self.accession)

    def output(self):
        # only pass on the bam and bai files (i.e. trim off any .log files)
        return self.input()[:2]


class SampleBAM(utils.DatabaseTask):
    """
    Wrapper task to return a fully processed BAM file for a given sample.

    Align raw reads to a reference genome, deduplicate, realign indels, rescale bases, and merge multiple libraries.

    :type species: str
    :type population: str
    :type sample: str
    """
    species = luigi.Parameter()
    population = luigi.Parameter()
    sample = luigi.Parameter()

    def complete(self):
        # this task is complete if our input task is complete (necessary to force ValidateBamFile to run)
        return self.requires().complete()

    def requires(self):
        if self.sample_data.get('path'):
            # validate external BAM files before using them
            return ValidateBamFile(self.species, self.population, self.sample)
        else:
            accessions = self.list_accessions()
            if len(accessions) > 1:
                # merge all accessions together, and deduplicate a second time (just in case)
                return DeduplicatedBAM(self.species, self.sample)
            else:
                return AccessionBAM(self.species, self.sample, accessions[0])

    def output(self):
        # only pass on the bam and bai files (i.e. trim off any .log and .err files from ValidateBamFile)
        return self.input()[:2]


if __name__ == '__main__':
    luigi.run()
