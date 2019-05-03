#!/usr/bin/env python
# -*- coding: utf-8 -*-

import luigi
import pysam
import os
import random
import glob

from collections import defaultdict

from alleletraj.qtl.qtls import PopulateAllLoci, MIN_DAF
from alleletraj.modern.load_snps import ModernSNPsPipeline
from alleletraj.pipeline_ensembl import EnsemblPipeline
from alleletraj.modern.alignment import ReferenceFASTA
from alleletraj.modern.snp_call import ReferencePloidy, MIN_GENO_QUAL
from samples import LoadSamples
from alleletraj.utils import PipelineTask, PipelineWrapperTask, run_cmd

# minimum depth of coverage to call diploid genotypes
MIN_GENO_DEPTH = 10

# number of bases to hard clip
HARD_CLIP_DIST = 5

# minimum mapping quality (hard filtered)
HARD_MAPQ_CUTOFF = 30

# minimum base quality (hard filtered)
HARD_BASEQ_CUTOFF = 30


class MergeAllLoci(PipelineTask):
    """
    Merge all the overlapping QTL and pseudo-QTL windows.

    :type species: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        yield PopulateAllLoci(self.species)

    def output(self):
        return luigi.LocalTarget('bed/{}-loci.bed'.format(self.basename))

    def run(self):
        # open a db connection
        dbc = self.db_conn()

        # get all the unique QTL windows
        loci = dbc.get_records_sql("""
            SELECT DISTINCT q.chrom, q.start, q.end
              FROM qtls q
             WHERE q.chrom = '{chrom}' 
               AND q.valid = 1
          ORDER BY q.start, q.end""".format(chrom=self.chrom), key=None)

        tmp_loci = 'bed/{}-tmp-loci.bed'.format(self.basename)

        # write all the QTL regions to a BED file
        with open(tmp_loci, 'w') as fout:
            for locus in loci:
                # NOTE that BED starts are zero-based and BED ends are one-based
                fout.write('{}\t{}\t{}\n'.format(locus['chrom'], locus['start']-1, locus['end']))

        # now merge overlapping loci
        bed = run_cmd(['bedtools', 'merge', '-i', tmp_loci])

        # tidy up the tmp file
        os.remove(tmp_loci)

        # save the merged loci
        with self.output().open('w') as fout:
            fout.write(bed)


# TODO align the ancient data


class LoadAncientSNPs(PipelineTask):
    """
    Load all the ancient data for SNPs that fall within the loci of interest.

    :type species: str
    :type chrom: str
    """
    species = luigi.Parameter()
    chrom = luigi.Parameter()

    def requires(self):
        yield ReferenceFASTA(self.species)
        yield ReferencePloidy(self.species)
        yield MergeAllLoci(self.species, self.chrom)
        yield LoadSamples(self.species)
        yield ModernSNPsPipeline(self.species)
        yield EnsemblPipeline(self.species)

    def output(self):
        return luigi.LocalTarget('db/{}-{}.log'.format(self.basename, self.classname))

    # noinspection SqlResolve
    def run(self):
        # unpack the params
        (ref_file, _), pld_file, bed_file, _, _, _ = self.input()
        log_file = self.output()

        # open a db connection
        dbc = self.db_conn()

        # get all the samples and their BAM files
        samples = dbc.get_records_sql("""
            SELECT s.*, GROUP_CONCAT(sf.path) paths
             FROM samples s
             JOIN sample_files sf
               ON sf.sample_id = s.id
            WHERE s.valid = 1
         GROUP BY s.id""")

        # fix issue with weird chars in accession code
        for sid in samples:
            samples[sid]['accession'] = samples[sid]['accession'].encode('utf-8')

        log = log_file.open('w')
        fin = bed_file.open('r')

        # iterate over the loci in the BED file
        for locus in fin:
            chrom, start, end = locus.split()

            # get all the modern SNPs in this locus
            snps = dbc.get_records_sql("""
                SELECT ms.site, ms.ancestral, ms.derived, ev.ref, ev.alt
                  FROM modern_snps ms
                  JOIN modern_snp_daf msd
                    ON msd.modsnp_id = ms.id
             LEFT JOIN ensembl_variants ev
                    ON ev.id = ms.variant_id
                 WHERE ms.chrom = '{chrom}'
                   AND ms.site BETWEEN {start} AND {end}
                   AND msd.daf >= {daf}""".format(chrom=chrom, start=start, end=end, daf=MIN_DAF), key='site')

            log.write("INFO: Scanning locus chr{}:{}-{} for {:,} SNPs".format(chrom, start, end, len(snps)))

            # not all SNPs have a dbsnp entry, so we need to scan the reference to find which alleles are REF/ALT
            # because bcftools needs this info to constrain the diploid genotype calls
            with pysam.FastaFile(ref_file.path) as fasta_file:
                for site in snps:
                    if not snps[site]['ref']:
                        # fetch the ref and snp alleles
                        ref_allele = fasta_file.fetch(chrom, site-1, site)
                        snp_alleles = [snps[site]['ancestral'], snps[site]['derived']]

                        if ref_allele not in snp_alleles:
                            log.write("WARNING: chr{}:{} REF allele {} not found in SNP alleles {}"
                                      .format(chrom, site, ref_allele, snp_alleles))

                            # remove this site
                            snps.pop(site)
                            continue

                        snp_alleles.remove(ref_allele)

                        # set the REF/ALT alleles
                        snps[site]['ref'] = ref_allele
                        snps[site]['alt'] = snp_alleles.pop()

            # the column headers for batch inserting into the db
            fields = ('sample_id', 'chrom', 'site', 'base', 'mapq', 'baseq', 'dist')

            num_reads = 0

            # randomise the order of samples to reduce disk I/O for parallel jobs
            sample_ids = samples.keys()
            random.shuffle(sample_ids)

            # check all the samples for coverage in this locus
            for sample_id in sample_ids:

                # get the sample record
                sample = samples[sample_id]

                log.write("INFO: Scanning locus chr{}:{}-{} in sample {}"
                          .format(chrom, start, end, sample['accession']))

                # buffer the reads so we can bulk insert them into the db
                reads = defaultdict(list)

                # there may be multiple BAM files for each sample
                for path in sample['paths'].split(','):

                    # open the BAM file for reading
                    with pysam.AlignmentFile(path, 'rb') as bam_file:

                        for pileup_column in bam_file.pileup(chrom, int(start), int(end)):

                            # NOTE PileupColumn.reference_pos is 0 based
                            # see http://pysam.readthedocs.io/en/latest/api.html#pysam.PileupColumn.reference_pos
                            site = pileup_column.reference_pos + 1

                            # skip all non-SNP sites
                            if site not in snps:
                                continue

                            # iterate over all the reads for this site
                            for pileup_read in pileup_column.pileups:

                                # skip alignments that don't have a base at this site (i.e. indels)
                                if pileup_read.is_del or pileup_read.is_refskip:
                                    continue

                                # get the read position
                                read_pos = pileup_read.query_position

                                # get the aligned base for this read
                                base = pileup_read.alignment.query_sequence[read_pos]

                                # get the map quality
                                mapq = pileup_read.alignment.mapping_quality

                                # get the base quality
                                baseq = pileup_read.alignment.query_qualities[read_pos]

                                # get the overall length of the read
                                read_length = len(pileup_read.alignment.query_sequence)

                                # how close is the base to the edge of the read
                                dist = min(read_pos, read_length - read_pos)

                                # setup the record to insert, in this order
                                read = (sample_id, chrom, site, base, mapq, baseq, dist)

                                # store the read so we can batch insert later
                                reads[(chrom, site)].append(read)

                # now we're buffered all the reads for this sample at this locus, lets call diploid genotypes on those
                # that pass our depth threshold
                diploid = [idx for idx in reads if len(reads[idx]) >= MIN_GENO_DEPTH]

                if diploid:
                    log.write("INFO: Calling diploid bases in {:,} sites for sample {}"
                              .format(len(diploid), sample_id))

                    suffix = 'luigi-tmp-{:010}'.format(random.randrange(0, 1e10))

                    # make some temp files
                    vcf_file, sex_file, tsv_file, tgz_file, rgs_file,  = [
                        'vcf/diploid-sample{}-{}.{}'.format(sample_id, suffix, ext) for ext in
                        ['vcf', 'sex', 'tsv', 'tsv.gz', 'rgs']]

                    # sort the diploid positions
                    diploid.sort()

                    # save all the callable positions to a file
                    with open(tsv_file, 'w') as fout:
                        fout.write("\n".join("{}\t{}\t{},{}".format(chrom, site, snps[site]['ref'], snps[site]['alt'])
                                             for (chrom, site) in diploid))

                    # bgzip and index the target file
                    run_cmd(["bgzip -c {} > {}".format(tsv_file, tgz_file)], shell=True, verbose=False)
                    run_cmd(["tabix -s1 -b2 -e2 {}".format(tgz_file)], shell=True, verbose=False)

                    # sample names in the BAM file(s) may not be consistent, so override the @SM code with accession
                    with open(rgs_file, 'w') as fout:
                        for path in sample['paths'].split(','):
                            fout.write('*\t{}\t{}'.format(path, sample['accession']))

                    # bcftools needs the sex specified in a separate file
                    with open(sex_file, 'w') as fout:
                        fout.write('{}\t{}\n'.format(sample['accession'], sample['sex']))

                    params = {
                        'ref': ref_file.path,
                        'reg': '{}:{}-{}'.format(chrom, int(start) + 1, end),  # restrict the callable region
                        'tgz': tgz_file,                                        # only call the specified SNPs
                        'rgs': rgs_file,
                        'bam': ' '.join(sample['paths'].split(',')),            # use all the BAM files
                        'pld': pld_file.path,
                        'sex': sex_file,
                        'vcf': vcf_file
                    }

                    # call bases with bcftools (and drop indels and other junk, but keen non-variant sites)
                    # uses both --region (random access) and --targets (streaming) for optimal speed
                    # see https://samtools.github.io/bcftools/bcftools.html#mpileup
                    cmd = "bcftools mpileup --fasta-ref {ref} --regions {reg} --targets-file {tgz} --read-groups {rgs}"\
                          " --output-type u {bam} | " \
                          "bcftools call --multiallelic-caller --ploidy-file {pld} --samples-file {sex} " \
                          " --targets-file {tgz} --constrain alleles --output-type u | " \
                          "bcftools view --exclude-types indels,bnd,other --exclude INFO/INDEL=1 --output-file {vcf} " \
                          .format(**params)

                    # run the base calling
                    run_cmd([cmd], shell=True, verbose=False)

                    # parse the results with pysam
                    for rec in pysam.VariantFile(vcf_file).fetch():

                        # NOTE VariantRecord.pos is 1 based
                        # https://pysam.readthedocs.io/en/latest/api.html#pysam.VariantRecord.pos
                        site = rec.pos

                        # get the genotype call for this site
                        geno = rec.samples[sample['accession']]['GT']

                        # get the genotype quality
                        genoq = int(rec.qual)

                        if genoq >= MIN_GENO_QUAL:
                            # the genotype is good, so drop the raw reads
                            reads.pop((chrom, site))

                            # decode the GT notation into allele calls (e.g. 0/0, 0/1, 1/1)
                            alleles = [rec.alleles[idx] for idx in geno if idx is not None]

                            for allele in alleles:
                                # compose the read records
                                read = {
                                    'sample_id':   sample_id,
                                    'chrom':       chrom,
                                    'site':        site,
                                    'genoq':       genoq,
                                    'base':        allele
                                }

                                dbc.save_record('sample_reads', read)

                    # delete the temp files
                    for tmp in glob.glob("vcf/*{}*".format(suffix)):
                        os.remove(tmp)

                randcall = []

                for chrom, site in reads:
                    # apply hard filters
                    qual = [read for read in reads[(chrom, site)]
                            if read[fields.index('dist')] >= HARD_CLIP_DIST
                            and read[fields.index('mapq')] >= HARD_MAPQ_CUTOFF
                            and read[fields.index('baseq')] >= HARD_BASEQ_CUTOFF]

                    # call a pseudo-haploid genotype at random
                    if len(qual):
                        randcall.append(random.choice(qual))

                # count the total number of reads
                num_reads += len(randcall)

                # bulk insert all the reads for this sample
                if randcall:
                    dbc.save_records('sample_reads', fields, randcall)

            log.write("INFO: Found {:,} reads for locus chr{}:{}-{}".format(num_reads, chrom, start, end))

        log.close()
        fin.close()


class AncientSNPsPipeline(PipelineWrapperTask):
    """
    Populate the modern_snps table.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        for chrom in self.chromosomes:
            yield LoadAncientSNPs(self.species, chrom)


if __name__ == '__main__':
    luigi.run()
