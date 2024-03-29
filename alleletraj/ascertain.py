#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
from collections import OrderedDict

# third party modules
import luigi
import unicodecsv as csv

# local modules
from alleletraj import utils
from alleletraj.ancient.snps import AncientSNPsPipeline
from alleletraj.ensembl.load import LoadEnsemblGenes, LoadEnsemblVariants, FlagSNPsNearIndels
from alleletraj.qtl.load import MC1R_GENE
from alleletraj.ref import ReferenceFASTA
from alleletraj.snpchip.load import LoadSNPChipVariants

# the number of flanking SNPs (on either side) to include
QTL_FLANK_NUM_SNPS = 3

# the number of selective sweep SNPs to include
SWEEP_NUM_SNPS = 5

# the distance between sweep peaks
SWEEP_PEAK_WIDTH = 1000

# the number of "neutral" SNPs to include in the ascertainment
NUM_NEUTRAL_SNPS = 60000

# the number of "ancestral" SNPs to include in the ascertainment
NUM_ANCESTRAL_SNPS = 30000


class FetchGWASPeaks(utils.MySQLTask):
    """
    Fetch all the GWAS peaks from the QTL db.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return AncientSNPsPipeline(self.species)

    def queries(self):
        self.dbc.execute_sql("""
            INSERT 
              INTO ascertainment (qtl_id, type, rsnumber, chrom, site, ref, alt, chip_name, snp_name)
            SELECT q.id, 'peak',
                   ev.rsnumber, ev.chrom, ev.start AS site, ev.ref, ev.alt,
                   sc.chip_name, GROUP_CONCAT(sc.snp_name) AS snp_name
              FROM (
                      # get a unique list of GWAS peaks     
                      SELECT min(id) AS id, peak, site
                        FROM qtls
                       WHERE associationType = 'Association'
                         AND valid = 1
                    GROUP BY peak 
    
                   ) AS q
              JOIN ensembl_variants ev
                ON ev.rsnumber = q.peak
         LEFT JOIN snpchip sc
                ON sc.rsnumber = ev.rsnumber
          GROUP BY ev.rsnumber""")


class FetchGWASFlankingSNPs(utils.MySQLTask):
    """
    Fetch the best flanking SNPs for each GWAS peak.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return FlagSNPsNearIndels(self.species, self.chrom)

    # noinspection SqlResolve
    def queries(self):
        # get the best flanking SNPs for each QTL
        qtls = self.dbc.get_records_sql("""
            SELECT q.id, q.chrom, q.site,
                    
                    SUBSTRING_INDEX(
                        GROUP_CONCAT(
                            IF(ms.site < q.site, ev.id, NULL) 
                              ORDER BY ms.snpchip_id IS NULL, qs.num_reads DESC, RAND()
                        ), ',',  {num_snps}) left_flank,
                    
                    SUBSTRING_INDEX(
                        GROUP_CONCAT(
                            IF(ms.site > q.site, ev.id, NULL) 
                              ORDER BY ms.snpchip_id IS NULL, qs.num_reads DESC, RAND()
                        ), ',',  {num_snps}) right_flank
              FROM (
                        # get a unique list of GWAS peaks
                        SELECT min(id) AS id, chrom, peak, site
                          FROM qtls
                         WHERE associationType = 'Association'
                           AND valid = 1
                      GROUP BY peak
                      
                   ) AS q
              JOIN qtl_snps qs
                ON qs.qtl_id = q.id
              JOIN modern_snps ms
                ON ms.id = qs.modsnp_id
              JOIN ensembl_variants ev
                ON ev.id = ms.variant_id
               AND (ev.indel IS NULL OR ms.snpchip_id IS NOT NULL)
          GROUP BY q.id""".format(num_snps=QTL_FLANK_NUM_SNPS), key=None)

        # we have to do this iteratively, as FIND_IN_SET() performs terribly
        for qtl in qtls:
            # merge the flanking SNP modsnp_ids
            modsnps = ','.join(flank for flank in [qtl['left_flank'], qtl['right_flank']] if flank)

            self.dbc.execute_sql("""
                INSERT 
                  INTO ascertainment (qtl_id, type, rsnumber, chrom, site, ref, alt, chip_name, snp_name)
                SELECT {qtl_id}, IF(ev.start > {site}, 'right', 'left'),
                       ev.rsnumber, ev.chrom, ev.start AS site, ev.ref, ev.alt,
                       sc.chip_name, GROUP_CONCAT(sc.snp_name) AS snp_name
                  FROM ensembl_variants ev
             LEFT JOIN snpchip sc
                    ON sc.rsnumber = ev.rsnumber
                 WHERE ev.type = 'SNV'
                   AND ev.id IN ({modsnps})
              GROUP BY ev.rsnumber
                   """.format(qtl_id=qtl['id'], site=qtl['site'], modsnps=modsnps))


class FetchSelectiveSweepSNPs(utils.MySQLTask):
    """
    Fetch the best SNPs from the selective sweep loci.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        for chrom in self.chromosomes:
            yield FlagSNPsNearIndels(self.species, chrom)

    # noinspection SqlResolve
    def queries(self):
        qtls = self.dbc.get_records_sql("""
            SELECT qtl_id AS id,
                   SUBSTRING_INDEX(
                       GROUP_CONCAT(modsnp_id ORDER BY ss.p, ABS(ss.site-near.site)), 
                       ',',  {num_snps}) AS snps
              FROM (
                        # find the nearest sweep SNP to each QTL SNP in dbsnp 
                        SELECT ms.id AS modsnp_id,
                               ms.site,
                               SUBSTRING_INDEX(
                                   GROUP_CONCAT(ss.id ORDER BY ABS(ss.site-ms.site)), 
                                   ',',  1) AS ss_id
                          FROM sweep_snps ss
                          JOIN qtl_snps qs
                            ON qs.qtl_id = ss.qtl_id 
                          JOIN modern_snps ms
                            ON ms.id = qs.modsnp_id
                          JOIN ensembl_variants ev
                            ON ev.id = ms.variant_id
                           AND (ev.indel IS NULL OR ms.snpchip_id IS NOT NULL)
                      GROUP BY ms.id
                      
                   ) AS near
              JOIN sweep_snps ss
                ON ss.id = near.ss_id
          GROUP BY ss.qtl_id
               """.format(num_snps=SWEEP_NUM_SNPS, offset=SWEEP_PEAK_WIDTH / 2), key=None)

        # we have to do this iteratively, as FIND_IN_SET() performs terribly
        for qtl in qtls:
            self.dbc.execute_sql("""
                INSERT 
                  INTO ascertainment (qtl_id, type, rsnumber, chrom, site, ref, alt, chip_name, snp_name)
                SELECT {qtl_id}, 'sweep', 
                       ev.rsnumber, ev.chrom, ev.start AS site, ev.ref, ev.alt,
                       sc.chip_name, GROUP_CONCAT(sc.snp_name) AS snp_name
                  FROM modern_snps ms
             LEFT JOIN ensembl_variants ev
                    ON ev.id = ms.variant_id
             LEFT JOIN snpchip sc
                    ON sc.rsnumber = ev.rsnumber
                 WHERE ms.id IN ({modsnps})
              GROUP BY ms.id
                   """.format(qtl_id=qtl['id'], modsnps=qtl['snps']))


class FetchMC1RSNPs(utils.MySQLTask):
    """
    Get all the dnsnp SNPs which fall within the MC1R gene.

    See https://www.ensembl.org/sus_scrofa/Gene/Summary?g=ENSSSCG00000020924&db=core

    Also http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000341
         http://journals.plos.org/plosgenetics/article/file?id=10.1371/journal.pgen.1000341.s002&type=supplementary

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        yield LoadEnsemblGenes(self.species)
        yield LoadEnsemblVariants(self.species)
        yield AncientSNPsPipeline(self.species)

    # noinspection SqlResolve
    def queries(self):
        # get the MC1R qtl
        qtl = self.dbc.get_record('qtls', {'associationType': 'MC1R'})

        self.dbc.execute_sql("""
            INSERT 
              INTO ascertainment (qtl_id, type, rsnumber, chrom, site, ref, alt, chip_name, snp_name)
            SELECT {qtl_id}, 'mc1r', 
                   ev.rsnumber, ev.chrom, ev.start AS site, ev.ref, ev.alt,
                   sc.chip_name, GROUP_CONCAT(sc.snp_name) AS snp_name
              FROM ensembl_genes eg
              JOIN ensembl_variants ev
                ON ev.chrom = eg.chrom
               AND ev.start BETWEEN eg.start AND eg.end
               AND ev.type = 'SNV'
               AND LENGTH(ev.alt) = 1
         LEFT JOIN snpchip sc
                ON sc.rsnumber = ev.rsnumber
             WHERE eg.gene_name = '{gene_name}'
          GROUP BY ev.rsnumber
               """.format(qtl_id=qtl['id'], gene_name=MC1R_GENE))


class FetchNeutralSNPs(utils.MySQLTask):
    """
    Get neutral SNPs (excluding all QTLs and gene regions, w/ buffer).

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        yield ReferenceFASTA(self.species)
        yield LoadEnsemblVariants(self.species)
        yield AncientSNPsPipeline(self.species)

    # noinspection SqlResolve
    def queries(self):
        # unpack the inputs
        (_, fai_file), _, _ = self.input()

        # get the sizes of the chromosomes
        sizes = utils.get_chrom_sizes(fai_file)

        # get the total size of the autosomes
        total = sum(sizes[chrom] for chrom in self.autosomes)

        # calculate the proportional size of each autosomal chrom
        percent = OrderedDict((chrom, float(sizes[chrom]) / total) for chrom in self.autosomes)

        for chrom in percent:
            # get the weighted number of SNPs for this chrom
            num_snps = int(round(NUM_NEUTRAL_SNPS * percent[chrom]))

            self.dbc.execute_sql("""
                INSERT
                  INTO ascertainment (qtl_id, type, rsnumber, chrom, site, ref, alt, chip_name, snp_name)
                SELECT q.id, q.associationType,
                       ev.rsnumber, ev.chrom, ev.start AS site, ev.ref, ev.alt,
                       sc.chip_name, GROUP_CONCAT(DISTINCT sc.snp_name) AS snp_name
                  FROM qtls q
                  JOIN qtl_snps qs
                    ON q.id = qs.qtl_id
                  JOIN modern_snps ms
                    ON ms.id = qs.modsnp_id
                  JOIN ensembl_variants ev
                    ON ev.id = ms.variant_id
                   AND (ev.indel IS NULL OR ms.snpchip_id IS NOT NULL)
             LEFT JOIN snpchip sc
                    ON sc.rsnumber = ev.rsnumber
                 WHERE q.chrom = '{chrom}'
                   AND q.associationType = 'Neutral'
                   AND q.valid = 1
              GROUP BY ms.id
              ORDER BY ms.snpchip_id IS NULL, qs.num_reads DESC, rand()
                 LIMIT {num_snps}
                   """.format(chrom=chrom, num_snps=num_snps))


class FetchAncestralSNPs(utils.MySQLTask):
    """
    Get ancestral SNPs which are variable in ASD (Asian domestic) and SUM (Sumatran Sus scrofa) populations.

    Prefer SNPs which are on the SNPchip and those already chosen in a different category.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        yield ReferenceFASTA(self.species)
        yield AncientSNPsPipeline(self.species)

    # noinspection SqlResolve
    def queries(self):
        # unpack the inputs
        (_, fai_file), _ = self.input()

        # get the sizes of the chromosomes
        sizes = utils.get_chrom_sizes(fai_file)

        # get the total size of the autosomes
        total = sum(sizes[chrom] for chrom in self.autosomes)

        # calculate the proportional size of each autosomal chrom
        percent = OrderedDict((chrom, float(sizes[chrom]) / total) for chrom in self.autosomes)

        for chrom in percent:
            # get the weighted number of SNPs for this chrom
            num_snps = int(round(NUM_ANCESTRAL_SNPS * percent[chrom]))

            # TODO doesn't work anymore, because modern_snps doesn't have a population column
            self.dbc.execute_sql("""
                INSERT
                  INTO ascertainment (qtl_id, type, rsnumber, chrom, site, ref, alt, chip_name, snp_name)
                SELECT NULL, 'Ancestral',
                       ev.rsnumber, ev.chrom, ev.start AS site, ev.ref, ev.alt,
                       sc.chip_name, GROUP_CONCAT(DISTINCT sc.snp_name) AS snp_name
                  FROM modern_snps asd
                  JOIN modern_snps sum
                    ON sum.population = 'SUM'  # TODO fix me
                   AND sum.chrom = asd.chrom
                   AND sum.variant_id = asd.variant_id
                  JOIN ensembl_variants ev
                    ON ev.id = asd.variant_id
                   AND (ev.indel IS NULL OR asd.snpchip_id IS NOT NULL)
             LEFT JOIN snpchip sc
                    ON sc.rsnumber = ev.rsnumber 
             LEFT JOIN ascertainment a
                    ON asd.chrom = a.chrom
                   AND asd.site = a.site
                 WHERE asd.population = 'ASD'  # TODO fix me
                   AND asd.chrom = '{chrom}'
                   AND asd.derived_count > 1
              GROUP BY asd.id
              ORDER BY asd.snpchip_id IS NULL, a.id IS NULL, rand()
                 LIMIT {num_snps}
                   """.format(chrom=chrom, num_snps=num_snps))


class FetchRemainingSNPChipSNPs(utils.MySQLTask):
    """
    Include any SNPchip SNPs which were not already included in a previous ascertainment category.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        yield LoadEnsemblVariants(self.species)
        yield LoadSNPChipVariants(self.species)

    def queries(self):
        self.dbc.execute_sql("""
            INSERT
              INTO ascertainment (qtl_id, type, rsnumber, chrom, site, ref, alt, chip_name, snp_name)
            SELECT NULL, 'snpchip',
                   ev.rsnumber, ev.chrom, ev.start AS site, ev.ref, ev.alt,
                   sc.chip_name, GROUP_CONCAT(DISTINCT sc.snp_name) AS snp_name
              FROM snpchip sc
              JOIN ensembl_variants ev
                ON ev.rsnumber = sc.rsnumber
         LEFT JOIN ascertainment a
                ON a.chrom = ev.chrom
               AND a.site = ev.start
             WHERE a.id IS NULL
               AND ev.type = 'SNV'
               AND length(ev.alt) = 1
          GROUP BY ev.rsnumber""")


class PerformAscertainment(utils.PipelineWrapperTask):
    """
    Ascertain the best SNPs for our capture array.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        # fetch all the GWAS peaks from the QTL db
        yield FetchGWASPeaks(self.species)

        # fetch the best flanking SNPs for each GWAS peak
        yield FetchGWASFlankingSNPs(self.species)

        # fetch the best SNPs from the selective sweep loci
        yield FetchSelectiveSweepSNPs(self.species)

        # get all MC1R snps
        yield FetchMC1RSNPs(self.species)

        # get neutral SNPs (excluding all QTLs and gene regions, w/ buffer)
        yield FetchNeutralSNPs(self.species)

        # get ancestral SNPs which are in variable in ASD and Sumatran scrofa
        yield FetchAncestralSNPs(self.species)

        # include any snpchip SNPs which are not already included
        yield FetchRemainingSNPChipSNPs(self.species)


class ExportAscertainedSNPs(utils.DatabaseTask):
    """
    Export all the ascertained SNPs to a TSV.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return PerformAscertainment(self.species)

    def output(self):
        return luigi.LocalTarget('data/tsv/{}.candidate-snps.tsv'.format(self.species))

    def run(self):
        # get all the unique SNPs in the ascertainment
        reads = self.dbc.get_records_sql("""
            SELECT GROUP_CONCAT(DISTINCT qtl_id ORDER BY qtl_id) AS qtls,
                   GROUP_CONCAT(DISTINCT type ORDER BY type) AS types,
                   rsnumber, chrom, site, ref, alt, chip_name, 
                   GROUP_CONCAT(DISTINCT snp_name ORDER BY snp_name) AS snp_name
              FROM ascertainment
          GROUP BY rsnumber
          ORDER BY chrom, site
               """.format(), key=None)

        with self.output().open('w') as fout:
            writer = csv.DictWriter(fout, fieldnames=reads[0].keys(), delimiter='\t')
            writer.writeheader()

            # write the data to disk
            for read in reads:
                writer.writerow(read)


class AscertainmentPipeline(utils.PipelineWrapperTask):
    """
    Perform the ascertainment for the pig SNP capture.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return ExportAscertainedSNPs(self.species)


if __name__ == '__main__':
    luigi.run()
