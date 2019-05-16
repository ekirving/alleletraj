#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
import csv
import os
import re
import subprocess
from collections import Iterable, OrderedDict, defaultdict

# third party modules
import luigi

# import my libraries
from alleletraj.const import CHROMOSOMES, REF_ASSEMBLY, BINOMIAL_NAME, OUTGROUP_POP
from alleletraj.db.conn import Database

# enforce max interval size of 1 Gb
MAX_INTERVAL_SIZE = int(1e6)


def run_cmd(cmd, shell=False, background=False, stdout=None, stderr=None, verbose=True):
    """
    Executes the given command in a system subprocess

    :param cmd: The system command to run (list|string)
    :param shell: Use the native shell
    :param background: Flag to tun the process in the background
    :param stdout: File handle to redirect stdout
    :param stderr: File handle to redirect stderr
    :param verbose: Print the command before running it
    :return: The stdout stream
    """
    # subprocess only accepts strings
    cmd = [str(args) for args in cmd]

    if verbose:
        print(u' '.join(cmd))

    if background:
        subprocess.Popen(cmd, shell=shell)

    else:

        stdout = subprocess.PIPE if not stdout else stdout
        stderr = subprocess.PIPE if not stderr else stderr

        # run the command
        proc = subprocess.Popen(cmd, shell=shell, stdout=stdout, stderr=stderr)

        # fetch any output and error
        (out, err) = proc.communicate()

        # bail if something went wrong
        if proc.returncode != 0:

            # decode return codes
            if proc.returncode == 139:
                err = 'Segmentation fault (core dumped) ' + err

            raise RuntimeError(err)

        # some commands log progress to stderr
        return out if out else err


def merge_intervals(ranges, capped=True):
    """
    Merge overlapping intervals, so we only check each site once
    """

    try:
        saved = list(ranges[0])
    except IndexError:
        raise StopIteration

    # sort the intervals by start, and iterate over the list
    for start, end in sorted([sorted(t) for t in ranges]):
        if start <= saved[1]:
            # if the current interval overlaps the saved one then take the largest end point
            saved[1] = max(saved[1], end)
        else:
            if capped:
                # enforce max interval size
                for i in range(saved[0], saved[1], MAX_INTERVAL_SIZE):
                    yield (i, min(saved[1], i + MAX_INTERVAL_SIZE - 1))
            else:
                yield (saved[0], saved[1])

            saved[0] = start
            saved[1] = end

    # return the final interval
    if capped:
        for i in range(saved[0], saved[1], MAX_INTERVAL_SIZE):
            yield (i, min(saved[1], i + MAX_INTERVAL_SIZE - 1))
    else:
        yield (saved[0], saved[1])


def load_samples_csv(csv_file):
    """
    Load a samples csv file as a nested dictionary
    """
    populations = dict()
    with open(csv_file, 'r') as fin:
        data = csv.DictReader(fin)
        data._fieldnames = [re.sub(r'\W+', '', field).lower() for field in data.fieldnames]  # strip bad Excel chars
        for row in data:
            pop = populations.get(row['population'], dict())

            sample = dict()
            for column in row:
                sample[column] = row[column] if column != 'accessions' \
                    else [val for val in row[column].split(';') if val != '']

            pop[row['sample']] = sample
            populations[row['population']] = pop

    return populations


def trim_ext(full_path, n=1):
    return '.'.join(full_path.split('.')[:-n])


def trim_path_ext(full_path):
    return trim_ext(os.path.basename(full_path))


def insert_suffix(full_path, suffix):
    splitpath = full_path.split('.')
    splitpath.insert(-1, suffix)
    return '.'.join(splitpath)


def dump(obj):
    for attr in dir(obj):
        if hasattr(obj, attr):
            print("obj.%s = %s" % (attr, getattr(obj, attr)))


def curl_download(url, filename):
    """
    Downloads a remote url to a local file path using cURL
    """
    run_cmd(['curl', '-s', '--output', filename, url])


def get_chrom_sizes(fai_file, exclude_scaffolds=True):
    """
    Get the names and sizes of all the chromosomes in a reference by iterating over the index
    """
    chroms = {}

    with fai_file.open('r') as fin:
        for line in fin:
            # get the chrom name and size
            contig, size, _, _, _ = line.split()

            # strip any annoying `chr` prefix
            nochr = contig.replace('chr', '')

            # is this contig an actual chromosome
            is_chrom = nochr.isdigit() or nochr in ['X', 'Y', 'Z', 'W', 'MT']

            if exclude_scaffolds and not is_chrom:
                continue

            chroms[nochr] = int(size)

    return chroms


class PipelineTask(luigi.Task):
    """
    Pipeline task that implements several dynamic attributes
    """

    @property
    def resources(self):
        """
        Dynamically set task resource usage
        """
        resources = {'cpu-cores': 1}

        if hasattr(self, 'db_lock_tables'):
            for table in self.db_lock_tables:
                # resolve table names that contain task parameters
                if hasattr(self, 'chrom'):
                    table = table.format(chrom=self.chrom)
                resources[table] = 1

        return resources

    @property
    def priority(self):
        """
        Dynamically set task priority
        """

        # deprioritise large values of K, m or n
        offset = sum([getattr(self, name) for name in self.get_param_names() if name in ['k', 'm', 'n']])

        # prioritise chromosomes by number, as running the largest chroms first is more efficient for multithreading
        if hasattr(self, 'chrom') and hasattr(self, 'species'):
            offset = CHROMOSOMES[self.assembly].index(self.chrom) + 1

        return 1000 - offset if offset else 0

    @property
    def classname(self):
        """
        The name of the current class
        """
        return type(self).__name__

    @property
    def basename(self):
        """
        Collapse all the param values into a hyphen delimited list
        """
        params = []

        for name, value in self._all_params():
            if name == 'chrom':
                params.append('chr{}'.format(value))
            elif isinstance(value, str):
                params.append(value)
            elif isinstance(value, bool):
                if value:
                    params.append(name)
            elif isinstance(value, Iterable):
                params.append('{}({})'.format(name, u','.join([a for a in value])))
            else:
                if value is not None:
                    params.append('{}{}'.format(name, value))

        return '-'.join(params)

    @property
    def binomial(self):
        """
        Scientific binomial name of the species
        """
        return BINOMIAL_NAME[self.species]

    @property
    def assembly(self):
        """
        Identifier of the reference assembly for the species
        """
        return REF_ASSEMBLY[self.species]

    @property
    def chromosomes(self):
        """
        List of chromosome names (e.g. 1, 2, ..., X, Y, MT)
        """
        return CHROMOSOMES[self.assembly]

    @property
    def autosomes(self):
        """
        List of autosomal chromosome names (e.g. 1, 2, ..., 29)
        """
        return [chrom for chrom in self.chromosomes if chrom.isdigit()]

    @property
    def java_mem(self):
        """
        Memory to allocate to java processes
        """
        return "-Xmx{}G".format(self.resources['ram-gb'])

    def _all_params(self):
        """
        Get all the params as a (name, value) tuple
        """
        return [(name, getattr(self, name)) for name in self.get_param_names()]


class DatabaseTask(PipelineTask):
    """
    Pipeline task with a db connection to retrieve sample metadata

    :type species: str
    """
    species = luigi.Parameter()

    _dbc = None
    _outgroup = None

    @property
    def dbc(self):
        """
        Create a private connection to the db
        """
        if self._dbc is None:
            self._dbc = Database(self.species)

        return self._dbc

    @property
    def outgroup(self):
        """
        Name of the outgroup sample
        """
        if self._outgroup is None:
            self._outgroup = self.dbc.get_record('samples', {'populaion': OUTGROUP_POP})

        return self._outgroup['name']

    def list_samples(self, ancient=None, modern=None, outgroup=False):
        """
        List all the samples.

        :return: {(pop, sample): record, ...}
        """
        conds = {}
        if ancient or modern:
            conds['ancient'] = 1 if ancient else 0

        samples = self.dbc.get_records('samples', conds, sort='population, name', key=None)

        return OrderedDict([((sample['population'], sample['name']), sample) for sample in samples
                            if outgroup or sample['population'] != OUTGROUP_POP])

    def list_populations(self, ancient=None, modern=None, outgroup=False):
        """
        List all the populations, and their samples.

        :return: {pop: {sample: record, ...}, ...}
        """
        samples = self.list_samples(ancient, modern, outgroup)

        # group the samples by population
        data = defaultdict(dict)

        for pop, sample in samples:
            data[pop][sample] = samples[(pop, sample)]

        return data

    def list_accessions(self, sample):
        """
        List all the accession codes for the sample.

        :return: [accession, ...]
        """
        return self.dbc.get_records_sql("""
            SELECT sr.accession
              FROM samples s 
              JOIN sample_runs sr
                ON sr.sample_id = s.id
             WHERE s.name = '{sample}'
               """.format(sample=sample), key='accession').keys()


class PipelineExternalTask(luigi.ExternalTask, PipelineTask):
    """
    Let ExternalTasks access dynamic properties of the PipelineTask
    """
    pass


class PipelineWrapperTask(luigi.WrapperTask, DatabaseTask):
    """
    Let WrapperTasks access dynamic properties of the PipelineTask
    """
    pass
