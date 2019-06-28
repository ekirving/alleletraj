#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
import os
import subprocess
from collections import Iterable, OrderedDict, defaultdict
from datetime import timedelta
from time import time

# third party modules
import luigi

# import my libraries
from alleletraj.const import CHROMOSOMES, REF_ASSEMBLY, BINOMIAL_NAME, OUTGROUP_POP
from alleletraj.db.conn import Database

# enforce max interval size of 1 Gb
MAX_INTERVAL_SIZE = int(1e6)


def run_cmd(cmd, shell=False, stdout=None, stderr=None, verbose=True):
    """
    Executes the given command in a system subprocess.

    :param cmd: The system command to run (list|string)
    :param shell: Use the native shell
    :param stdout: File handle to redirect stdout
    :param stderr: File handle to redirect stderr
    :param verbose: Print the command before running it
    :return: The stdout stream
    """
    # subprocess only accepts strings
    cmd = [str(args) for args in cmd]

    if verbose:
        print(' '.join(cmd))

    stdout = stdout or subprocess.PIPE
    stderr = stderr or subprocess.PIPE

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
    return out or err


def merge_intervals(ranges, capped=True):
    """
    Merge overlapping intervals, so we only check each site once.
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
    Downloads a remote url to a local file path using cURL.
    """
    run_cmd(['curl', '-s', '--output', filename, url])


def get_chrom_sizes(fai_file, exclude_scaffolds=True):
    """
    Get the names and sizes of all the chromosomes in a reference by iterating over the index.
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
    Pipeline task that implements several dynamic attributes.
    """
    resources = {'cpu-cores': 1}

    # TODO consider recursing through the whole dependency tree
    # https://stackoverflow.com/questions/42465710/can-python-task-scheduler-luigi-detect-indirect-dependencies

    @property
    def priority(self):
        """
        Dynamically set task priority.

        NOTE luigi task inherit the priority from their children
        """

        # deprioritise large values of K, m or n
        offset = sum([getattr(self, name) for name in self.get_param_names() if name in ['k', 'm', 'n']])

        try:
            # prioritise chromosomes by number, as running the largest chroms first is more efficient for multithreading
            offset = CHROMOSOMES[self.assembly].index(self.chrom) + 1
        except AttributeError:
            pass

        return 10000 - offset if offset else 0

    @property
    def basename(self):
        """
        Collapse all the param values into a hyphen delimited list.
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
        Scientific binomial name of the species.
        """
        return BINOMIAL_NAME[self.species]

    @property
    def assembly(self):
        """
        Identifier of the reference assembly for the species.
        """
        return REF_ASSEMBLY[self.species]

    @property
    def chromosomes(self):
        """
        List of chromosome names (e.g. 1, 2, ..., X, Y, MT).
        """
        return CHROMOSOMES[self.assembly]

    @property
    def autosomes(self):
        """
        List of autosomal chromosome names (e.g. 1, 2, ..., 29).
        """
        return [chrom for chrom in self.chromosomes if chrom.isdigit()]

    @property
    def java_mem(self):
        """
        Memory to allocate to java processes.
        """
        return "-Xmx{}G".format(self.resources['ram-gb'])

    @property
    def java_gc_threads(self):
        """
        Threads to allocate to garbage collection in java processes.

        See 'Why does a Picard tool use so many threads?' in https://broadinstitute.github.io/picard/faq.html
        """
        return "-XX:ParallelGCThreads={}".format(self.resources['cpu-cores'])

    def _all_params(self):
        """
        Get all the params as a (name, value) tuple.
        """
        return [(name, getattr(self, name)) for name in self.get_param_names()]

    def all_params(self):
        """
        Get all the params as a dictionary, so we can use the **kwargs syntax.
        """
        return dict(self._all_params())

    def input_targets(self, ext=None):
        """
        Get all input targets, filtered by extension.
        """
        inputs = self.input()

        targets = []

        if isinstance(inputs, luigi.LocalTarget):
            targets.append(inputs)
        else:
            for task in inputs:
                if isinstance(task, luigi.LocalTarget):
                    targets.append(task)
                else:
                    for target in task:
                        targets.append(target)

        return [target for target in targets if target.path.endswith(ext)] if ext else targets


class DatabaseTask(PipelineTask):
    """
    Pipeline task with a database connection to retrieve sample metadata.

    :type species: str
    """
    species = luigi.Parameter()

    _outgroup = None
    _sample = None
    _accession = None

    @property
    def resources(self):
        """
        Handle table locking in the luigi scheduler so we don't crash the database by flooding it with locked queries.
        """
        resources = {'cpu-cores': 1}

        if hasattr(self, 'db_lock_tables'):
            for table in self.db_lock_tables:
                # resolve table names that contain task parameters (e.g. chrom)
                table = table.format(**self.param_kwargs)

                # table resources are species specific
                table = '{}_{}'.format(self.species, table)

                # the default quantity for resources is '1', so this locks the table for other tasks
                resources[table] = 1

        return resources

    @property
    def dbc(self):
        """
        Create a single use connection to the database.

        Connection is non-persistant becauase there can be thousands of tasks instantiated during scheduling.
        """
        return Database(self.species)

    @property
    def outgroup(self):
        """
        Name of the outgroup sample.
        """
        if self._outgroup is None:
            self._outgroup = self.dbc.get_record('samples', {'population': OUTGROUP_POP})

        return self._outgroup['name']

    @property
    def sample_data(self):
        """
        The sample record.
        """
        if self._sample is None:
            self._sample = self.dbc.get_record('samples', {'name': self.sample})

        return self._sample

    @property
    def accession_data(self):
        """
        The accession record.
        """
        if self._accession is None:
            self._accession = self.dbc.get_record('sample_runs', {'accession': self.accession})
            self._accession['paired'] = bool(int(self._accession['paired']))  # cast as boolean

        return self._accession

    def list_samples(self, ancient=False, modern=False, outgroup=False):
        """
        List of all the population and sample tuples.

        :return: {(pop, sample): record, ...}
        """
        conds = {}
        if ancient ^ modern:  # xor
            conds['ancient'] = 1 if ancient else 0

        samples = self.dbc.get_records('samples', conds, sort='population, name', key=None)

        return OrderedDict([((sample['population'], sample['name']), sample)
                            for sample in samples if outgroup or sample['population'] != OUTGROUP_POP])

    def list_populations(self, ancient=False, modern=False, outgroup=False):
        """
        Nested list of all the populations and their samples.

        :return: {pop: {sample: record, ...}, ...}
        """
        samples = self.list_samples(ancient, modern, outgroup)

        # group the samples by population
        data = defaultdict(dict)

        for pop, sample in samples:
            data[pop][sample] = samples[(pop, sample)]

        return data

    def list_accessions(self):
        """
        List of all the accession codes for the current sample.

        :return: [accession, ...]
        """
        accessions = self.dbc.get_records_sql("""
            SELECT sr.accession
              FROM samples s 
              JOIN sample_runs sr
                ON sr.sample_id = s.id
             WHERE s.name = '{sample}'
               """.format(sample=self.sample), key=None)

        return [row['accession'] for row in accessions]


class MySQLTask(DatabaseTask):
    """
    Pipeline task which only acts on the database.

    Automatically handles completion tracking.
    """
    _dbc = None

    _variables = {
        # reduce spurious locking errors when running lots of updates on the same table
        'innodb_lock_wait_timeout': 600
    }

    @property
    def dbc(self):
        """
        Create a persistent connection to the database.
        """
        if self._dbc is None:
            self._dbc = Database(self.species, self._variables)

        return self._dbc

    def output(self):
        return luigi.LocalTarget(
            'data/db/{}/{}/{}-{}.log'.format(Database.get_name(self.species), self.task_module, type(self).__name__,
                                             self.basename))

    def queries(self):
        """
        To be overridden in a subclass.
        """
        pass

    def run(self):
        """
        Wrapper for :ref:`MySQLTask.queries` to handle completion tracking for SQL queries.
        """
        start = time()

        # run all the sql queries
        self.queries()

        with self.output().open('w') as fout:
            fout.write('Execution took {}'.format(timedelta(seconds=time() - start)))


class PipelineExternalTask(luigi.ExternalTask, PipelineTask):
    """
    Let luigi.ExternalTasks share the properties of the PipelineTask.
    """
    pass


class PipelineWrapperTask(luigi.WrapperTask, DatabaseTask):
    """
    Let luigi.WrapperTasks share the properties of the DatabaseTask.
    """
    pass
