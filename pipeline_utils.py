#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess
import luigi
import os

# import common libraries
from collections import Iterable
from multiprocessing import Process

# import my libraries
from pipeline_consts import CHROM_SIZE, CPU_CORES_ONE, REF_ASSEMBLY, OUTGROUP, BINOMIAL_NAME, SAMPLES

from database import Database

# enforce max interval size of 1 Gb
MAX_INTERVAL_SIZE = int(1e6)


def run_cmd(cmd, shell=False, background=False, stdout=None, stderr=None):
    """
    Executes the given command in a system subprocess

    :param cmd: The system command to run (list|string)
    :param shell: Use the native shell
    :param background: Flag to tun the process in the background
    :param stdout: File handle to redirect stdout
    :param stderr: File handle to redirect stderr
    :return: The stdout stream
    """
    # subprocess only accepts strings
    cmd = [str(args) for args in cmd]

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

        return out


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


def run_in_parallel(*fns):
    """
    Simple wrapper for running functions in parallel.
    """
    proc = []

    for fn in fns:
        p = Process(target=fn)
        p.start()
        proc.append(p)

    for p in proc:
        p.join()


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


class PipelineTask(luigi.Task):
    """
    PrioritisedTask that implements a several dynamic attributes
    """
    @property
    def resources(self):
        """
        Dynamically set task resource usage.
        """
        resources = {'cpu-cores': CPU_CORES_ONE}

        if hasattr(self, 'db_lock_tables'):
            for table in self.db_lock_tables:
                resources[table] = 1

        return resources

    @property
    def priority(self):
        """
        Dynamically set task priority.
        """
        # TODO do CPU intensive tasks first

        # deprioritise large values of K or m
        offset = sum([getattr(self, name) for name in self.get_param_names() if name in ['k', 'm']])

        # prioritise chromosomes by size, as running the largest chroms first is more efficient for multithreading
        if hasattr(self, 'chrom') and hasattr(self, 'species'):
            chrom = self.chrom.replace('chr', '')
            sizes = CHROM_SIZE[self.assembly]
            offset = sorted(sizes.values(), reverse=True).index(sizes[chrom]) + 1

        return 100 - offset if offset else 0

    @property
    def basename(self):
        """
        Collapse all the param values into a hyphen delimited list
        """
        params = []

        for name, value in self.all_params():
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
    def assembly(self):
        """
        Identifier of the reference assembly for the species
        """
        return REF_ASSEMBLY[self.species]

    @property
    def binomial(self):
        """
        Scientific binomial name of the species
        """
        return BINOMIAL_NAME[self.species]

    @property
    def chromosomes(self):
        """
        List of chromosomes identifiers (e.g. 1, 2, ..., X, Y)
        """
        return CHROM_SIZE[self.assembly]

    @property
    def classname(self):
        """
        The name of the current class
        """
        return type(self).__name__

    @property
    def java_mem(self):
        """
        Memory to allocate to java processes
        """
        return "-Xmx{}G".format(self.resources['ram-gb'])

    @property
    def outgroup(self):
        """
        Identifier of the outgroup sample
        """
        return OUTGROUP[self.species]

    @property
    def populations(self):
        """
        List of the populations for the species
        """
        return SAMPLES[self.species]

    @property
    def samples(self):
        """
        List of the modern samples for the species
        """
        return SAMPLES[self.species][self.population]

    def all_params(self):
        """
        Get all the params as a (name, value) tuple
        """
        return [(name, getattr(self, name)) for name in self.get_param_names()]

    def db_conn(self):
        """
        Create a private connection to the database
        """
        return Database(self.species)


class PipelineExternalTask(luigi.ExternalTask, PipelineTask):
    """
    Let ExternalTasks access dynamic properties of the PipelineTask
    """
    pass


class PipelineWrapperTask(luigi.WrapperTask, PipelineTask):
    """
    Let WrapperTasks access dynamic properties of the PipelineTask
    """
    pass
