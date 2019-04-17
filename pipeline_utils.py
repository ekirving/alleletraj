#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess
import luigi
import os

# import common libraries
from collections import OrderedDict, defaultdict, Iterable
from multiprocessing import Process

# import my libraries
from db_conn import db_conn
from pipeline_consts import *


def run_cmd(cmd, shell=False, background=False, stdout=None, stderr=None):
    """
    Executes the given command in a system subprocess

    :param cmd: The system command to run (list|string)
    :param shell: Use the native shell
    :param background: Flag to tun the process in the background
    :param stdout: File handle to redirect stdout
    :param stderr: File handle to redirect stderr
    :param pwd: Handle commands than run natively from other locations
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


def extract_qtl_fields(dbfile, fields):
    data = defaultdict(list)

    # get all the QTL IDs for this scope
    with open(dbfile, 'rU') as fin:

        # get the column headers
        header = fin.readline().split('\t')

        # get the index of the target fields
        columns = [(idx, field) for idx, field in enumerate(header) if field in fields]

        for line in fin:
            try:
                line = line.split('\t')
                for idx, field in columns:
                    datum = line[idx].strip()
                    if datum:
                        data[field].append(datum)
            except IndexError:
                # ignore badly formatted lines
                pass

    return data


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


def unicode_truncate(s, length, encoding='utf-8'):
    encoded = s.encode(encoding)[:length]
    return encoded.decode(encoding, 'ignore')


def dump(obj):
    for attr in dir(obj):
        if hasattr(obj, attr):
            print("obj.%s = %s" % (attr, getattr(obj, attr)))


class PipelineTask(luigi.Task):
    """
    PrioritisedTask that implements a dynamic priority method
    """
    priority = PRIORITY_LOW
    resources = {'cpu-cores': 1}

    @property
    def priority(self):
        """
        Set a dynamic priority for tasks.
        """
        # TODO do CPU intensive tasks first

        # deprioritise large values of K or m
        offset = sum([getattr(self, name) for name in self.get_param_names() if name in ['k', 'm']])

        # prioritise chromosomes by size, as running the largest chroms first is more efficient for multithreading
        if hasattr(self, 'chrom') and hasattr(self, 'species'):
            chrom = self.chrom.replace('chr', '')
            sizes = CHROM_SIZE[self.species]
            offset = sorted(sizes.values(), reverse=True).index(sizes[chrom]) + 1

        return 100-offset if offset else 0

    @property
    def basename(self):
        """
        Collapse all the param values into a hyphen delimited list
        """
        params = []

        for name, value in self.all_params():
            if isinstance(value, str):
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

    def all_params(self):
        """
        Get all the params as a (name, value) tuple
        :return:
        """
        return [(name, getattr(self, name)) for name in self.get_param_names()]

    @property
    def java_mem(self):
        # memory to allocate to java
        return "-Xmx{}G".format(self.resources['ram-gb'])