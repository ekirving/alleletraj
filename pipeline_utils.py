#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess

# import common libraries
from collections import OrderedDict, defaultdict
from datetime import timedelta
from pprint import pprint
from time import time
from multiprocessing import Process

# import my libraries
from db_conn import db_conn
from pipeline_consts import *


def run_cmd(cmd, shell=False, background=False):
    """
    Executes the given command in a system subprocess

    :param cmd: The system command to run (list|string)
    :param shell: Use the native shell
    :param background: Run the process in the background
    :return: The stdout stream
    """
    # subprocess only accepts strings
    cmd = [str(args) for args in cmd]

    # print(u' '.join(cmd))

    if background:
        subprocess.Popen(cmd, shell=shell)

    else:
        # run the command
        proc = subprocess.Popen(cmd, shell=shell, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # fetch the output and error
        (stdout, stderr) = proc.communicate()

        # bail if something went wrong
        if proc.returncode:
            raise Exception(stderr)

        return stdout


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
