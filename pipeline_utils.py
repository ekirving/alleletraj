#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess

from db_conn import db_conn
from pipeline_consts import *


def run_cmd(cmd, shell=False, background=False):
    """
    Executes the given command in a system subprocess

    :param cmd: The system command to run (list|string)
    :param shell: Use the native shell
    :return: The stdout stream
    """
    # subprocess only accepts strings
    cmd = [str(args) for args in cmd]

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
    saved = list(ranges[0])

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
