# -*- coding: utf-8 -*-

import _io
import sys
import os
import subprocess as sp

from snakemake.utils import format


__author__ = "Johannes Köster"


STDOUT = sys.stdout
if not isinstance(sys.stdout, _io.TextIOWrapper):
    # workaround for nosetest since it overwrites sys.stdout
    # in a strange way that does not work with Popen
    STDOUT = None


class shell:
    _process_args = {}
    _process_prefix = ""

    @classmethod
    def executable(cls, cmd):
        cls._process_args["executable"] = cmd

    @classmethod
    def prefix(cls, prefix):
        cls._process_prefix = format(prefix, stepout=2)

    def __new__(
        cls, cmd, *args, async=False, iterable=False, read=False, **kwargs):
        cmd = format(cmd, *args, stepout=2, **kwargs)
        stdout = sp.PIPE if iterable or async or read else STDOUT
        proc = sp.Popen(cls._process_prefix + cmd, shell=True, stdout=stdout,
            close_fds=True, **cls._process_args)

        if iterable:
            return proc.stdout
        if read:
            return proc.stdout.read()
        if async:
            return proc
        proc.wait()
        return None


if "SHELL" in os.environ:
    shell.executable(os.environ["SHELL"])
